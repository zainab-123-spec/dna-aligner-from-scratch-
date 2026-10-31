"""Command line interface: ``dnaalign align``, ``dnaalign map``, ``dnaalign simulate``."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from dnaalign.align import global_align, local_align
from dnaalign.mapper import ReadMapper
from dnaalign.scoring import Scoring, blosum62
from dnaalign.seqio import Record, read_fasta, read_sequences, write_fastq
from dnaalign.simulate import simulate_reads


def _sequence_arg(value: str) -> str:
    """Accept either a literal sequence or a path to a FASTA/FASTQ file."""
    path = Path(value)
    if path.is_file():
        return next(iter(read_sequences(path))).sequence
    return value.upper()


def cmd_align(args: argparse.Namespace) -> int:
    a, b = _sequence_arg(args.a), _sequence_arg(args.b)
    if args.matrix == "blosum62":
        scoring = blosum62(args.gap_open, args.gap_extend)
    else:
        scoring = Scoring(args.match, args.mismatch, args.gap_open, args.gap_extend)
    aligner = local_align if args.mode == "local" else global_align
    aln = aligner(a, b, scoring)
    print(aln)
    print(f"cigar = {aln.cigar()}  identity = {aln.identity:.1%}")
    return 0


def cmd_map(args: argparse.Namespace) -> int:
    reference = next(iter(read_fasta(args.ref))).sequence
    t0 = time.perf_counter()
    mapper = ReadMapper(reference, k=args.k, band=args.band)
    t1 = time.perf_counter()
    out = open(args.out, "w") if args.out else sys.stdout
    out.write("read\tmapped\tstrand\tpos\tscore\tmapq\tcigar\n")
    n = n_mapped = 0
    for rec in read_sequences(args.reads):
        r = mapper.map_read(rec.name, rec.sequence)
        n += 1
        n_mapped += r.mapped
        pos = r.position + 1 if r.mapped else 0  # 1-based like SAM
        out.write(f"{r.name}\t{int(r.mapped)}\t{r.strand}\t{pos}\t{r.score:g}\t{r.mapq}\t{r.cigar}\n")
    t2 = time.perf_counter()
    if out is not sys.stdout:
        out.close()
    print(
        f"indexed {len(reference):,} bp in {t1 - t0:.2f}s; mapped {n_mapped}/{n} reads in {t2 - t1:.2f}s",
        file=sys.stderr,
    )
    return 0


def cmd_simulate(args: argparse.Namespace) -> int:
    reference = next(iter(read_fasta(args.ref))).sequence
    reads = simulate_reads(reference, args.n, args.length, args.error_rate, seed=args.seed)
    out = open(args.out, "w") if args.out else sys.stdout
    write_fastq((Record(r.name, r.sequence) for r in reads), out)
    if out is not sys.stdout:
        out.close()
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="dnaalign", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("align", help="align two sequences (strings or FASTA files)")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--mode", choices=["global", "local"], default="global")
    p.add_argument("--match", type=float, default=1)
    p.add_argument("--mismatch", type=float, default=-1)
    p.add_argument("--gap-open", type=float, default=-2)
    p.add_argument("--gap-extend", type=float, default=-2)
    p.add_argument("--matrix", choices=["none", "blosum62"], default="none")
    p.set_defaults(func=cmd_align)

    p = sub.add_parser("map", help="map reads (FASTA/FASTQ) to a reference FASTA")
    p.add_argument("--ref", required=True)
    p.add_argument("--reads", required=True)
    p.add_argument("-k", type=int, default=15)
    p.add_argument("--band", type=int, default=16)
    p.add_argument("--out", help="output TSV (default: stdout)")
    p.set_defaults(func=cmd_map)

    p = sub.add_parser("simulate", help="simulate reads with errors from a reference")
    p.add_argument("--ref", required=True)
    p.add_argument("-n", type=int, default=1000)
    p.add_argument("--length", type=int, default=150)
    p.add_argument("--error-rate", type=float, default=0.01)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", help="output FASTQ (default: stdout)")
    p.set_defaults(func=cmd_simulate)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
