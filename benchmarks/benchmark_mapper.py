"""Benchmark the read mapper on the real SARS-CoV-2 reference (NC_045512.2).

Simulates reads at several error rates, maps them, and checks each read
against its known true position (correct = right strand and within
``--tolerance`` bp). Also times the O(n*m) global alignment for a figure.

Run: python benchmarks/benchmark_mapper.py
"""

from __future__ import annotations

import argparse
import json
import platform
import random
import time
from pathlib import Path

from dnaalign import ReadMapper, dna_scoring, global_align
from dnaalign.seqio import read_fasta
from dnaalign.simulate import simulate_reads

HERE = Path(__file__).parent


def bench_mapper(reference: str, args: argparse.Namespace) -> list[dict]:
    t0 = time.perf_counter()
    mapper = ReadMapper(reference, k=args.k, band=args.band)
    index_s = time.perf_counter() - t0
    print(f"index: {len(mapper.index):,} distinct {args.k}-mers in {index_s:.3f}s")
    rows = []
    for rate in args.error_rates:
        reads = simulate_reads(reference, args.reads, args.length, rate, seed=args.seed)
        t0 = time.perf_counter()
        results = [mapper.map_read(r.name, r.sequence) for r in reads]
        secs = time.perf_counter() - t0
        mapped = sum(r.mapped for r in results)
        correct = sum(
            res.mapped
            and res.strand == truth.strand
            and abs(res.position - truth.true_position) <= args.tolerance
            for res, truth in zip(results, reads)
        )
        row = {
            "error_rate": rate,
            "reads": len(reads),
            "mapped_pct": round(100 * mapped / len(reads), 2),
            "correct_pct": round(100 * correct / len(reads), 2),
            "precision_pct": round(100 * correct / mapped, 2) if mapped else 0.0,
            "seconds": round(secs, 2),
            "reads_per_sec": round(len(reads) / secs, 1),
            "index_seconds": round(index_s, 3),
        }
        print(row)
        rows.append(row)
    return rows


def bench_dp(sizes: list[int]) -> list[dict]:
    rng = random.Random(1)
    sc = dna_scoring()
    rows = []
    for n in sizes:
        a = "".join(rng.choice("ACGT") for _ in range(n))
        b = "".join(rng.choice("ACGT") for _ in range(n))
        t0 = time.perf_counter()
        global_align(a, b, sc)
        secs = time.perf_counter() - t0
        rows.append({"n": n, "cells": (n + 1) ** 2, "seconds": round(secs, 4)})
        print(rows[-1])
    return rows


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--ref", default=str(HERE / "data" / "NC_045512.2.fasta"))
    p.add_argument("--reads", type=int, default=2000)
    p.add_argument("--length", type=int, default=150)
    p.add_argument("-k", type=int, default=15)
    p.add_argument("--band", type=int, default=16)
    p.add_argument("--tolerance", type=int, default=10)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--error-rates", type=float, nargs="+", default=[0.0, 0.01, 0.02, 0.05, 0.10])
    p.add_argument("--dp-sizes", type=int, nargs="+", default=[100, 200, 400, 800, 1600])
    p.add_argument("--out", default=str(HERE / "results" / "mapper_benchmark.json"))
    args = p.parse_args()

    rec = next(iter(read_fasta(args.ref)))
    print(f"reference {rec.name}: {len(rec.sequence):,} bp")
    result = {
        "reference": rec.name,
        "reference_bp": len(rec.sequence),
        "read_length": args.length,
        "k": args.k,
        "band": args.band,
        "tolerance_bp": args.tolerance,
        "python": platform.python_version(),
        "machine": platform.machine(),
        "mapper": bench_mapper(rec.sequence, args),
        "dp_scaling": bench_dp(args.dp_sizes),
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
