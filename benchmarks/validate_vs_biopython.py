"""Compare dnaalign scores with Biopython's PairwiseAligner on random cases.

Biopython is used ONLY here (and in tests) as an independent reference.
Run: python benchmarks/validate_vs_biopython.py --cases 2000
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from Bio.Align import PairwiseAligner, substitution_matrices

from dnaalign import Scoring, blosum62, global_align, local_align, score_alignment

PROTEIN = "ARNDCQEGHILKMFPSTWYV"


def biopython_score(a: str, b: str, sc: Scoring, mode: str, protein: bool) -> float:
    al = PairwiseAligner()
    al.mode = mode
    if protein:
        al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    else:
        al.match_score, al.mismatch_score = sc.match, sc.mismatch
    al.open_gap_score, al.extend_gap_score = sc.gap_open, sc.gap_extend
    return al.score(a, b)


def random_case(rng: random.Random) -> tuple[str, str, Scoring, str, bool]:
    protein = rng.random() < 0.25
    alphabet = PROTEIN if protein else "ACGT"
    a = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 80)))
    # make b a mutated copy of a part of a half the time, random otherwise
    if rng.random() < 0.5:
        b = "".join(c if rng.random() > 0.2 else rng.choice(alphabet) for c in a)
        b = b[rng.randint(0, len(b) // 3) :]
        b = b or rng.choice(alphabet)
    else:
        b = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 80)))
    if protein:
        sc = blosum62(-rng.randint(5, 12), -rng.randint(1, 2))
    else:
        ext = -rng.randint(1, 3)
        op = ext if rng.random() < 0.5 else ext - rng.randint(1, 5)
        sc = Scoring(rng.randint(1, 3), -rng.randint(1, 4), op, ext)
    return a, b, sc, rng.choice(["global", "local"]), protein


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--cases", type=int, default=3000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", default="benchmarks/results/validation.json")
    args = p.parse_args()
    rng = random.Random(args.seed)
    counts: dict[str, list[int]] = {}
    mismatches = []
    for _ in range(args.cases):
        a, b, sc, mode, protein = random_case(rng)
        aln = (global_align if mode == "global" else local_align)(a, b, sc)
        ref = biopython_score(a, b, sc, mode, protein)
        rescored = score_alignment(aln.aligned_a, aln.aligned_b, sc) if aln.aligned_a else 0.0
        ok = abs(aln.score - ref) < 1e-9 and abs(rescored - aln.score) < 1e-9
        key = f"{mode}/{'protein-BLOSUM62' if protein else 'DNA'}/{'linear' if sc.is_linear else 'affine'}"
        counts.setdefault(key, [0, 0])
        counts[key][0] += ok
        counts[key][1] += 1
        if not ok:
            mismatches.append({"a": a, "b": b, "mode": mode, "ours": aln.score, "biopython": ref})
    total_ok = sum(v[0] for v in counts.values())
    summary = {
        "cases": args.cases,
        "agree": total_ok,
        "by_category": {k: {"agree": v[0], "total": v[1]} for k, v in sorted(counts.items())},
        "mismatch_examples": mismatches[:5],
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
