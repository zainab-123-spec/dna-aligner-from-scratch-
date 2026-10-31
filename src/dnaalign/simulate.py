"""Simulate short sequencing reads with errors, so we know the right answer.

Each read records its true origin in its name, e.g.
``read7_pos1234_strand-`` means it came from reference position 1234
(0-based, forward strand) and was reverse-complemented.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

from dnaalign.kmer_index import reverse_complement

BASES = "ACGT"


@dataclass
class SimulatedRead:
    name: str
    sequence: str
    true_position: int
    strand: str
    n_errors: int


def mutate(seq: str, sub_rate: float, indel_rate: float, rng: random.Random) -> tuple[str, int]:
    """Add random substitutions, insertions and deletions to ``seq``.

    ``indel_rate`` is split evenly between insertions and deletions.
    Returns the mutated sequence and the number of errors added.
    """
    out, errors = [], 0
    for base in seq:
        r = rng.random()
        if r < sub_rate:
            out.append(rng.choice([b for b in BASES if b != base]))
            errors += 1
        elif r < sub_rate + indel_rate / 2:
            errors += 1  # deletion: drop this base
        elif r < sub_rate + indel_rate:
            out.append(base)
            out.append(rng.choice(BASES))  # insertion after this base
            errors += 1
        else:
            out.append(base)
    return "".join(out), errors


def simulate_reads(
    reference: str,
    n_reads: int,
    length: int = 150,
    error_rate: float = 0.01,
    indel_fraction: float = 0.1,
    seed: int = 0,
) -> list[SimulatedRead]:
    """Draw ``n_reads`` reads from random places/strands of ``reference``.

    ``error_rate`` is the total per-base error rate; ``indel_fraction`` of it
    is insertions/deletions and the rest substitutions.
    """
    if length > len(reference):
        raise ValueError("read length longer than reference")
    rng = random.Random(seed)
    reads = []
    for i in range(n_reads):
        pos = rng.randrange(0, len(reference) - length + 1)
        fragment = reference[pos : pos + length]
        strand = rng.choice("+-")
        sub = error_rate * (1 - indel_fraction)
        seq, n_err = mutate(fragment, sub, error_rate * indel_fraction, rng)
        if strand == "-":
            seq = reverse_complement(seq)
        reads.append(SimulatedRead(f"read{i}_pos{pos}_strand{strand}", seq, pos, strand, n_err))
    return reads
