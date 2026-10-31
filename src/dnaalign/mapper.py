"""Seed-and-extend read mapping.

For each read:

1. **Seed** - cut the read into k-mers and look each one up in the
   :class:`~dnaalign.kmer_index.KmerIndex` (fast hash lookups).
2. **Vote** - every hit says "the read starts near reference position
   ``ref_pos - read_pos``" (its *diagonal*). Hits from the true location agree,
   random hits do not, so we count votes per diagonal.
3. **Extend** - for the best few diagonals, run a *banded* Smith-Waterman
   alignment of the read against a small reference window (a strip of width
   ``2 * band + 1`` around the diagonal), then keep the best score.

Both strands are tried: the read itself and its reverse complement.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from dnaalign.align import local_align
from dnaalign.kmer_index import KmerIndex, reverse_complement
from dnaalign.scoring import Scoring, dna_scoring


@dataclass
class MappingResult:
    """Where a read landed. ``position`` is 0-based on the forward strand."""

    name: str
    mapped: bool
    strand: str = "*"
    position: int = -1
    score: float = 0.0
    cigar: str = "*"
    seeds: int = 0
    second_score: float = 0.0

    @property
    def mapq(self) -> int:
        """A rough confidence (0-60) from the gap between best and second-best score."""
        if not self.mapped:
            return 0
        if self.second_score <= 0:
            return 60
        return max(0, min(60, int(60 * (self.score - self.second_score) / self.score)))


class ReadMapper:
    """Map short reads to one reference sequence.

    Args:
        reference: the genome to map against.
        k: k-mer (seed) length.
        band: half-width of the band used in the extension alignment.
        max_candidates: how many top diagonals to extend per read.
        min_score_fraction: a read is "mapped" only if its score is at least this
            fraction of the perfect score (``len(read) * match``).
        scoring: alignment scoring (default :func:`dna_scoring`).
    """

    def __init__(
        self,
        reference: str,
        k: int = 15,
        band: int = 16,
        max_candidates: int = 3,
        min_score_fraction: float = 0.3,
        scoring: Scoring | None = None,
    ) -> None:
        self.reference = reference.upper()
        self.index = KmerIndex(self.reference, k)
        self.band = band
        self.max_candidates = max_candidates
        self.min_score_fraction = min_score_fraction
        self.scoring = scoring or dna_scoring()

    def _candidates(self, read: str) -> list[tuple[int, int]]:
        """Vote for diagonals; return ``(votes, diagonal)`` best first."""
        votes: Counter[int] = Counter()
        for read_pos, ref_pos in self.index.seeds(read):
            votes[ref_pos - read_pos] += 1
        # Merge nearby diagonals (indels shift the diagonal a little).
        merged: list[tuple[int, int]] = []
        for diag, _count in sorted(votes.items(), key=lambda kv: -kv[1]):
            if all(abs(diag - d) > self.band for _, d in merged):
                total = sum(c for d2, c in votes.items() if abs(d2 - diag) <= self.band)
                merged.append((total, diag))
        merged.sort(key=lambda t: -t[0])
        return merged[: self.max_candidates]

    def _extend(self, read: str, diag: int) -> tuple[float, int, str]:
        """Banded local alignment around ``diag``; returns (score, ref_start, cigar)."""
        w0 = max(0, diag - self.band)
        window = self.reference[w0 : diag + len(read) + self.band]
        aln = local_align(read, window, self.scoring, band=self.band, band_offset=diag - w0)
        # Position of the read's first base, counting soft-clipped read bases.
        return aln.score, w0 + aln.start_b - aln.start_a, aln.cigar()

    def map_read(self, name: str, read: str) -> MappingResult:
        """Map one read; tries both strands and keeps the best alignment."""
        read = read.upper()
        tried: list[tuple[float, str, int, str, int]] = []
        for strand, seq in (("+", read), ("-", reverse_complement(read))):
            for votes, diag in self._candidates(seq):
                score, pos, cigar = self._extend(seq, diag)
                tried.append((score, strand, pos, cigar, votes))
        if not tried:
            return MappingResult(name, mapped=False)
        tried.sort(key=lambda t: -t[0])
        score, strand, pos, cigar, votes = tried[0]
        second = tried[1][0] if len(tried) > 1 else 0.0
        if score < self.min_score_fraction * len(read) * self.scoring.match:
            return MappingResult(name, mapped=False, seeds=votes)
        return MappingResult(name, True, strand, pos, score, cigar, votes, second)

    def map_reads(self, reads: list[tuple[str, str]]) -> list[MappingResult]:
        """Map a list of ``(name, sequence)`` pairs."""
        return [self.map_read(name, seq) for name, seq in reads]
