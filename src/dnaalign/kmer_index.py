"""A k-mer index: a hash map from every k-letter word to where it occurs.

Example with k = 3 and reference ``ACGTACG``::

    ACG -> [0, 4]
    CGT -> [1]
    GTA -> [2]
    TAC -> [3]

Building it is O(L) for a reference of length L (one dictionary insert per
position). Looking a k-mer up is O(1) on average, which is why read mappers
use it to find *where to look* before running the slow O(n*m) alignment.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterator

DNA = set("ACGT")
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    """Return the other DNA strand read 5'->3' (A<->T, C<->G, reversed)."""
    return seq.translate(_COMPLEMENT)[::-1]


def iter_kmers(seq: str, k: int) -> Iterator[tuple[int, str]]:
    """Yield ``(position, kmer)`` for every k-mer made only of A/C/G/T."""
    if k <= 0:
        raise ValueError("k must be positive")
    for pos in range(len(seq) - k + 1):
        kmer = seq[pos : pos + k]
        if set(kmer) <= DNA:
            yield pos, kmer


class KmerIndex:
    """Hash map ``kmer -> list of positions`` over a reference sequence.

    Args:
        reference: the genome (upper-case A/C/G/T; other letters are skipped).
        k: word length. Longer k = fewer random hits but more sensitive to errors.
    """

    def __init__(self, reference: str, k: int = 15) -> None:
        self.k = k
        self.reference = reference.upper()
        self.table: dict[str, list[int]] = defaultdict(list)
        for pos, kmer in iter_kmers(self.reference, k):
            self.table[kmer].append(pos)
        self.table = dict(self.table)

    def __len__(self) -> int:
        """Number of distinct k-mers stored."""
        return len(self.table)

    def __contains__(self, kmer: str) -> bool:
        return kmer in self.table

    def lookup(self, kmer: str) -> list[int]:
        """All reference positions where ``kmer`` starts (empty list if none)."""
        return self.table.get(kmer.upper(), [])

    def seeds(self, read: str, step: int = 1, max_hits: int = 50) -> list[tuple[int, int]]:
        """Exact k-mer matches between ``read`` and the reference.

        Returns ``(read_pos, ref_pos)`` pairs. k-mers that occur more than
        ``max_hits`` times (repeats) are ignored because they are uninformative.
        """
        hits = []
        for read_pos, kmer in iter_kmers(read.upper(), self.k):
            if read_pos % step:
                continue
            positions = self.table.get(kmer)
            if positions and len(positions) <= max_hits:
                hits.extend((read_pos, ref_pos) for ref_pos in positions)
        return hits
