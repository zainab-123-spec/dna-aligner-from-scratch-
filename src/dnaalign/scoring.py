"""Scoring schemes: how many points a column of an alignment is worth.

An alignment is a list of columns. Each column is one of:

* two letters (a match or a mismatch), scored by :meth:`Scoring.pair`
* a letter against a gap ``-``, scored with the gap penalties

Gap scores follow the same convention as Biopython's ``PairwiseAligner``:
a gap of length ``L`` costs ``gap_open + (L - 1) * gap_extend``.
When ``gap_open == gap_extend`` every gap letter costs the same, which is
called a *linear* gap penalty. Otherwise it is an *affine* gap penalty.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from importlib import resources


@dataclass(frozen=True)
class Scoring:
    """Points for matches, mismatches and gaps.

    Args:
        match: score when two letters are equal (DNA mode).
        mismatch: score when two letters differ (DNA mode).
        gap_open: score of the first letter of a gap (should be <= 0).
        gap_extend: score of every further letter of the same gap (<= 0).
        matrix: optional substitution matrix ``{(x, y): score}``. When given,
            it is used instead of ``match`` / ``mismatch`` (e.g. BLOSUM62).
    """

    match: float = 1.0
    mismatch: float = -1.0
    gap_open: float = -2.0
    gap_extend: float = -2.0
    matrix: dict[tuple[str, str], float] | None = field(default=None, compare=False)

    def __post_init__(self) -> None:
        if self.gap_open > 0 or self.gap_extend > 0:
            raise ValueError("gap scores must be zero or negative")

    @property
    def is_linear(self) -> bool:
        """True when every gap letter costs the same (no extra opening cost)."""
        return self.gap_open == self.gap_extend

    def pair(self, x: str, y: str) -> float:
        """Score of putting letter ``x`` on top of letter ``y``."""
        if self.matrix is not None:
            try:
                return self.matrix[(x.upper(), y.upper())]
            except KeyError as err:
                raise KeyError(f"letter pair {x!r}/{y!r} not in substitution matrix") from err
        return self.match if x.upper() == y.upper() else self.mismatch

    def gap(self, length: int) -> float:
        """Score of one gap that is ``length`` letters long."""
        if length <= 0:
            return 0.0
        return self.gap_open + (length - 1) * self.gap_extend


def dna_scoring(
    match: float = 2, mismatch: float = -4, gap_open: float = -6, gap_extend: float = -1
) -> Scoring:
    """A sensible DNA scheme for read mapping (similar in spirit to short-read mappers)."""
    return Scoring(match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend)


def parse_matrix(text: str) -> dict[tuple[str, str], float]:
    """Parse a substitution matrix in the classic NCBI text format.

    The first non-comment line lists the column letters; every next line is a
    row letter followed by one number per column.
    """
    rows = [line.split() for line in text.splitlines() if line.strip() and not line.startswith("#")]
    header = rows[0]
    matrix: dict[tuple[str, str], float] = {}
    for row in rows[1:]:
        letter, values = row[0], row[1:]
        if len(values) != len(header):
            raise ValueError(f"row {letter!r} has {len(values)} values, expected {len(header)}")
        for col, value in zip(header, values):
            matrix[(letter, col)] = float(value)
    return matrix


@lru_cache(maxsize=1)
def _blosum62_table() -> dict[tuple[str, str], float]:
    text = resources.files("dnaalign").joinpath("data/BLOSUM62.txt").read_text()
    return parse_matrix(text)


def blosum62(gap_open: float = -11, gap_extend: float = -1) -> Scoring:
    """Protein scoring with the BLOSUM62 matrix (BLASTP default gaps: 11 / 1)."""
    return Scoring(gap_open=gap_open, gap_extend=gap_extend, matrix=_blosum62_table())
