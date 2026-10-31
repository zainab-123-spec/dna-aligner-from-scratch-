"""Pairwise alignment with dynamic programming (DP).

Two classic algorithms:

* **Needleman-Wunsch** (:func:`global_align`): align the *whole* of ``a``
  against the *whole* of ``b``.
* **Smith-Waterman** (:func:`local_align`): find the best-scoring *piece* of
  ``a`` that matches a piece of ``b``. Scores are never allowed to drop below
  zero, so a bad region simply "restarts" the alignment.

Both come in two flavours:

* **linear gaps** - one table ``S``; every gap letter costs the same.
* **affine gaps (Gotoh, 1982)** - three tables, because opening a gap costs
  more than making it longer:

  - ``M[i][j]``: best score where ``a[i-1]`` is paired with ``b[j-1]``
  - ``X[i][j]``: best score where ``a[i-1]`` is paired with a gap
  - ``Y[i][j]``: best score where ``b[j-1]`` is paired with a gap

Filling an ``(n+1) x (m+1)`` table takes **O(n*m) time and memory**.
"""

from __future__ import annotations

from dataclasses import dataclass

from dnaalign.scoring import Scoring

NEG_INF = float("-inf")

# Traceback arrows for the one-table (linear gap) version.
STOP, DIAG, UP, LEFT = 0, 1, 2, 3
# State names for the three-table (affine gap) version.
M_STATE, X_STATE, Y_STATE, START = 0, 1, 2, 3


@dataclass
class Alignment:
    """Result of an alignment.

    ``aligned_a`` and ``aligned_b`` have the same length; ``-`` marks a gap.
    ``start_*``/``end_*`` are 0-based, half-open coordinates of the aligned
    part of each input (for a global alignment they cover everything).
    """

    score: float
    aligned_a: str
    aligned_b: str
    start_a: int
    end_a: int
    start_b: int
    end_b: int

    def cigar(self) -> str:
        """CIGAR string with ``a`` as the query and ``b`` as the reference.

        ``M`` = letter vs letter, ``I`` = letter in ``a`` only (insertion),
        ``D`` = letter in ``b`` only (deletion). Example: ``"10M2I5M"``.
        """
        ops = []
        for x, y in zip(self.aligned_a, self.aligned_b):
            ops.append("I" if y == "-" else "D" if x == "-" else "M")
        out, i = [], 0
        while i < len(ops):
            j = i
            while j < len(ops) and ops[j] == ops[i]:
                j += 1
            out.append(f"{j - i}{ops[i]}")
            i = j
        return "".join(out)

    @property
    def identity(self) -> float:
        """Fraction of alignment columns where the two letters are identical."""
        if not self.aligned_a:
            return 0.0
        same = sum(1 for x, y in zip(self.aligned_a, self.aligned_b) if x == y and x != "-")
        return same / len(self.aligned_a)

    def __str__(self) -> str:
        middle = "".join(
            "|" if x == y and x != "-" else " " if "-" in (x, y) else "."
            for x, y in zip(self.aligned_a, self.aligned_b)
        )
        return f"{self.aligned_a}\n{middle}\n{self.aligned_b}\nscore = {self.score:g}"


def score_alignment(aligned_a: str, aligned_b: str, scoring: Scoring) -> float:
    """Re-score an existing alignment column by column (useful for checking)."""
    if len(aligned_a) != len(aligned_b):
        raise ValueError("aligned strings must have equal length")
    total = 0.0
    prev = None  # which side had a gap in the previous column
    for x, y in zip(aligned_a, aligned_b):
        if x == "-" and y == "-":
            raise ValueError("a column cannot be gap/gap")
        if x == "-" or y == "-":
            side = "a" if x == "-" else "b"
            total += scoring.gap_extend if prev == side else scoring.gap_open
            prev = side
        else:
            total += scoring.pair(x, y)
            prev = None
    return total


def global_align(a: str, b: str, scoring: Scoring | None = None) -> Alignment:
    """Needleman-Wunsch global alignment of ``a`` and ``b``.

    Uses the simple one-table algorithm for linear gaps and Gotoh's
    three-table algorithm for affine gaps.

    >>> global_align("GATTACA", "GCATGCU", Scoring(1, -1, -1, -1)).score
    0.0
    """
    scoring = scoring or Scoring()
    if scoring.is_linear:
        return _linear(a, b, scoring, local=False)
    return _affine(a, b, scoring, local=False)


def local_align(
    a: str,
    b: str,
    scoring: Scoring | None = None,
    band: int | None = None,
    band_offset: int = 0,
) -> Alignment:
    """Smith-Waterman local alignment of ``a`` and ``b``.

    Args:
        a, b: the two sequences.
        scoring: scoring scheme (default: match 1, mismatch -1, gap -2).
        band: if given, only fill cells with ``|j - i - band_offset| <= band``.
            This is a *banded* alignment: it costs O(n * band) instead of
            O(n * m) and is what the read mapper uses after a seed hit tells it
            roughly which diagonal the read lies on.
        band_offset: the expected diagonal ``j - i`` (only used with ``band``).
    """
    scoring = scoring or Scoring()
    return (
        _affine(a, b, scoring, local=True, band=band, band_offset=band_offset)
        if (band is not None or not scoring.is_linear)
        else _linear(a, b, scoring, local=True)
    )


# ---------------------------------------------------------------------------
# Linear gaps: one score table S and one arrow table T.
# ---------------------------------------------------------------------------
def _linear(a: str, b: str, sc: Scoring, local: bool) -> Alignment:
    n, m, gap = len(a), len(b), sc.gap_open
    S = [[0.0] * (m + 1) for _ in range(n + 1)]
    T = [[STOP] * (m + 1) for _ in range(n + 1)]
    if not local:  # first row/column = "everything so far aligned to gaps"
        for i in range(1, n + 1):
            S[i][0], T[i][0] = i * gap, UP
        for j in range(1, m + 1):
            S[0][j], T[0][j] = j * gap, LEFT

    best, best_i, best_j = 0.0, 0, 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = S[i - 1][j - 1] + sc.pair(a[i - 1], b[j - 1])
            up = S[i - 1][j] + gap  # a[i-1] against a gap
            left = S[i][j - 1] + gap  # b[j-1] against a gap
            score, arrow = diag, DIAG
            if up > score:
                score, arrow = up, UP
            if left > score:
                score, arrow = left, LEFT
            if local and score <= 0:  # Smith-Waterman: never go below zero
                score, arrow = 0.0, STOP
            S[i][j], T[i][j] = score, arrow
            if local and score > best:
                best, best_i, best_j = score, i, j

    i, j = (best_i, best_j) if local else (n, m)
    end_a, end_b = i, j
    top, bottom = [], []
    while (i > 0 or j > 0) and T[i][j] != STOP:
        arrow = T[i][j]
        if arrow == DIAG:
            top.append(a[i - 1])
            bottom.append(b[j - 1])
            i, j = i - 1, j - 1
        elif arrow == UP:
            top.append(a[i - 1])
            bottom.append("-")
            i -= 1
        else:
            top.append("-")
            bottom.append(b[j - 1])
            j -= 1
    final = best if local else S[n][m]
    return Alignment(final, "".join(reversed(top)), "".join(reversed(bottom)), i, end_a, j, end_b)


# ---------------------------------------------------------------------------
# Affine gaps (Gotoh): three score tables M, X, Y plus one arrow table each.
# An arrow stores which table we came FROM in the previous cell.
# ---------------------------------------------------------------------------
def _best_of(*pairs: tuple[float, int]) -> tuple[float, int]:
    """Return the (score, state) pair with the highest score (first wins ties)."""
    best = pairs[0]
    for p in pairs[1:]:
        if p[0] > best[0]:
            best = p
    return best


def _affine(
    a: str,
    b: str,
    sc: Scoring,
    local: bool,
    band: int | None = None,
    band_offset: int = 0,
) -> Alignment:
    n, m = len(a), len(b)
    go, ge = sc.gap_open, sc.gap_extend
    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    X = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    tM = [[START] * (m + 1) for _ in range(n + 1)]
    tX = [[START] * (m + 1) for _ in range(n + 1)]
    tY = [[START] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0.0
    if not local:
        for i in range(1, n + 1):  # a[:i] against one long gap
            X[i][0] = go + (i - 1) * ge
            tX[i][0] = M_STATE if i == 1 else X_STATE
        for j in range(1, m + 1):
            Y[0][j] = go + (j - 1) * ge
            tY[0][j] = M_STATE if j == 1 else Y_STATE

    best, best_i, best_j = 0.0, 0, 0
    for i in range(1, n + 1):
        if band is None:
            j_lo, j_hi = 1, m
        else:
            j_lo = max(1, i + band_offset - band)
            j_hi = min(m, i + band_offset + band)
        Mi, Xi, Yi, Mp, Xp, Yp = M[i], X[i], Y[i], M[i - 1], X[i - 1], Y[i - 1]
        ai = a[i - 1]
        for j in range(j_lo, j_hi + 1):
            # M: a[i-1] paired with b[j-1], coming diagonally from any table.
            prev, state = _best_of((Mp[j - 1], M_STATE), (Xp[j - 1], X_STATE), (Yp[j - 1], Y_STATE))
            if local and prev <= 0:
                prev, state = 0.0, START  # start a fresh local alignment here
            Mi[j] = prev + sc.pair(ai, b[j - 1])
            tM[i][j] = state
            # X: a[i-1] paired with a gap (move down). Open from M/Y, extend X.
            Xi[j], tX[i][j] = _best_of((Mp[j] + go, M_STATE), (Xp[j] + ge, X_STATE), (Yp[j] + go, Y_STATE))
            # Y: b[j-1] paired with a gap (move right). Open from M/X, extend Y.
            Yi[j], tY[i][j] = _best_of(
                (Mi[j - 1] + go, M_STATE), (Yi[j - 1] + ge, Y_STATE), (Xi[j - 1] + go, X_STATE)
            )
            if local and Mi[j] > best:
                best, best_i, best_j = Mi[j], i, j

    if local:
        if best <= 0:
            return Alignment(0.0, "", "", 0, 0, 0, 0)
        i, j, state, final = best_i, best_j, M_STATE, best
    else:
        final, state = _best_of((M[n][m], M_STATE), (X[n][m], X_STATE), (Y[n][m], Y_STATE))
        i, j = n, m
        if n == 0 and m == 0:
            return Alignment(0.0, "", "", 0, 0, 0, 0)

    end_a, end_b = i, j
    top, bottom = [], []
    while state != START and (i > 0 or j > 0):
        if state == M_STATE:
            top.append(a[i - 1])
            bottom.append(b[j - 1])
            state = tM[i][j]
            i, j = i - 1, j - 1
        elif state == X_STATE:
            top.append(a[i - 1])
            bottom.append("-")
            state = tX[i][j]
            i -= 1
        else:
            top.append("-")
            bottom.append(b[j - 1])
            state = tY[i][j]
            j -= 1
    return Alignment(final, "".join(reversed(top)), "".join(reversed(bottom)), i, end_a, j, end_b)
