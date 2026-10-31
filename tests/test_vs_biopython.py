"""Our scores must equal Biopython's PairwiseAligner (an independent implementation)."""

import random

import pytest

from dnaalign import Scoring, blosum62, global_align, local_align

Align = pytest.importorskip("Bio.Align")


def _bio(mode, sc, protein=False):
    al = Align.PairwiseAligner()
    al.mode = mode
    if protein:
        al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    else:
        al.match_score, al.mismatch_score = sc.match, sc.mismatch
    al.open_gap_score, al.extend_gap_score = sc.gap_open, sc.gap_extend
    return al


SCHEMES = [Scoring(1, -1, -2, -2), Scoring(2, -3, -5, -2), Scoring(1, -2, -3, -1), Scoring(5, -4, -10, -0.5)]


@pytest.mark.parametrize("mode", ["global", "local"])
@pytest.mark.parametrize("sc", SCHEMES, ids=["linear", "affine1", "affine2", "affine_frac"])
def test_dna_random_cases(mode, sc):
    rng = random.Random(hash((mode, sc.gap_open)) % 1000)
    bio = _bio(mode, sc)
    fn = global_align if mode == "global" else local_align
    for _ in range(60):
        a = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 60)))
        b = "".join(c if rng.random() > 0.25 else rng.choice("ACGT") for c in a)
        b = b[rng.randint(0, 10) :] or "A"
        assert fn(a, b, sc).score == pytest.approx(bio.score(a, b))


@pytest.mark.parametrize("mode", ["global", "local"])
def test_protein_blosum62_random_cases(mode):
    rng = random.Random(7)
    sc = blosum62()
    bio = _bio(mode, sc, protein=True)
    fn = global_align if mode == "global" else local_align
    for _ in range(40):
        a = "".join(rng.choice("ARNDCQEGHILKMFPSTWYV") for _ in range(rng.randint(5, 60)))
        b = "".join(rng.choice("ARNDCQEGHILKMFPSTWYV") for _ in range(rng.randint(5, 60)))
        assert fn(a, b, sc).score == pytest.approx(bio.score(a, b))
