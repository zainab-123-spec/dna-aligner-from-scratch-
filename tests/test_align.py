import pytest

from dnaalign.align import global_align, local_align, score_alignment
from dnaalign.scoring import Scoring, blosum62

LINEAR = Scoring(1, -1, -2, -2)
AFFINE = Scoring(2, -3, -5, -1)


def test_identical_global():
    aln = global_align("ACGTACGT", "ACGTACGT", LINEAR)
    assert aln.score == 8
    assert aln.aligned_a == aln.aligned_b == "ACGTACGT"
    assert aln.cigar() == "8M"
    assert aln.identity == 1.0


def test_textbook_example():
    # Classic Wikipedia Needleman-Wunsch example: match 1, mismatch -1, gap -1.
    assert global_align("GATTACA", "GCATGCU", Scoring(1, -1, -1, -1)).score == 0


def test_global_single_gap():
    aln = global_align("ACGTTACG", "ACGTACG", LINEAR)
    assert aln.score == 7 - 2
    assert aln.aligned_b.count("-") == 1


def test_global_empty_inputs():
    assert global_align("", "", LINEAR).score == 0
    assert global_align("ACG", "", LINEAR).score == -6
    assert global_align("", "ACGT", AFFINE).score == AFFINE.gap(4)


def test_global_covers_everything():
    aln = global_align("AAACCC", "AAAGGGCCC", AFFINE)
    assert aln.aligned_a.replace("-", "") == "AAACCC"
    assert aln.aligned_b.replace("-", "") == "AAAGGGCCC"
    assert (aln.start_a, aln.end_a, aln.start_b, aln.end_b) == (0, 6, 0, 9)


def test_affine_prefers_one_long_gap():
    # One gap of 3 (-5 -1 -1 = -7) beats three separate gaps (3 * -5 = -15).
    aln = global_align("AAAACCCCGGGG", "AAAAGGGG", AFFINE)
    assert "---" in aln.aligned_b
    assert aln.cigar() == "4M4I4M"


def test_local_finds_common_core():
    aln = local_align("TTTTGATTACATTTT", "CCCGATTACACCC", LINEAR)
    assert aln.score == 7
    assert aln.aligned_a == "GATTACA"
    assert (aln.start_a, aln.end_a) == (4, 11)
    assert (aln.start_b, aln.end_b) == (3, 10)


def test_local_no_similarity_gives_zero():
    aln = local_align("AAAA", "CCCC", LINEAR)
    assert aln.score == 0
    assert aln.aligned_a == ""


def test_local_affine_no_similarity():
    assert local_align("AAAA", "TTTT", AFFINE).score == 0


def test_local_score_never_negative(rng):
    for _ in range(20):
        a = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 30)))
        b = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 30)))
        assert local_align(a, b, AFFINE).score >= 0


def test_linear_and_gotoh_agree_when_open_equals_extend(rng):
    from dnaalign.align import _affine, _linear

    sc = Scoring(1, -1, -2, -2)
    for _ in range(30):
        a = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 25)))
        b = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 25)))
        assert _linear(a, b, sc, local=False).score == _affine(a, b, sc, local=False).score
        assert _linear(a, b, sc, local=True).score == _affine(a, b, sc, local=True).score


@pytest.mark.parametrize("sc", [LINEAR, AFFINE, Scoring(1, -2, -4, -1)])
@pytest.mark.parametrize("mode", ["global", "local"])
def test_traceback_rescores_to_reported_score(rng, sc, mode):
    fn = global_align if mode == "global" else local_align
    for _ in range(25):
        a = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 40)))
        b = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 40)))
        aln = fn(a, b, sc)
        if aln.aligned_a:
            assert score_alignment(aln.aligned_a, aln.aligned_b, sc) == aln.score
            assert aln.aligned_a.replace("-", "") == a[aln.start_a : aln.end_a]
            assert aln.aligned_b.replace("-", "") == b[aln.start_b : aln.end_b]


def test_protein_blosum62_local():
    aln = local_align("HEAGAWGHEE", "PAWHEAE", blosum62(-8, -8))
    assert aln.score > 0
    assert "AW" in aln.aligned_a


def test_banded_equals_full_when_band_is_wide(rng):
    for _ in range(20):
        a = "".join(rng.choice("ACGT") for _ in range(30))
        b = "".join(rng.choice("ACGT") for _ in range(35))
        full = local_align(a, b, AFFINE)
        banded = local_align(a, b, AFFINE, band=100)
        assert full.score == banded.score


def test_banded_only_uses_band():
    # The match lies far off the main diagonal, so a narrow band cannot see it.
    a = "GATTACAGATTACA"
    b = "C" * 40 + a
    assert local_align(a, b, AFFINE, band=3, band_offset=0).score < local_align(a, b, AFFINE).score
    assert local_align(a, b, AFFINE, band=3, band_offset=40).score == 2 * len(a)


def test_score_alignment_rejects_bad_input():
    with pytest.raises(ValueError):
        score_alignment("AC", "A", LINEAR)
    with pytest.raises(ValueError):
        score_alignment("A-", "A-", LINEAR)


def test_str_shows_bars():
    text = str(global_align("ACGT", "ACCT", LINEAR))
    assert "||.|" in text


def test_cigar_with_deletion():
    aln = global_align("AAAAGGGG", "AAAACCCCGGGG", AFFINE)
    assert aln.cigar() == "4M4D4M"
