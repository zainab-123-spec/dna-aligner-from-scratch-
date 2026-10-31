import pytest

from dnaalign.scoring import Scoring, blosum62, dna_scoring, parse_matrix


def test_match_and_mismatch():
    sc = Scoring(match=2, mismatch=-3)
    assert sc.pair("A", "A") == 2
    assert sc.pair("A", "C") == -3


def test_pair_is_case_insensitive():
    assert Scoring().pair("a", "A") == Scoring().match


def test_gap_lengths_affine():
    sc = Scoring(gap_open=-5, gap_extend=-1)
    assert sc.gap(0) == 0
    assert sc.gap(1) == -5
    assert sc.gap(4) == -8
    assert not sc.is_linear


def test_linear_flag():
    assert Scoring(gap_open=-2, gap_extend=-2).is_linear


def test_positive_gap_rejected():
    with pytest.raises(ValueError):
        Scoring(gap_open=1)


def test_blosum62_known_values():
    sc = blosum62()
    assert sc.pair("W", "W") == 11
    assert sc.pair("A", "A") == 4
    assert sc.pair("W", "C") == -2
    assert sc.pair("E", "Q") == sc.pair("Q", "E") == 2


def test_blosum62_matches_biopython():
    sm = pytest.importorskip("Bio.Align.substitution_matrices")
    ref = sm.load("BLOSUM62")
    ours = blosum62()
    for x in ref.alphabet:
        for y in ref.alphabet:
            assert ours.pair(x, y) == ref[x][y]


def test_unknown_letter_in_matrix():
    with pytest.raises(KeyError):
        blosum62().pair("A", "J")


def test_parse_matrix_bad_row():
    with pytest.raises(ValueError):
        parse_matrix("   A C\nA 1 0\nC 0\n")


def test_dna_scoring_defaults():
    sc = dna_scoring()
    assert (sc.match, sc.mismatch, sc.gap_open, sc.gap_extend) == (2, -4, -6, -1)
