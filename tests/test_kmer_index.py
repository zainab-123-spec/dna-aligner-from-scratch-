import pytest

from dnaalign.kmer_index import KmerIndex, iter_kmers, reverse_complement


def test_small_index():
    idx = KmerIndex("ACGTACG", k=3)
    assert idx.lookup("ACG") == [0, 4]
    assert idx.lookup("GTA") == [2]
    assert idx.lookup("TTT") == []
    assert len(idx) == 4
    assert "CGT" in idx


def test_kmers_with_n_are_skipped():
    assert [k for _, k in iter_kmers("ACNGT", 2)] == ["AC", "GT"]


def test_bad_k():
    with pytest.raises(ValueError):
        list(iter_kmers("ACGT", 0))


def test_reverse_complement():
    assert reverse_complement("AACG") == "CGTT"
    assert reverse_complement(reverse_complement("GATTACA")) == "GATTACA"


def test_every_position_is_indexed(reference):
    idx = KmerIndex(reference, k=12)
    for pos in range(0, len(reference) - 12, 97):
        assert pos in idx.lookup(reference[pos : pos + 12])


def test_seeds_report_diagonal(reference):
    idx = KmerIndex(reference, k=11)
    read = reference[1000:1100]
    diagonals = {ref - rp for rp, ref in idx.seeds(read)}
    assert 1000 in diagonals


def test_repeats_are_filtered():
    idx = KmerIndex("A" * 100, k=5)
    assert idx.seeds("AAAAAAAA", max_hits=10) == []
