from dnaalign.kmer_index import reverse_complement
from dnaalign.mapper import ReadMapper
from dnaalign.simulate import mutate, simulate_reads


def test_exact_read_forward(reference):
    m = ReadMapper(reference, k=13)
    r = m.map_read("r", reference[2000:2150])
    assert r.mapped and r.strand == "+" and r.position == 2000
    assert r.cigar == "150M"
    assert r.score == 300
    assert r.mapq == 60


def test_exact_read_reverse(reference):
    m = ReadMapper(reference, k=13)
    r = m.map_read("r", reverse_complement(reference[3100:3250]))
    assert r.mapped and r.strand == "-" and r.position == 3100


def test_read_with_insertion_and_deletion(reference):
    frag = reference[500:650]
    read = frag[:40] + "GG" + frag[40:100] + frag[103:]
    r = ReadMapper(reference, k=13).map_read("r", read)
    assert r.mapped and abs(r.position - 500) <= 2
    assert "I" in r.cigar and "D" in r.cigar


def test_random_read_is_unmapped(reference, rng):
    junk = "".join(rng.choice("ACGT") for _ in range(150))
    assert not ReadMapper(reference, k=15).map_read("junk", junk).mapped


def test_simulated_reads_map_correctly(reference):
    m = ReadMapper(reference, k=13)
    reads = simulate_reads(reference, 60, 120, error_rate=0.03, seed=3)
    correct = 0
    for read in reads:
        res = m.map_read(read.name, read.sequence)
        correct += res.mapped and res.strand == read.strand and abs(res.position - read.true_position) <= 10
    assert correct >= 57


def test_simulator_truth_without_errors(reference):
    for read in simulate_reads(reference, 20, 100, error_rate=0.0, seed=1):
        frag = reference[read.true_position : read.true_position + 100]
        expected = frag if read.strand == "+" else reverse_complement(frag)
        assert read.sequence == expected
        assert read.n_errors == 0


def test_mutate_counts_errors(rng):
    seq = "ACGT" * 250
    out, n = mutate(seq, 0.05, 0.0, rng)
    assert len(out) == len(seq)
    assert n == sum(a != b for a, b in zip(seq, out))
    assert 20 <= n <= 90
