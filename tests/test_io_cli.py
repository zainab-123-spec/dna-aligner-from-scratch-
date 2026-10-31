import io

import pytest

from dnaalign.cli import main
from dnaalign.seqio import Record, read_fasta, read_fastq, read_sequences, write_fasta, write_fastq


def test_fasta_multiline():
    text = ">seq1 first one\nACGT\nacgt\n\n>seq2\nTTTT\n"
    recs = list(read_fasta(io.StringIO(text)))
    assert [r.name for r in recs] == ["seq1", "seq2"]
    assert recs[0].sequence == "ACGTACGT"
    assert recs[0].description == "first one"


def test_fasta_without_header_fails():
    with pytest.raises(ValueError):
        list(read_fasta(io.StringIO("ACGT\n")))


def test_fastq_roundtrip():
    buf = io.StringIO()
    write_fastq([Record("r1", "ACGT", quality="IIII"), Record("r2", "GG")], buf)
    recs = list(read_fastq(io.StringIO(buf.getvalue())))
    assert [(r.name, r.sequence, r.quality) for r in recs] == [("r1", "ACGT", "IIII"), ("r2", "GG", "II")]


def test_fastq_quality_length_checked():
    with pytest.raises(ValueError):
        list(read_fastq(io.StringIO("@r\nACGT\n+\nII\n")))


def test_fastq_truncated():
    with pytest.raises(ValueError):
        list(read_fastq(io.StringIO("@r\nACGT\n")))


def test_fasta_roundtrip_and_autodetect(tmp_path):
    path = tmp_path / "x.fa"
    with open(path, "w") as fh:
        write_fasta([Record("chr", "ACGT" * 50, "test")], fh, width=60)
    rec = next(read_sequences(path))
    assert rec.sequence == "ACGT" * 50


def test_cli_align(capsys):
    assert main(["align", "GATTACA", "GCATGCU", "--gap-open", "-1", "--gap-extend", "-1"]) == 0
    assert "score = 0" in capsys.readouterr().out


def test_cli_align_protein(capsys):
    main(["align", "HEAGAWGHEE", "PAWHEAE", "--mode", "local", "--matrix", "blosum62"])
    assert "score" in capsys.readouterr().out


def test_cli_simulate_and_map(tmp_path, reference):
    ref = tmp_path / "ref.fa"
    with open(ref, "w") as fh:
        write_fasta([Record("ref", reference)], fh)
    reads = tmp_path / "reads.fq"
    out = tmp_path / "hits.tsv"
    main(
        [
            "simulate",
            "--ref",
            str(ref),
            "-n",
            "25",
            "--length",
            "100",
            "--error-rate",
            "0.01",
            "--out",
            str(reads),
        ]
    )
    assert len(list(read_fastq(reads))) == 25
    main(["map", "--ref", str(ref), "--reads", str(reads), "-k", "13", "--out", str(out)])
    lines = out.read_text().splitlines()
    assert lines[0].startswith("read\tmapped")
    assert len(lines) == 26
    mapped = [line for line in lines[1:] if line.split("\t")[1] == "1"]
    for line in mapped:  # names look like read3_pos1234_strand+
        name, _, strand, pos = line.split("\t")[:4]
        true_pos = int(name.split("_pos")[1].split("_")[0])
        assert abs(int(pos) - 1 - true_pos) <= 10
        assert strand == name[-1]
    assert len(mapped) >= 24
