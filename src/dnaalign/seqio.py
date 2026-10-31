"""Tiny FASTA and FASTQ readers/writers.

FASTA::

    >name optional description
    ACGTACGT...
    ACGT...          (sequence may continue on many lines)

FASTQ (always 4 lines per read)::

    @name
    ACGT...
    +
    IIII...          (one quality letter per base)
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import IO


@dataclass
class Record:
    """One sequence. ``quality`` is only set for FASTQ records."""

    name: str
    sequence: str
    description: str = ""
    quality: str | None = None


def _open(path_or_handle: str | Path | IO[str]) -> IO[str]:
    if isinstance(path_or_handle, (str, Path)):
        return open(path_or_handle, encoding="ascii")
    return path_or_handle


def read_fasta(source: str | Path | IO[str]) -> Iterator[Record]:
    """Yield every record of a FASTA file (path or open text handle)."""
    handle = _open(source)
    try:
        name, desc, chunks = None, "", []
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if name is not None:
                    yield Record(name, "".join(chunks).upper(), desc)
                header = line[1:].split(maxsplit=1)
                name = header[0] if header else ""
                desc = header[1] if len(header) > 1 else ""
                chunks = []
            else:
                if name is None:
                    raise ValueError("FASTA file must start with a '>' header line")
                chunks.append(line)
        if name is not None:
            yield Record(name, "".join(chunks).upper(), desc)
    finally:
        if handle is not source:
            handle.close()


def read_fastq(source: str | Path | IO[str]) -> Iterator[Record]:
    """Yield every record of a FASTQ file (4 lines per record)."""
    handle = _open(source)
    try:
        lines = (line.rstrip("\n\r") for line in handle)
        for header in lines:
            if not header.strip():
                continue
            if not header.startswith("@"):
                raise ValueError(f"expected '@' header line, got {header[:30]!r}")
            seq = next(lines, None)
            plus = next(lines, None)
            qual = next(lines, None)
            if seq is None or plus is None or qual is None:
                raise ValueError(f"truncated FASTQ record {header!r}")
            if not plus.startswith("+"):
                raise ValueError(f"expected '+' line in record {header!r}")
            if len(qual) != len(seq):
                raise ValueError(f"sequence and quality lengths differ in {header!r}")
            parts = header[1:].split(maxsplit=1)
            yield Record(parts[0], seq.upper(), parts[1] if len(parts) > 1 else "", qual)
    finally:
        if handle is not source:
            handle.close()


def read_sequences(path: str | Path) -> Iterator[Record]:
    """Read FASTA or FASTQ, guessing the format from the first character."""
    with open(path, encoding="ascii") as handle:
        first = handle.read(1)
    if first == "@":
        return read_fastq(path)
    if first == ">":
        return read_fasta(path)
    raise ValueError(f"{path}: not a FASTA ('>') or FASTQ ('@') file")


def write_fasta(records: Iterable[Record], handle: IO[str], width: int = 70) -> None:
    """Write records as FASTA, wrapping sequences at ``width`` letters."""
    for rec in records:
        header = f">{rec.name} {rec.description}".rstrip()
        handle.write(header + "\n")
        for i in range(0, len(rec.sequence), width):
            handle.write(rec.sequence[i : i + width] + "\n")


def write_fastq(records: Iterable[Record], handle: IO[str]) -> None:
    """Write records as FASTQ (missing qualities become 'I', i.e. Phred 40)."""
    for rec in records:
        qual = rec.quality if rec.quality is not None else "I" * len(rec.sequence)
        header = f"@{rec.name} {rec.description}".rstrip()
        handle.write(f"{header}\n{rec.sequence}\n+\n{qual}\n")
