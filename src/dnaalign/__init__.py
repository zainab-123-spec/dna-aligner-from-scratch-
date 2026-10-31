"""dnaalign: sequence alignment and read mapping written from scratch.

Everything here is plain Python so that each step of the algorithms can be
read, printed and understood. The main entry points are:

* :func:`dnaalign.align.global_align` - Needleman-Wunsch (whole sequences)
* :func:`dnaalign.align.local_align` - Smith-Waterman (best matching piece)
* :class:`dnaalign.kmer_index.KmerIndex` - hash map from k-mer to positions
* :class:`dnaalign.mapper.ReadMapper` - seed-and-extend short-read mapper
"""

from dnaalign.align import Alignment, global_align, local_align, score_alignment
from dnaalign.kmer_index import KmerIndex
from dnaalign.mapper import MappingResult, ReadMapper
from dnaalign.scoring import Scoring, blosum62, dna_scoring

__version__ = "0.1.0"

__all__ = [
    "Alignment",
    "KmerIndex",
    "MappingResult",
    "ReadMapper",
    "Scoring",
    "blosum62",
    "dna_scoring",
    "global_align",
    "local_align",
    "score_alignment",
]
