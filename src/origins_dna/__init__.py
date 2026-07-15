"""Origins-DNA public API."""

from .encoding import EncodedDNA, decode_2bit, encode_2bit, normalize_dna, sequence_sha256
from .features import base_frequencies, kmer_counts, shannon_entropy_bits
from .simulation import MutationRecord, hamming_distance, mutate_substitutions

__all__ = [
    "EncodedDNA",
    "MutationRecord",
    "base_frequencies",
    "decode_2bit",
    "encode_2bit",
    "hamming_distance",
    "kmer_counts",
    "mutate_substitutions",
    "normalize_dna",
    "sequence_sha256",
]
