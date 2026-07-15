"""Origins-DNA public API."""

from .encoding import EncodedDNA, decode_2bit, encode_2bit, normalize_dna, sequence_sha256
from .features import base_frequencies, kmer_counts, shannon_entropy_bits
from .orbital_atomic import (
    ConsolidatedWindow,
    OrbitalAtomicState,
    complement,
    consolidate_windows,
    orbital_atomic_encode,
    orbital_atomic_summary,
    reverse_complement,
)
from .simulation import MutationRecord, hamming_distance, mutate_substitutions

__all__ = [
    "ConsolidatedWindow",
    "EncodedDNA",
    "MutationRecord",
    "OrbitalAtomicState",
    "base_frequencies",
    "complement",
    "consolidate_windows",
    "decode_2bit",
    "encode_2bit",
    "hamming_distance",
    "kmer_counts",
    "mutate_substitutions",
    "normalize_dna",
    "orbital_atomic_encode",
    "orbital_atomic_summary",
    "reverse_complement",
    "sequence_sha256",
    "shannon_entropy_bits",
]
