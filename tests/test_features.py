import math

import pytest

from origins_dna import base_frequencies, kmer_counts, shannon_entropy_bits


def test_base_frequencies() -> None:
    assert base_frequencies("AACG") == {
        "A": 0.5,
        "C": 0.25,
        "G": 0.25,
        "T": 0.0,
    }


def test_empty_base_frequencies() -> None:
    assert base_frequencies("") == {base: 0.0 for base in "ACGT"}


def test_overlapping_kmer_counts() -> None:
    assert kmer_counts("AAAA", 2) == {"AA": 3}


def test_invalid_k_is_rejected() -> None:
    with pytest.raises(ValueError):
        kmer_counts("ACGT", 0)


def test_entropy_uniform_distribution() -> None:
    assert math.isclose(shannon_entropy_bits("ACGT"), 2.0)
