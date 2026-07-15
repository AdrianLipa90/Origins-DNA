import pytest

from origins_dna import EncodedDNA, decode_2bit, encode_2bit, normalize_dna, sequence_sha256


@pytest.mark.parametrize(
    "sequence",
    ["", "A", "AC", "ACG", "ACGT", "ACGTA", "ACGTACGTAC"],
)
def test_2bit_round_trip(sequence: str) -> None:
    encoded = encode_2bit(sequence)
    assert decode_2bit(encoded) == sequence


def test_normalization_removes_whitespace_and_uppercases() -> None:
    assert normalize_dna(" acg\nt ") == "ACGT"


def test_invalid_symbol_is_rejected() -> None:
    with pytest.raises(ValueError, match="invalid DNA symbols"):
        normalize_dna("ACNT")


def test_payload_length_is_checked() -> None:
    with pytest.raises(ValueError, match="expected"):
        EncodedDNA(length=5, payload=b"\x00")


def test_sequence_hash_uses_canonical_sequence() -> None:
    assert sequence_sha256(" acgt ") == sequence_sha256("ACGT")
