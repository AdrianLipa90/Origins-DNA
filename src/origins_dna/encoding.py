"""Strict normalization and reversible 2-bit encoding for DNA."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib

_BASE_TO_BITS = {"A": 0b00, "C": 0b01, "G": 0b10, "T": 0b11}
_BITS_TO_BASE = ("A", "C", "G", "T")


@dataclass(frozen=True, slots=True)
class EncodedDNA:
    """Packed DNA sequence with explicit original nucleotide length."""

    length: int
    payload: bytes

    def __post_init__(self) -> None:
        if self.length < 0:
            raise ValueError("length must be non-negative")
        expected_bytes = (self.length + 3) // 4
        if len(self.payload) != expected_bytes:
            raise ValueError(
                f"payload has {len(self.payload)} bytes; expected {expected_bytes} "
                f"for {self.length} nucleotides"
            )


def normalize_dna(sequence: str) -> str:
    """Upper-case a sequence, remove whitespace and reject non-ACGT symbols."""

    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")

    normalized = "".join(sequence.split()).upper()
    invalid = sorted(set(normalized) - set(_BASE_TO_BITS))
    if invalid:
        raise ValueError(f"invalid DNA symbols: {''.join(invalid)}")
    return normalized


def encode_2bit(sequence: str) -> EncodedDNA:
    """Pack four DNA bases into each byte, preserving exact sequence length."""

    normalized = normalize_dna(sequence)
    payload = bytearray((len(normalized) + 3) // 4)

    for index, base in enumerate(normalized):
        byte_index = index // 4
        shift = 6 - 2 * (index % 4)
        payload[byte_index] |= _BASE_TO_BITS[base] << shift

    return EncodedDNA(length=len(normalized), payload=bytes(payload))


def decode_2bit(encoded: EncodedDNA) -> str:
    """Decode an :class:`EncodedDNA` value back to its canonical sequence."""

    if not isinstance(encoded, EncodedDNA):
        raise TypeError("encoded must be an EncodedDNA instance")

    bases: list[str] = []
    for index in range(encoded.length):
        byte_index = index // 4
        shift = 6 - 2 * (index % 4)
        value = (encoded.payload[byte_index] >> shift) & 0b11
        bases.append(_BITS_TO_BASE[value])
    return "".join(bases)


def sequence_sha256(sequence: str) -> str:
    """Return SHA-256 of the canonical normalized DNA sequence."""

    normalized = normalize_dna(sequence)
    return hashlib.sha256(normalized.encode("ascii")).hexdigest()
