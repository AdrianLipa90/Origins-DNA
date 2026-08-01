"""Three-quarter-scaled causal readouts for Origins-DNA v0.9.1."""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import math
from typing import Callable, Iterable, Sequence

from .causal_dhse import (
    WINDOWS,
    causal_trajectory,
    prefix_products,
    projective_image,
    terminal_product,
)
from .dhse_exact import LENGTHS, exact_forcing_response, exact_profile
from .encoding import normalize_dna

THREE_QUARTER = 3.0 / 4.0
PRIMARY_WINDOW = 4
POOL_ORDER = ("mean", "minimum", "maximum", "signed_absmax")


def three_quarter_aggregate(
    values: Iterable[float], *, exponent: float = THREE_QUARTER
) -> float:
    """Return len(values)^(-exponent) times the signed sum."""
    data = tuple(float(value) for value in values)
    if not data:
        return 0.0
    if not math.isfinite(exponent) or exponent <= 0.0:
        raise ValueError("exponent must be finite and positive")
    if any(not math.isfinite(value) for value in data):
        raise ValueError("values must be finite")
    return math.fsum(data) / (len(data) ** exponent)


def signed_absmax(values: Iterable[float]) -> float:
    """Return the signed element with greatest absolute magnitude."""
    data = tuple(float(value) for value in values)
    if not data:
        return 0.0
    if any(not math.isfinite(value) for value in data):
        raise ValueError("values must be finite")
    return max(data, key=lambda value: (abs(value), value))


def signed_three_quarter_amplitude(value: float) -> float:
    """Dimensionless signed |x|^(3/4) transform."""
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ValueError("value must be finite")
    if numeric == 0.0:
        return 0.0
    return math.copysign(abs(numeric) ** THREE_QUARTER, numeric)


@lru_cache(maxsize=1)
def _self_dual_references() -> tuple[int, ...]:
    return tuple(exact_profile(length).evaluate(Fraction(1)) for length in LENGTHS)


def self_dual_response(
    response: Sequence[int], *, references: Sequence[int] | None = None
) -> tuple[float, ...]:
    """Center exact Stage-M responses at q=1 and apply 3/4 amplitude."""
    if len(response) != len(LENGTHS):
        raise ValueError("response must contain one value per Stage-M length")
    anchors = (
        _self_dual_references()
        if references is None
        else tuple(int(value) for value in references)
    )
    if len(anchors) != len(LENGTHS) or any(value <= 0 for value in anchors):
        raise ValueError("references must be positive and match Stage-M lengths")
    return tuple(
        signed_three_quarter_amplitude((float(value) - anchor) / anchor)
        for value, anchor in zip(response, anchors)
    )


def phase_window_features(
    sequence: str, *, exponent: float = THREE_QUARTER
) -> tuple[float, ...]:
    """Return the frozen six-coordinate reciprocal-odd window readout."""
    products = prefix_products(sequence)
    if not products:
        return (0.0,) * 6
    logs = tuple(
        math.log(q.numerator) - math.log(q.denominator)
        for q in (projective_image(product) for product in products)
    )
    deltas = tuple(right - left for left, right in zip(logs, logs[1:]))
    return (
        three_quarter_aggregate(logs, exponent=exponent),
        three_quarter_aggregate(deltas, exponent=exponent),
        logs[-1],
        signed_absmax(logs),
        min(logs),
        max(logs),
    )


def even_window_features(sequence: str) -> tuple[float, ...]:
    """Return the frozen reciprocal-even DHSE geometry window readout."""
    states = causal_trajectory(sequence)
    if not states:
        return (0.0,) * 17
    transformed = tuple(self_dual_response(state.response) for state in states)
    output: list[float] = []
    for index in range(len(LENGTHS)):
        channel = tuple(row[index] for row in transformed)
        output.extend(
            (
                three_quarter_aggregate(channel),
                channel[-1],
                min(channel),
                max(channel),
            )
        )
    residuals = tuple(float(state.residual) for state in states)
    output.append(three_quarter_aggregate(residuals))
    return tuple(output)


def _pool_feature_rows(rows: Sequence[Sequence[float]], width: int) -> tuple[float, ...]:
    if not rows:
        return (0.0,) * (width * len(POOL_ORDER))
    if any(len(row) != width for row in rows):
        raise ValueError("feature rows have inconsistent width")
    output: list[float] = []
    for index in range(width):
        column = tuple(float(row[index]) for row in rows)
        output.extend(
            (
                math.fsum(column) / len(column),
                min(column),
                max(column),
                signed_absmax(column),
            )
        )
    return tuple(output)


def _pooled_length_block(
    sequence: str,
    length: int,
    builder: Callable[[str], tuple[float, ...]],
    width: int,
) -> tuple[float, ...]:
    normalized = normalize_dna(sequence)
    if length <= 0:
        raise ValueError("window length must be positive")
    rows = tuple(
        builder(normalized[start : start + length])
        for start in range(0, max(0, len(normalized) - length + 1))
    )
    return _pool_feature_rows(rows, width)


def _multiscale_vector(
    sequence: str,
    builder: Callable[[str], tuple[float, ...]],
    width: int,
) -> tuple[float, ...]:
    normalized = normalize_dna(sequence)
    output = list(_pooled_length_block(normalized, PRIMARY_WINDOW, builder, width))
    for length in WINDOWS:
        output.extend(_pooled_length_block(normalized, length, builder, width))
    return tuple(output)


def phase_three_quarter_vector(sequence: str) -> tuple[float, ...]:
    return _multiscale_vector(sequence, phase_window_features, 6)


def phase_ordinary_mean_vector(sequence: str) -> tuple[float, ...]:
    return _multiscale_vector(
        sequence,
        lambda window: phase_window_features(window, exponent=1.0),
        6,
    )


def even_three_quarter_vector(sequence: str) -> tuple[float, ...]:
    return _multiscale_vector(sequence, even_window_features, 17)


def exact_three_quarter_gate_receipt(
    motif: str = "CGCG", decoy: str = "GCGC"
) -> dict[str, object]:
    motif_normalized = normalize_dna(motif)
    decoy_normalized = normalize_dna(decoy)
    if sorted(motif_normalized) != sorted(decoy_normalized):
        raise ValueError("motif and decoy must be composition matched")
    q_motif = projective_image(terminal_product(motif_normalized))
    q_decoy = projective_image(terminal_product(decoy_normalized))
    response_motif = exact_forcing_response(q_motif)
    response_decoy = exact_forcing_response(q_decoy)
    reciprocal = q_motif * q_decoy == 1
    opposite_phase = reciprocal and q_motif != 1 and q_decoy != 1
    even_equal = response_motif == response_decoy
    passed = reciprocal and opposite_phase and even_equal
    return {
        "schema": "origins-dna.causal-three-quarter-gates.v1",
        "technical_status": "PASS" if passed else "FAIL",
        "motif": motif_normalized,
        "decoy": decoy_normalized,
        "q_motif": [q_motif.numerator, q_motif.denominator],
        "q_decoy": [q_decoy.numerator, q_decoy.denominator],
        "reciprocal_terminal_coordinates": reciprocal,
        "opposite_nonzero_signed_phase": opposite_phase,
        "equal_stage_m_response": even_equal,
        "stage_m_response": list(response_motif),
        "three_quarter_exponent": THREE_QUARTER,
    }
