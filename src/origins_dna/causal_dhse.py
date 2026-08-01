"""Noncommutative 5'->3' causal DHSE trajectories for Origins-DNA v0.9."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import math
from typing import Iterable, Sequence

from .dhse_exact import exact_forcing_response, primitive_maps
from .encoding import normalize_dna
from .orbital_atomic import complement, reverse_complement

Matrix = tuple[int, int, int, int]
BASES = "ACGT"
WINDOWS = (3, 4, 6, 12, 24)
IDENTITY: Matrix = (1, 0, 0, 1)


def reciprocal_conjugate(matrix: Matrix) -> Matrix:
    a, b, c, d = matrix
    return (d, c, b, a)


def multiply(left: Matrix, right: Matrix) -> Matrix:
    """Return left @ right with exact integer coefficients."""
    a, b, c, d = left
    A, B, C, D = right
    return (
        a * A + b * C,
        a * B + b * D,
        c * A + d * C,
        c * B + d * D,
    )


def determinant(matrix: Matrix) -> int:
    a, b, c, d = matrix
    return a * d - b * c


def trace(matrix: Matrix) -> int:
    a, _, _, d = matrix
    return a + d


def projective_image(matrix: Matrix, q: Fraction | int = Fraction(1)) -> Fraction:
    value = Fraction(q)
    a, b, c, d = matrix
    denominator = c * value + d
    if denominator == 0:
        raise ZeroDivisionError("projective image denominator is zero")
    return (a * value + b) / denominator


@lru_cache(maxsize=1)
def base_operator_assignment() -> dict[str, Matrix]:
    """Apply the preregistered lexicographic reciprocal-pair assignment."""
    matrices = sorted(tuple(int(value) for value in row) for row in primitive_maps().tolist())
    matrix_set = set(matrices)
    pairs: list[tuple[Matrix, Matrix]] = []
    seen: set[Matrix] = set()
    for matrix in matrices:
        if matrix in seen:
            continue
        partner = reciprocal_conjugate(matrix)
        if partner == matrix or partner not in matrix_set:
            seen.add(matrix)
            continue
        pair = tuple(sorted((matrix, partner)))
        pairs.append(pair)  # type: ignore[arg-type]
        seen.update(pair)
    pairs.sort(key=lambda pair: pair[0])
    first, second = pairs[:2]
    assignment = {"A": first[0], "T": first[1], "C": second[0], "G": second[1]}
    if len(set(assignment.values())) != 4:
        raise RuntimeError("causal operator assignment is not four-way distinct")
    return assignment


def operator_word(sequence: str) -> tuple[Matrix, ...]:
    normalized = normalize_dna(sequence)
    assignment = base_operator_assignment()
    return tuple(assignment[base] for base in normalized)


def prefix_products(sequence: str) -> tuple[Matrix, ...]:
    """Return M_xk ... M_x1 for all prefixes read 5'->3'."""
    current = IDENTITY
    output: list[Matrix] = []
    for operator in operator_word(sequence):
        current = multiply(operator, current)
        output.append(current)
    return tuple(output)


@dataclass(frozen=True, slots=True)
class CausalState:
    index: int
    base: str
    product: Matrix
    q: Fraction
    signed_log_q: float
    residual: Fraction
    response: tuple[int, int, int, int]


def _residual(q: Fraction) -> Fraction:
    return abs(q - 1) / (q + 1)


def causal_trajectory(sequence: str) -> tuple[CausalState, ...]:
    normalized = normalize_dna(sequence)
    products = prefix_products(normalized)
    states = []
    for index, (base, product) in enumerate(zip(normalized, products), start=1):
        q = projective_image(product)
        states.append(
            CausalState(
                index=index,
                base=base,
                product=product,
                q=q,
                signed_log_q=math.log(q.numerator) - math.log(q.denominator),
                residual=_residual(q),
                response=exact_forcing_response(q),
            )
        )
    return tuple(states)


def terminal_product(sequence: str) -> Matrix:
    products = prefix_products(sequence)
    return IDENTITY if not products else products[-1]


def conjugate_product_identity(sequence: str) -> bool:
    """Complement-only identity: star reverses multiplication order.

    star(M_n...M_1) = star(M_1)...star(M_n). Therefore complement-only is
    checked against explicit operatorwise conjugation rather than a false
    terminal equality with star(total).
    """
    normalized = normalize_dna(sequence)
    complemented = complement(normalized)
    expected = tuple(reciprocal_conjugate(matrix) for matrix in operator_word(normalized))
    return operator_word(complemented) == expected


def reverse_complement_terminal_identity(sequence: str) -> bool:
    """For reverse complement, terminal(rc(w)) == star(terminal(w))."""
    return terminal_product(reverse_complement(sequence)) == reciprocal_conjugate(
        terminal_product(sequence)
    )


def double_reverse_complement_identity(sequence: str) -> bool:
    normalized = normalize_dna(sequence)
    return reverse_complement(reverse_complement(normalized)) == normalized and operator_word(
        reverse_complement(reverse_complement(normalized))
    ) == operator_word(normalized)


def order_is_noncommutative(left: str, right: str) -> bool:
    if len(left) != len(right):
        raise ValueError("comparison words must have equal length")
    if sorted(left) != sorted(right):
        raise ValueError("comparison words must be composition matched")
    return terminal_product(left) != terminal_product(right)


def _mean(values: Iterable[float]) -> float:
    data = tuple(values)
    return 0.0 if not data else math.fsum(data) / len(data)


def _trajectory_summary(states: Sequence[CausalState], *, include_dhse: bool) -> tuple[float, ...]:
    if not states:
        return (0.0,) * (20 if include_dhse else 12)
    logs = tuple(state.signed_log_q for state in states)
    residuals = tuple(float(state.residual) for state in states)
    log_delta = tuple(b - a for a, b in zip(logs, logs[1:]))
    residual_delta = tuple(b - a for a, b in zip(residuals, residuals[1:]))
    terminal = states[-1].product
    a, b, c, d = terminal
    scale = max(abs(a), abs(b), abs(c), abs(d), 1)
    base = (
        _mean(logs),
        _mean(abs(value) for value in logs),
        logs[-1],
        _mean(residuals),
        max(residuals),
        residuals[-1],
        _mean(log_delta),
        _mean(abs(value) for value in log_delta),
        _mean(residual_delta),
        _mean(abs(value) for value in residual_delta),
        determinant(terminal) / (scale * scale),
        trace(terminal) / scale,
    )
    if not include_dhse:
        return base
    response_means = tuple(
        _mean(state.response[index] for state in states) for index in range(4)
    )
    response_deltas = tuple(
        _mean(
            right.response[index] - left.response[index]
            for left, right in zip(states, states[1:])
        )
        for index in range(4)
    )
    return base + response_means + response_deltas


def causal_dhse_vector(sequence: str, *, windows: Sequence[int] = WINDOWS) -> tuple[float, ...]:
    normalized = normalize_dna(sequence)
    output: list[float] = []
    for length in windows:
        chunks = (
            normalized[start : start + length]
            for start in range(0, max(0, len(normalized) - length + 1))
        )
        summaries = tuple(_trajectory_summary(causal_trajectory(chunk), include_dhse=True) for chunk in chunks)
        width = 20
        if not summaries:
            output.extend((0.0,) * width)
        else:
            output.extend(_mean(row[index] for row in summaries) for index in range(width))
    return tuple(output)


def causal_projective_only_vector(sequence: str, *, windows: Sequence[int] = WINDOWS) -> tuple[float, ...]:
    normalized = normalize_dna(sequence)
    output: list[float] = []
    for length in windows:
        summaries = tuple(
            _trajectory_summary(causal_trajectory(normalized[start : start + length]), include_dhse=False)
            for start in range(0, max(0, len(normalized) - length + 1))
        )
        width = 12
        if not summaries:
            output.extend((0.0,) * width)
        else:
            output.extend(_mean(row[index] for row in summaries) for index in range(width))
    return tuple(output)


def terminal_only_vector(sequence: str, *, windows: Sequence[int] = WINDOWS) -> tuple[float, ...]:
    normalized = normalize_dna(sequence)
    output: list[float] = []
    for length in windows:
        matrices = tuple(
            terminal_product(normalized[start : start + length])
            for start in range(0, max(0, len(normalized) - length + 1))
        )
        if not matrices:
            output.extend((0.0,) * 6)
            continue
        images = tuple(projective_image(matrix) for matrix in matrices)
        output.extend(
            (
                _mean(math.log(q.numerator) - math.log(q.denominator) for q in images),
                _mean(float(_residual(q)) for q in images),
                _mean(determinant(matrix) for matrix in matrices),
                _mean(trace(matrix) for matrix in matrices),
                _mean(matrix[0] / max(matrix) for matrix in matrices),
                _mean(matrix[3] / max(matrix) for matrix in matrices),
            )
        )
    return tuple(output)
