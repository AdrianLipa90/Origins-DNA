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
    """Return J @ matrix @ J for J=((0, 1), (1, 0))."""
    a, b, c, d = matrix
    return (d, c, b, a)


def multiply(left: Matrix, right: Matrix) -> Matrix:
    a, b, c, d = left
    A, B, C, D = right
    return (
        a * A + b * C,
        a * B + b * D,
        c * A + d * C,
        c * B + d * D,
    )


def reciprocal_conjugation_identity(left: Matrix, right: Matrix) -> bool:
    """Check star(left @ right) == star(left) @ star(right) exactly."""
    return reciprocal_conjugate(multiply(left, right)) == multiply(
        reciprocal_conjugate(left), reciprocal_conjugate(right)
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


def _terminal_from_word(word: str, assignment: dict[str, Matrix]) -> Matrix:
    current = IDENTITY
    for base in word:
        current = multiply(assignment[base], current)
    return current


@lru_cache(maxsize=1)
def base_operator_assignment() -> dict[str, Matrix]:
    """Return the amended deterministic reciprocal-pair assignment."""
    matrices = sorted(
        tuple(int(value) for value in row) for row in primitive_maps().tolist()
    )
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
    first = pairs[0]
    chosen: tuple[Matrix, Matrix] | None = None
    for candidate in pairs[1:]:
        trial = {
            "A": first[0],
            "C": candidate[0],
            "G": candidate[1],
            "T": first[1],
        }
        if len(set(trial.values())) != 4:
            continue
        if _terminal_from_word("CGCG", trial) != _terminal_from_word("CCGG", trial):
            chosen = candidate
            break
    if chosen is None:
        raise RuntimeError("no admissible noncommutative C/G reciprocal pair found")
    assignment = {
        "A": first[0],
        "C": chosen[0],
        "G": chosen[1],
        "T": first[1],
    }
    if tuple(assignment) != tuple(BASES):
        raise RuntimeError("causal operator assignment is not in canonical base order")
    return assignment


def operator_word(sequence: str) -> tuple[Matrix, ...]:
    normalized = normalize_dna(sequence)
    assignment = base_operator_assignment()
    return tuple(assignment[base] for base in normalized)


def prefix_products(sequence: str) -> tuple[Matrix, ...]:
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
                index,
                base,
                product,
                q,
                math.log(q.numerator) - math.log(q.denominator),
                _residual(q),
                exact_forcing_response(q),
            )
        )
    return tuple(states)


def terminal_product(sequence: str) -> Matrix:
    products = prefix_products(sequence)
    return IDENTITY if not products else products[-1]


def conjugate_product_identity(sequence: str) -> bool:
    """Check complement at both operator-word and terminal-product levels."""
    normalized = normalize_dna(sequence)
    complemented = complement(normalized)
    expected_word = tuple(
        reciprocal_conjugate(matrix) for matrix in operator_word(normalized)
    )
    word_identity = operator_word(complemented) == expected_word
    terminal_identity = terminal_product(complemented) == reciprocal_conjugate(
        terminal_product(normalized)
    )
    return word_identity and terminal_identity


def reverse_complement_expected_terminal(sequence: str) -> Matrix:
    """Return the exact terminal implied by reverse-complement order.

    The involution is conjugation by J and therefore an automorphism:
    star(A @ B) = star(A) @ star(B).

    For w=x1...xn, terminal(w)=Mn...M1. Reverse-complement changes the word to
    comp(xn)...comp(x1), hence its terminal is star(M1)...star(Mn). This is
    generally different from star(terminal(w)).
    """
    current = IDENTITY
    for operator in operator_word(sequence):
        current = multiply(current, reciprocal_conjugate(operator))
    return current


def reverse_complement_terminal_identity(sequence: str) -> bool:
    return terminal_product(
        reverse_complement(sequence)
    ) == reverse_complement_expected_terminal(sequence)


def double_reverse_complement_identity(sequence: str) -> bool:
    normalized = normalize_dna(sequence)
    twice = reverse_complement(reverse_complement(normalized))
    return twice == normalized and operator_word(twice) == operator_word(normalized)


def order_is_noncommutative(left: str, right: str) -> bool:
    if len(left) != len(right):
        raise ValueError("comparison words must have equal length")
    if sorted(left) != sorted(right):
        raise ValueError("comparison words must be composition matched")
    return terminal_product(left) != terminal_product(right)


def _mean(values: Iterable[float]) -> float:
    data = tuple(values)
    return 0.0 if not data else math.fsum(data) / len(data)


def _trajectory_summary(
    states: Sequence[CausalState], *, include_dhse: bool
) -> tuple[float, ...]:
    if not states:
        return (0.0,) * (20 if include_dhse else 12)
    logs = tuple(state.signed_log_q for state in states)
    residuals = tuple(float(state.residual) for state in states)
    log_delta = tuple(b - a for a, b in zip(logs, logs[1:]))
    residual_delta = tuple(b - a for a, b in zip(residuals, residuals[1:]))
    terminal = states[-1].product
    scale = max(abs(value) for value in terminal) or 1
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


def causal_dhse_vector(
    sequence: str, *, windows: Sequence[int] = WINDOWS
) -> tuple[float, ...]:
    normalized = normalize_dna(sequence)
    output: list[float] = []
    for length in windows:
        summaries = tuple(
            _trajectory_summary(
                causal_trajectory(normalized[start : start + length]),
                include_dhse=True,
            )
            for start in range(0, max(0, len(normalized) - length + 1))
        )
        output.extend(
            (0.0,) * 20
            if not summaries
            else (_mean(row[index] for row in summaries) for index in range(20))
        )
    return tuple(output)


def causal_projective_only_vector(
    sequence: str, *, windows: Sequence[int] = WINDOWS
) -> tuple[float, ...]:
    normalized = normalize_dna(sequence)
    output: list[float] = []
    for length in windows:
        summaries = tuple(
            _trajectory_summary(
                causal_trajectory(normalized[start : start + length]),
                include_dhse=False,
            )
            for start in range(0, max(0, len(normalized) - length + 1))
        )
        output.extend(
            (0.0,) * 12
            if not summaries
            else (_mean(row[index] for row in summaries) for index in range(12))
        )
    return tuple(output)


def terminal_only_vector(
    sequence: str, *, windows: Sequence[int] = WINDOWS
) -> tuple[float, ...]:
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
                _mean(
                    math.log(q.numerator) - math.log(q.denominator) for q in images
                ),
                _mean(float(_residual(q)) for q in images),
                _mean(determinant(matrix) for matrix in matrices),
                _mean(trace(matrix) for matrix in matrices),
                _mean(matrix[0] / max(matrix) for matrix in matrices),
                _mean(matrix[3] / max(matrix) for matrix in matrices),
            )
        )
    return tuple(output)
