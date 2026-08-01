"""Exact DHSE-001 Stage-M forcing evaluator embedded for Origins-DNA v0.8.

This is a narrow, provenance-preserving port of the primitive K=6 positive
integer Möbius census.  It builds the exact closed-interval step profile once
per process and evaluates rational DNA projective coordinates by binary search.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import gcd
from typing import DefaultDict

import numpy as np

LENGTHS = (1, 2, 3, 4)
SCALE = 6
EXPECTED_PRIMITIVE_MAP_COUNT = 952
EXPECTED = {
    1: {"interval_count": 188496, "breakpoint_count": 44, "q1_count": 24752, "maximum": 39984},
    2: {"interval_count": 1741094, "breakpoint_count": 3359, "q1_count": 314690, "maximum": 314690},
    3: {"interval_count": 5068044, "breakpoint_count": 150239, "q1_count": 943740, "maximum": 943740},
    4: {"interval_count": 11594096, "breakpoint_count": 1595693, "q1_count": 2219236, "maximum": 2224570},
}
EndpointCounts = DefaultDict[tuple[int, int], int]


@dataclass(frozen=True, slots=True)
class ExactForcingProfile:
    """Piecewise-constant closed-interval profile for one word length."""

    length: int
    breakpoints: tuple[Fraction, ...]
    point_counts: tuple[int, ...]
    open_counts: tuple[int, ...]
    interval_count: int
    maximum: int

    def evaluate(self, q: Fraction | int) -> int:
        value = Fraction(q)
        if value <= 0:
            raise ValueError("q must be positive")
        # Local import avoids another module dependency and keeps exact ordering.
        from bisect import bisect_left

        index = bisect_left(self.breakpoints, value)
        if index < len(self.breakpoints) and self.breakpoints[index] == value:
            return self.point_counts[index]
        if index == 0 or index == len(self.breakpoints):
            return 0
        return self.open_counts[index - 1]


def admissible_maps(scale: int = SCALE) -> np.ndarray:
    if scale < 1:
        raise ValueError("scale must be positive")
    rows = []
    for a in range(1, scale + 1):
        for b in range(0, scale + 1):
            for c in range(0, scale + 1):
                for d in range(1, scale + 1):
                    if a * d - b * c > 0:
                        rows.append((a, b, c, d))
    return np.asarray(rows, dtype=np.int64)


def coefficient_gcd(matrix: tuple[int, int, int, int]) -> int:
    value = 0
    for coefficient in matrix:
        value = gcd(value, coefficient)
    return value


@lru_cache(maxsize=1)
def primitive_maps() -> np.ndarray:
    rows = [
        tuple(int(value) for value in row)
        for row in admissible_maps().tolist()
        if coefficient_gcd(tuple(int(value) for value in row)) == 1
    ]
    output = np.asarray(rows, dtype=np.int64)
    if len(output) != EXPECTED_PRIMITIVE_MAP_COUNT:
        raise AssertionError("primitive DHSE map count changed")
    return output


def words_of_length(length: int) -> tuple[str, ...]:
    if length < 1:
        raise ValueError("word length must be positive")
    return tuple("".join(bits) for bits in product("LR", repeat=length))


def _left_multiply(
    matrix: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    current: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    a, b, c, d = matrix
    A, B, C, D = current
    return (
        a * A + b * C,
        a * B + b * D,
        c * A + d * C,
        c * B + d * D,
    )


def _word_matrix(
    word: str,
    left: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    right: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    current = (
        np.ones((1, 1), dtype=np.int64),
        np.zeros((1, 1), dtype=np.int64),
        np.zeros((1, 1), dtype=np.int64),
        np.ones((1, 1), dtype=np.int64),
    )
    for letter in word:
        current = _left_multiply(left if letter == "L" else right, current)
    return tuple(np.broadcast_to(entry, (size, size)) for entry in current)


def _reduced_pairs(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    common = np.gcd(numerator, denominator)
    pairs = np.empty(len(numerator), dtype=[("n", "<i8"), ("d", "<i8")])
    pairs["n"] = numerator // common
    pairs["d"] = denominator // common
    return pairs


def _accumulate(target: EndpointCounts, pairs: np.ndarray) -> None:
    unique, counts = np.unique(pairs, return_counts=True)
    for numerator, denominator, count in zip(unique["n"], unique["d"], counts):
        target[(int(numerator), int(denominator))] += int(count)


def endpoint_census(
    length: int,
) -> tuple[dict[tuple[int, int], int], dict[tuple[int, int], int], int]:
    maps = primitive_maps()
    size = len(maps)
    a, b, c, d = (maps[:, index] for index in range(4))
    left = (a[None, :], b[None, :], c[None, :], d[None, :])
    right = (a[:, None], b[:, None], c[:, None], d[:, None])
    starts: EndpointCounts = defaultdict(int)
    ends: EndpointCounts = defaultdict(int)
    interval_count = 0
    for word in words_of_length(length):
        A, B, C, D = _word_matrix(word, left, right, size)
        valid = (B > 0) & (C > 0) & (81 * A * D <= 121 * B * C)
        A, B, C, D = A[valid], B[valid], C[valid], D[valid]
        lower = _reduced_pairs(9 * A, 11 * C)
        upper = _reduced_pairs(11 * B, 9 * D)
        _accumulate(starts, lower)
        _accumulate(ends, upper)
        interval_count += len(lower)
    return dict(starts), dict(ends), interval_count


def _reciprocal_endpoint_symmetry(
    starts: dict[tuple[int, int], int],
    ends: dict[tuple[int, int], int],
) -> bool:
    return all(
        ends.get((denominator, numerator), 0) == count
        for (numerator, denominator), count in starts.items()
    ) and all(
        starts.get((denominator, numerator), 0) == count
        for (numerator, denominator), count in ends.items()
    )


@lru_cache(maxsize=4)
def exact_profile(length: int) -> ExactForcingProfile:
    if length not in LENGTHS:
        raise ValueError(f"length must be one of {LENGTHS}")
    starts, ends, interval_count = endpoint_census(length)
    breakpoints = tuple(
        sorted(Fraction(numerator, denominator) for numerator, denominator in set(starts) | set(ends))
    )
    active = 0
    point_counts: list[int] = []
    open_counts: list[int] = []
    maximum = 0
    for index, point in enumerate(breakpoints):
        key = (point.numerator, point.denominator)
        point_count = active + starts.get(key, 0)
        point_counts.append(point_count)
        maximum = max(maximum, point_count)
        active = point_count - ends.get(key, 0)
        if index < len(breakpoints) - 1:
            open_counts.append(active)
            maximum = max(maximum, active)
    if active != 0:
        raise AssertionError("DHSE endpoint sweep did not close")
    expected = EXPECTED[length]
    profile = ExactForcingProfile(
        length=length,
        breakpoints=breakpoints,
        point_counts=tuple(point_counts),
        open_counts=tuple(open_counts),
        interval_count=interval_count,
        maximum=maximum,
    )
    if (
        interval_count != expected["interval_count"]
        or len(breakpoints) != expected["breakpoint_count"]
        or profile.evaluate(Fraction(1)) != expected["q1_count"]
        or maximum != expected["maximum"]
        or not _reciprocal_endpoint_symmetry(starts, ends)
    ):
        raise AssertionError(f"DHSE Stage-M certificate mismatch at length {length}")
    return profile


@lru_cache(maxsize=None)
def exact_forcing_response(q: Fraction) -> tuple[int, int, int, int]:
    value = Fraction(q)
    if value <= 0:
        raise ValueError("q must be positive")
    return tuple(exact_profile(length).evaluate(value) for length in LENGTHS)  # type: ignore[return-value]


def exact_stage_m_receipt() -> dict[str, object]:
    rows = []
    for length in LENGTHS:
        profile = exact_profile(length)
        rows.append(
            {
                "length": length,
                "interval_count": profile.interval_count,
                "breakpoint_count": len(profile.breakpoints),
                "q1_count": profile.evaluate(Fraction(1)),
                "maximum": profile.maximum,
                "reciprocal_probe_passed": all(
                    profile.evaluate(q) == profile.evaluate(1 / q)
                    for q in (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1), Fraction(4, 3), Fraction(2), Fraction(4))
                ),
            }
        )
    return {
        "schema": "origins-dna.dhse-stage-m-port.v1",
        "source_repository": "AdrianLipa90/secret-of-a-half",
        "source_stage": "DHSE-001 Stage M",
        "parameters": {"scale": SCALE, "primitive_map_count": EXPECTED_PRIMITIVE_MAP_COUNT, "lengths": list(LENGTHS), "radius": [1, 10]},
        "lengths": rows,
        "technical_status": "PASS",
    }
