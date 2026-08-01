"""Exact DHSE projective-balance representation for DNA windows.

The module maps complementary nucleotide partitions to positive rational odds.
Watson--Crick complementation acts exactly as q -> 1/q on dual channels.  The
DHSE response is evaluated against the frozen Stage-M rational breakpoint
certificate from ``secret-of-a-half``; no fitted parameters or grid
interpolation are used.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import math
from typing import Iterable, Sequence

from .encoding import normalize_dna
from .orbital_atomic import complement, reverse_complement

DEFAULT_WINDOWS = (3, 6, 12, 24)
DEFAULT_STRIDE = 1
PRIMARY_OFFSET = Fraction(1, 1)
DUAL_CHANNELS = ("AT", "CG", "RY", "MK")
INVARIANT_CHANNEL = "WS"
DHSE_LENGTHS = (1, 2, 3, 4)

# Exact maximizer certificate imported from DHSE-001 Stage M.  These constants
# are provenance checks and are not used to fit the DNA representation.
DHSE_STAGE_M_MAXIMA = {
    1: ((Fraction(9, 11), Fraction(11, 12)), (Fraction(12, 11), Fraction(11, 9))),
    2: ((Fraction(1), Fraction(1)),),
    3: ((Fraction(1), Fraction(1)),),
    4: ((Fraction(9882, 9911), Fraction(341, 342)), (Fraction(342, 341), Fraction(9911, 9882))),
}
DHSE_STAGE_M_COUNTS_AT_ONE = {1: 24752, 2: 314690, 3: 943740, 4: 2219236}
DHSE_STAGE_M_MAXIMUM_COUNTS = {1: 39984, 2: 314690, 3: 943740, 4: 2224570}


@dataclass(frozen=True, slots=True)
class ProjectiveWindow:
    start: int
    stop: int
    odds: tuple[tuple[str, Fraction], ...]

    def odds_dict(self) -> dict[str, Fraction]:
        return dict(self.odds)


@dataclass(frozen=True, slots=True)
class DHSEWindowState:
    start: int
    stop: int
    channel: str
    q: Fraction
    signed_log_odds: float
    self_dual_residual: Fraction
    responses: tuple[int, int, int, int]


def _validate_offset(offset: Fraction | int) -> Fraction:
    value = Fraction(offset)
    if value <= 0:
        raise ValueError("offset must be positive")
    return value


def _count_partition(sequence: str, channel: str) -> tuple[int, int]:
    counts = {base: sequence.count(base) for base in "ACGT"}
    if channel == "AT":
        return counts["A"], counts["T"]
    if channel == "CG":
        return counts["C"], counts["G"]
    if channel == "RY":
        return counts["A"] + counts["G"], counts["C"] + counts["T"]
    if channel == "MK":
        return counts["A"] + counts["C"], counts["G"] + counts["T"]
    if channel == "WS":
        return counts["A"] + counts["T"], counts["C"] + counts["G"]
    raise ValueError(f"unknown projective channel: {channel}")


def projective_odds(
    sequence: str,
    channel: str,
    *,
    offset: Fraction | int = PRIMARY_OFFSET,
) -> Fraction:
    """Return an exact positive rational odds coordinate for one DNA sequence."""

    normalized = normalize_dna(sequence)
    alpha = _validate_offset(offset)
    numerator, denominator = _count_partition(normalized, channel)
    return (Fraction(numerator) + alpha) / (Fraction(denominator) + alpha)


def projective_window_odds(
    sequence: str,
    window_size: int,
    *,
    stride: int = DEFAULT_STRIDE,
    offset: Fraction | int = PRIMARY_OFFSET,
) -> tuple[ProjectiveWindow, ...]:
    if not isinstance(window_size, int) or not isinstance(stride, int):
        raise TypeError("window_size and stride must be integers")
    if window_size <= 0 or stride <= 0:
        raise ValueError("window_size and stride must be positive")
    normalized = normalize_dna(sequence)
    if window_size > len(normalized):
        return ()
    channels = DUAL_CHANNELS + (INVARIANT_CHANNEL,)
    output = []
    for start in range(0, len(normalized) - window_size + 1, stride):
        chunk = normalized[start : start + window_size]
        output.append(
            ProjectiveWindow(
                start=start,
                stop=start + window_size,
                odds=tuple(
                    (channel, projective_odds(chunk, channel, offset=offset))
                    for channel in channels
                ),
            )
        )
    return tuple(output)


def self_dual_residual(q: Fraction | int) -> Fraction:
    value = Fraction(q)
    if value <= 0:
        raise ValueError("q must be positive")
    return abs(value - 1) / (value + 1)


def signed_log_odds(q: Fraction | int) -> float:
    value = Fraction(q)
    if value <= 0:
        raise ValueError("q must be positive")
    return math.log(value.numerator) - math.log(value.denominator)


def _closed_interval_contains(q: Fraction, interval: tuple[Fraction, Fraction]) -> bool:
    return interval[0] <= q <= interval[1]


def dhse_certificate_response(q: Fraction | int) -> tuple[int, int, int, int]:
    """Evaluate the frozen Stage-M maximum-membership certificate.

    Each component is the exact Stage-M maximum count when ``q`` belongs to a
    certified maximizer component and the certified q=1 count otherwise.  This
    compact response is deliberately certificate-based.  A later heavyweight
    evaluator may import the complete breakpoint census without changing this
    preregistered public contract.
    """

    value = Fraction(q)
    if value <= 0:
        raise ValueError("q must be positive")
    output = []
    for length in DHSE_LENGTHS:
        on_maximum = any(
            _closed_interval_contains(value, interval)
            for interval in DHSE_STAGE_M_MAXIMA[length]
        )
        output.append(
            DHSE_STAGE_M_MAXIMUM_COUNTS[length]
            if on_maximum
            else DHSE_STAGE_M_COUNTS_AT_ONE[length]
        )
    return tuple(output)  # type: ignore[return-value]


def dhse_window_states(
    sequence: str,
    window_size: int,
    *,
    stride: int = DEFAULT_STRIDE,
    offset: Fraction | int = PRIMARY_OFFSET,
) -> tuple[DHSEWindowState, ...]:
    output = []
    for window in projective_window_odds(
        sequence, window_size, stride=stride, offset=offset
    ):
        for channel, q in window.odds:
            output.append(
                DHSEWindowState(
                    start=window.start,
                    stop=window.stop,
                    channel=channel,
                    q=q,
                    signed_log_odds=signed_log_odds(q),
                    self_dual_residual=self_dual_residual(q),
                    responses=dhse_certificate_response(q),
                )
            )
    return tuple(output)


def _mean(values: Iterable[float]) -> float:
    materialized = tuple(values)
    return 0.0 if not materialized else math.fsum(materialized) / len(materialized)


def _normalized_response(value: int, length: int) -> float:
    return value / DHSE_STAGE_M_MAXIMUM_COUNTS[length]


def _channel_summary(states: Sequence[DHSEWindowState], channel: str) -> tuple[float, ...]:
    selected = tuple(state for state in states if state.channel == channel)
    if not selected:
        return (0.0,) * 12
    logs = tuple(state.signed_log_odds for state in selected)
    residuals = tuple(float(state.self_dual_residual) for state in selected)
    response_means = tuple(
        _mean(
            _normalized_response(state.responses[index], length)
            for state in selected
        )
        for index, length in enumerate(DHSE_LENGTHS)
    )
    log_differences = tuple(right - left for left, right in zip(logs, logs[1:]))
    residual_differences = tuple(
        right - left for left, right in zip(residuals, residuals[1:])
    )
    return (
        _mean(logs),
        _mean(abs(value) for value in logs),
        _mean(residuals),
        max(residuals, default=0.0),
        *response_means,
        _mean(log_differences),
        _mean(abs(value) for value in log_differences),
        _mean(residual_differences),
        _mean(abs(value) for value in residual_differences),
    )


def dhse_projective_vector(
    sequence: str,
    *,
    windows: Sequence[int] = DEFAULT_WINDOWS,
    stride: int = DEFAULT_STRIDE,
    offset: Fraction | int = PRIMARY_OFFSET,
    include_signed_log: bool = True,
    include_dhse: bool = True,
    channels: Sequence[str] = DUAL_CHANNELS,
) -> tuple[float, ...]:
    """Return a deterministic fixed-size multiscale DHSE DNA representation."""

    normalized = normalize_dna(sequence)
    materialized_windows = tuple(windows)
    materialized_channels = tuple(channels)
    if not materialized_windows or any(value <= 0 for value in materialized_windows):
        raise ValueError("windows must contain positive integers")
    if not materialized_channels:
        raise ValueError("channels must not be empty")
    output: list[float] = []
    for window_size in materialized_windows:
        states = dhse_window_states(
            normalized, window_size, stride=stride, offset=offset
        )
        for channel in materialized_channels:
            summary = _channel_summary(states, channel)
            # summary layout: log(2), residual(2), DHSE(4), local changes(4)
            if include_signed_log and include_dhse:
                output.extend(summary)
            elif include_signed_log:
                output.extend(summary[:4] + summary[8:])
            elif include_dhse:
                output.extend(summary[2:8] + summary[10:])
            else:
                output.extend(summary[2:4] + summary[10:])
    return tuple(output)


def projective_only_vector(sequence: str) -> tuple[float, ...]:
    return dhse_projective_vector(sequence, include_dhse=False)


def dhse_without_signed_log_vector(sequence: str) -> tuple[float, ...]:
    return dhse_projective_vector(sequence, include_signed_log=False)


def invariant_ws_control_vector(sequence: str) -> tuple[float, ...]:
    return dhse_projective_vector(sequence, channels=(INVARIANT_CHANNEL,))


def verify_complement_covariance(sequence: str) -> dict[str, bool]:
    normalized = normalize_dna(sequence)
    complemented = complement(normalized)
    dual_inverse = all(
        projective_odds(complemented, channel)
        == 1 / projective_odds(normalized, channel)
        for channel in DUAL_CHANNELS
    )
    ws_invariant = projective_odds(complemented, INVARIANT_CHANNEL) == projective_odds(
        normalized, INVARIANT_CHANNEL
    )
    residual_invariant = all(
        self_dual_residual(projective_odds(complemented, channel))
        == self_dual_residual(projective_odds(normalized, channel))
        for channel in DUAL_CHANNELS
    )
    response_invariant = all(
        dhse_certificate_response(projective_odds(complemented, channel))
        == dhse_certificate_response(projective_odds(normalized, channel))
        for channel in DUAL_CHANNELS
    )
    signed_log_antisymmetric = all(
        math.isclose(
            signed_log_odds(projective_odds(complemented, channel)),
            -signed_log_odds(projective_odds(normalized, channel)),
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for channel in DUAL_CHANNELS
    )
    return {
        "dual_inverse": dual_inverse,
        "ws_invariant": ws_invariant,
        "residual_invariant": residual_invariant,
        "response_invariant": response_invariant,
        "signed_log_antisymmetric": signed_log_antisymmetric,
    }


def reverse_complement_window_covariance(
    sequence: str,
    window_size: int,
) -> bool:
    forward = projective_window_odds(sequence, window_size)
    reversed_states = projective_window_odds(reverse_complement(sequence), window_size)
    expected = tuple(
        tuple(
            (channel, value if channel == INVARIANT_CHANNEL else 1 / value)
            for channel, value in window.odds
        )
        for window in reversed(forward)
    )
    observed = tuple(window.odds for window in reversed_states)
    return observed == expected
