"""Experimental orbital/atomic DNA representation with explicit symmetries.

This module defines a mathematical feature map. It does not claim that the map
is a biological mechanism or that DNA has been decrypted.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .encoding import normalize_dna

_PHASE_BY_BASE = {
    "A": 0.0,
    "C": math.pi / 2.0,
    "T": math.pi,
    "G": 3.0 * math.pi / 2.0,
}
_RING_CLASS = {"A": 1.0, "G": 1.0, "C": -1.0, "T": -1.0}
_BOND_CLASS = {"A": -1.0, "T": -1.0, "C": 1.0, "G": 1.0}
_COMPLEMENT_TABLE = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True, slots=True)
class OrbitalAtomicState:
    """One nucleotide mapped to phase, ring class and bond-strength channels."""

    base: str
    phase_cos: float
    phase_sin: float
    ring_class: float
    bond_class: float

    def as_tuple(self) -> tuple[float, float, float, float]:
        """Return the numeric channels in stable public order."""

        return (self.phase_cos, self.phase_sin, self.ring_class, self.bond_class)


@dataclass(frozen=True, slots=True)
class ConsolidatedWindow:
    """Deterministic local consolidation over a fixed sequence window."""

    start: int
    stop: int
    mean_phase_cos: float
    mean_phase_sin: float
    phase_coherence: float
    mean_ring_class: float
    mean_bond_class: float
    signed_phase_flow: float


def complement(sequence: str) -> str:
    """Return the Watson--Crick complement of a normalized sequence."""

    return normalize_dna(sequence).translate(_COMPLEMENT_TABLE)


def reverse_complement(sequence: str) -> str:
    """Return the reverse Watson--Crick complement."""

    return complement(sequence)[::-1]


def orbital_atomic_encode(sequence: str) -> tuple[OrbitalAtomicState, ...]:
    """Map A/C/G/T to a four-channel experimental representation.

    Complementary bases occupy antipodal points on the unit circle.
    ``ring_class`` is +1 for purines and -1 for pyrimidines.
    ``bond_class`` is -1 for A/T and +1 for C/G.
    """

    normalized = normalize_dna(sequence)
    return tuple(_state_for_base(base) for base in normalized)


def _state_for_base(base: str) -> OrbitalAtomicState:
    phase = _PHASE_BY_BASE[base]
    return OrbitalAtomicState(
        base=base,
        phase_cos=math.cos(phase),
        phase_sin=math.sin(phase),
        ring_class=_RING_CLASS[base],
        bond_class=_BOND_CLASS[base],
    )


def _wrapped_phase_delta(left: float, right: float) -> float:
    """Return the shortest signed angular displacement from left to right.

    The exactly-antipodal tie is represented as +pi, giving a stable convention
    for transitions whose shortest direction is otherwise mathematically ambiguous.
    """

    delta = (right - left + math.pi) % (2.0 * math.pi) - math.pi
    if math.isclose(delta, -math.pi, rel_tol=0.0, abs_tol=1e-12):
        return math.pi
    return delta


def consolidate_windows(
    sequence: str,
    window_size: int,
    *,
    stride: int = 1,
) -> tuple[ConsolidatedWindow, ...]:
    """Consolidate local orbital/atomic states into overlapping windows."""

    _validate_window_parameters(window_size, stride)
    normalized = normalize_dna(sequence)
    if window_size > len(normalized):
        return ()

    states = orbital_atomic_encode(normalized)
    return tuple(
        _consolidate_chunk(states, start, start + window_size)
        for start in range(0, len(states) - window_size + 1, stride)
    )


def _validate_window_parameters(window_size: int, stride: int) -> None:
    if not isinstance(window_size, int):
        raise TypeError("window_size must be an integer")
    if not isinstance(stride, int):
        raise TypeError("stride must be an integer")
    if window_size <= 0:
        raise ValueError("window_size must be positive")
    if stride <= 0:
        raise ValueError("stride must be positive")


def _consolidate_chunk(
    states: tuple[OrbitalAtomicState, ...],
    start: int,
    stop: int,
) -> ConsolidatedWindow:
    chunk = states[start:stop]
    size = len(chunk)
    mean_cos = math.fsum(state.phase_cos for state in chunk) / size
    mean_sin = math.fsum(state.phase_sin for state in chunk) / size
    phases = [_PHASE_BY_BASE[state.base] for state in chunk]

    return ConsolidatedWindow(
        start=start,
        stop=stop,
        mean_phase_cos=mean_cos,
        mean_phase_sin=mean_sin,
        phase_coherence=math.hypot(mean_cos, mean_sin),
        mean_ring_class=math.fsum(state.ring_class for state in chunk) / size,
        mean_bond_class=math.fsum(state.bond_class for state in chunk) / size,
        signed_phase_flow=_phase_flow(phases),
    )


def _phase_flow(phases: list[float]) -> float:
    return math.fsum(
        _wrapped_phase_delta(left, right)
        for left, right in zip(phases, phases[1:])
    )


def orbital_atomic_summary(sequence: str) -> dict[str, float]:
    """Return a fixed-size summary for baseline model comparisons."""

    states = orbital_atomic_encode(sequence)
    if not states:
        return _zero_summary()

    count = len(states)
    mean_cos = math.fsum(state.phase_cos for state in states) / count
    mean_sin = math.fsum(state.phase_sin for state in states) / count
    phases = [_PHASE_BY_BASE[state.base] for state in states]

    return {
        "mean_phase_cos": mean_cos,
        "mean_phase_sin": mean_sin,
        "phase_coherence": math.hypot(mean_cos, mean_sin),
        "mean_ring_class": math.fsum(state.ring_class for state in states) / count,
        "mean_bond_class": math.fsum(state.bond_class for state in states) / count,
        "signed_phase_flow_per_step": _phase_flow(phases) / max(1, count - 1),
    }


def _zero_summary() -> dict[str, float]:
    return {
        "mean_phase_cos": 0.0,
        "mean_phase_sin": 0.0,
        "phase_coherence": 0.0,
        "mean_ring_class": 0.0,
        "mean_bond_class": 0.0,
        "signed_phase_flow_per_step": 0.0,
    }
