import math

import pytest

from origins_dna.orbital_atomic import (
    complement,
    consolidate_windows,
    orbital_atomic_encode,
    orbital_atomic_summary,
    reverse_complement,
)


def test_complement_and_reverse_complement() -> None:
    assert complement("ACGT") == "TGCA"
    assert reverse_complement("AACG") == "CGTT"


def test_complement_is_involution() -> None:
    sequence = "ACGTTGCA"
    assert complement(complement(sequence)) == sequence


def test_antipodal_complement_symmetry() -> None:
    original = orbital_atomic_encode("ACGT")
    paired = orbital_atomic_encode(complement("ACGT"))

    for left, right in zip(original, paired):
        assert right.phase_cos == pytest.approx(-left.phase_cos, abs=1e-12)
        assert right.phase_sin == pytest.approx(-left.phase_sin, abs=1e-12)
        assert right.ring_class == -left.ring_class
        assert right.bond_class == left.bond_class


def test_unit_circle_phase_channel() -> None:
    for state in orbital_atomic_encode("ACGT"):
        assert math.hypot(state.phase_cos, state.phase_sin) == pytest.approx(1.0)


def test_state_tuple_order_is_stable() -> None:
    state = orbital_atomic_encode("A")[0]
    assert state.as_tuple() == (1.0, 0.0, 1.0, -1.0)


def test_window_consolidation_coordinates_and_count() -> None:
    windows = consolidate_windows("ACGTAC", 3, stride=2)
    assert [(window.start, window.stop) for window in windows] == [(0, 3), (2, 5)]
    assert all(0.0 <= window.phase_coherence <= 1.0 for window in windows)


def test_window_larger_than_sequence_is_empty() -> None:
    assert consolidate_windows("ACG", 4) == ()


def test_window_validation() -> None:
    with pytest.raises(ValueError):
        consolidate_windows("ACGT", 0)
    with pytest.raises(ValueError):
        consolidate_windows("ACGT", 2, stride=0)
    with pytest.raises(TypeError):
        consolidate_windows("ACGT", 2.0)


def test_balanced_sequence_summary() -> None:
    summary = orbital_atomic_summary("ACGT")
    assert summary["mean_phase_cos"] == pytest.approx(0.0, abs=1e-12)
    assert summary["mean_phase_sin"] == pytest.approx(0.0, abs=1e-12)
    assert summary["phase_coherence"] == pytest.approx(0.0, abs=1e-12)
    assert summary["mean_ring_class"] == pytest.approx(0.0)
    assert summary["mean_bond_class"] == pytest.approx(0.0)


def test_empty_summary_is_zero_vector() -> None:
    assert set(orbital_atomic_summary("").values()) == {0.0}
