from __future__ import annotations

import math

import pytest

from origins_dna.benchmark_types import SyntheticMotifConfig
from origins_dna.causal_three_quarter import (
    THREE_QUARTER,
    even_three_quarter_vector,
    even_window_features,
    exact_three_quarter_gate_receipt,
    phase_ordinary_mean_vector,
    phase_three_quarter_vector,
    phase_window_features,
    self_dual_response,
    signed_three_quarter_amplitude,
    three_quarter_aggregate,
)
from origins_dna.causal_three_quarter_benchmark import (
    generate_paired_motif_dataset,
    paired_split_indices,
)


def test_three_quarter_scaling_is_exactly_frozen() -> None:
    assert THREE_QUARTER == 0.75
    assert three_quarter_aggregate([1.0] * 16) == pytest.approx(2.0)
    assert three_quarter_aggregate([1.0] * 16, exponent=1.0) == pytest.approx(1.0)
    assert signed_three_quarter_amplitude(16.0) == pytest.approx(8.0)
    assert signed_three_quarter_amplitude(-16.0) == pytest.approx(-8.0)


def test_exact_motif_decoy_reciprocal_gate() -> None:
    receipt = exact_three_quarter_gate_receipt()
    assert receipt["technical_status"] == "PASS"
    assert receipt["q_motif"] == [8, 5]
    assert receipt["q_decoy"] == [5, 8]
    assert receipt["reciprocal_terminal_coordinates"]
    assert receipt["opposite_nonzero_signed_phase"]
    assert receipt["equal_stage_m_response"]


def test_phase_is_odd_and_even_geometry_is_equal_for_frozen_pair() -> None:
    left = phase_window_features("CGCG")
    right = phase_window_features("GCGC")
    assert left[0] == pytest.approx(-right[0])
    assert left[1] == pytest.approx(-right[1])
    assert left[2] == pytest.approx(-right[2])
    assert left[3] == pytest.approx(-right[3])
    assert left[4] == pytest.approx(-right[5])
    assert left[5] == pytest.approx(-right[4])
    assert even_window_features("CGCG") == pytest.approx(even_window_features("GCGC"))


def test_self_dual_response_is_dimensionless_and_centered() -> None:
    references = (10, 20, 30, 40)
    assert self_dual_response(references, references=references) == (0.0,) * 4
    transformed = self_dual_response((20, 10, 30, 80), references=references)
    assert transformed[0] > 0.0
    assert transformed[1] < 0.0
    assert transformed[2] == 0.0
    assert transformed[3] > 0.0


def test_multiscale_vectors_are_finite_fixed_and_order_sensitive() -> None:
    left_sequence = "ACGT" * 8 + "CGCG"
    right_sequence = "ACGT" * 8 + "GCGC"
    phase_left = phase_three_quarter_vector(left_sequence)
    phase_right = phase_three_quarter_vector(right_sequence)
    mean_left = phase_ordinary_mean_vector(left_sequence)
    even_left = even_three_quarter_vector(left_sequence)
    assert len(phase_left) == len(phase_right) == len(mean_left) == 144
    assert len(even_left) == 408
    assert phase_left != phase_right
    assert phase_left != mean_left
    assert all(math.isfinite(value) for value in phase_left + phase_right + even_left)


def test_paired_background_generation_and_group_split_have_no_leakage() -> None:
    config = SyntheticMotifConfig(
        sample_count=20,
        sequence_length=32,
        motif="CGCG",
        seeds=(11,),
        bootstrap_resamples=10,
    )
    dataset = generate_paired_motif_dataset(config, seed=11)
    assert len(dataset.sequences) == 20
    assert len(set(dataset.pair_ids)) == 10
    for pair_id in set(dataset.pair_ids):
        indices = [i for i, value in enumerate(dataset.pair_ids) if value == pair_id]
        positive = next(dataset.sequences[i] for i in indices if dataset.labels[i] == 1)
        negative = next(dataset.sequences[i] for i in indices if dataset.labels[i] == 0)
        start_positive = positive.index("CGCG")
        start_negative = negative.index("GCGC")
        assert start_positive == start_negative
        assert positive[:start_positive] == negative[:start_negative]
        assert positive[start_positive + 4 :] == negative[start_negative + 4 :]
        assert "GCGC" not in positive
        assert "CGCG" not in negative
    train, test = paired_split_indices(
        dataset.pair_ids, train_fraction=0.7, seed=37
    )
    train_pairs = {dataset.pair_ids[index] for index in train}
    test_pairs = {dataset.pair_ids[index] for index in test}
    assert train_pairs.isdisjoint(test_pairs)
    assert {dataset.labels[index] for index in train} == {0, 1}
    assert {dataset.labels[index] for index in test} == {0, 1}
