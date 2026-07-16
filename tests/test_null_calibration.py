import json
from pathlib import Path
import subprocess
import sys

import pytest

from origins_dna.benchmark_types import SyntheticMotifConfig
from origins_dna.null_calibration import (
    NullCalibratedAggregate,
    build_null_calibrated_decisions,
    calibrate_projection_sensitivity_report,
    cluster_sign_flip_p_value,
    null_calibration_data_seed_schedule,
    null_calibration_projection_seed_schedule,
    prior_projection_seeds,
    run_null_calibrated_analysis,
    sequence_shuffle_composition_diagnostics,
)
from origins_dna.projection_sensitivity import (
    EXPLORATORY_PILOT_PROJECTION_SEEDS,
    ProjectionSeedResult,
    ProjectionSensitivityReport,
    projection_seed_schedule,
)


def _source_report() -> ProjectionSensitivityReport:
    config = SyntheticMotifConfig(
        sample_count=40,
        sequence_length=16,
        motif="CGCG",
        seeds=(3, 5),
        bootstrap_resamples=50,
    )
    rows = []
    for projection_seed, offset in ((101, 0.0), (202, 0.02), (303, 0.04)):
        rows.extend(
            [
                ProjectionSeedResult(
                    projection_seed,
                    8,
                    "observed_labels",
                    0.25 + offset,
                    (0.20 + offset, 0.30 + offset),
                ),
                ProjectionSeedResult(
                    projection_seed,
                    8,
                    "label_permutation_null",
                    0.075,
                    (0.05, 0.10),
                ),
                ProjectionSeedResult(
                    projection_seed,
                    8,
                    "sequence_shuffle_null",
                    0.10,
                    (0.10, 0.10),
                ),
            ]
        )
    return ProjectionSensitivityReport(
        schema="origins-dna.projection-seed-sensitivity.v1",
        benchmark_status="synthetic_projection_seed_sensitivity_result",
        config=config,
        budgets=(8,),
        projection_master_seed=80801,
        projection_seeds=(101, 202, 303),
        excluded_exploratory_projection_seeds=EXPLORATORY_PILOT_PROJECTION_SEEDS,
        hierarchical_bootstrap_resamples=200,
        seed_results=tuple(rows),
        aggregate_results=(),
        interpretation_boundary="test fixture",
    )


def test_independent_schedules_are_deterministic_and_disjoint() -> None:
    projection_first = null_calibration_projection_seed_schedule()
    projection_second = null_calibration_projection_seed_schedule()
    data_first = null_calibration_data_seed_schedule()
    data_second = null_calibration_data_seed_schedule()
    assert projection_first == projection_second
    assert data_first == data_second
    assert len(projection_first) == len(set(projection_first)) == 32
    assert len(data_first) == len(set(data_first)) == 8
    assert set(projection_first).isdisjoint(prior_projection_seeds())
    assert set(data_first).isdisjoint({11, 23, 37, 41, 53, 67, 79, 83})
    assert set(prior_projection_seeds()).issuperset(
        projection_seed_schedule(master_seed=70_701, count=32)
    )


def test_calibrated_arithmetic_preserves_nulls_separately() -> None:
    report = calibrate_projection_sensitivity_report(
        _source_report(),
        cluster_sign_flip_resamples=500,
    )
    assert report.schema == "origins-dna.null-calibrated-contrast.v1"
    assert len(report.seed_results) == 3 * 2
    first_permutation = next(
        item
        for item in report.seed_results
        if item.projection_seed == 101
        and item.null_condition == "label_permutation_null"
    )
    first_shuffle = next(
        item
        for item in report.seed_results
        if item.projection_seed == 101
        and item.null_condition == "sequence_shuffle_null"
    )
    assert first_permutation.per_data_seed_calibrated_differences == pytest.approx(
        (0.15, 0.20)
    )
    assert first_shuffle.per_data_seed_calibrated_differences == pytest.approx(
        (0.10, 0.20)
    )
    assert {item.null_condition for item in report.aggregate_results} == {
        "label_permutation_null",
        "sequence_shuffle_null",
    }


def test_cluster_sign_flip_is_deterministic_and_validated() -> None:
    values = (0.1, 0.2, 0.3, 0.4)
    first = cluster_sign_flip_p_value(values, resamples=2_000, seed=17)
    second = cluster_sign_flip_p_value(values, resamples=2_000, seed=17)
    assert first == second
    assert 0.0 < first <= 1.0
    assert cluster_sign_flip_p_value((0.0, 0.0), resamples=10, seed=1) == 1.0
    with pytest.raises(ValueError, match="must not be empty"):
        cluster_sign_flip_p_value((), resamples=10, seed=1)


def _aggregate(null_condition: str, direction: str, passed: bool) -> NullCalibratedAggregate:
    positive = direction == "tetranucleotide_advantage"
    return NullCalibratedAggregate(
        budget=8,
        null_condition=null_condition,
        projection_seed_count=32,
        data_seed_count=8,
        mean_calibrated_difference=0.1 if positive else -0.1,
        median_projection_mean_difference=0.1 if positive else -0.1,
        hierarchical_ci_low=0.05 if positive else -0.15,
        hierarchical_ci_high=0.15 if positive else -0.05,
        fraction_projection_means_positive=0.9 if positive else 0.1,
        fraction_projection_means_negative=0.1 if positive else 0.9,
        minimum_projection_mean_difference=-0.01,
        maximum_projection_mean_difference=0.2,
        cluster_sign_flip_p_value=0.01,
        direction=direction,
        interval_gate_passed=passed,
        sign_stability_gate_passed=passed,
        randomization_gate_passed=passed,
        calibrated_gate_passed=passed,
    )


def test_decision_requires_both_nulls_and_same_direction() -> None:
    passed = build_null_calibrated_decisions(
        (
            _aggregate("label_permutation_null", "tetranucleotide_advantage", True),
            _aggregate("sequence_shuffle_null", "tetranucleotide_advantage", True),
        ),
        budgets=(8,),
    )
    assert passed[0].promotion_gate_passed
    disagreed = build_null_calibrated_decisions(
        (
            _aggregate("label_permutation_null", "tetranucleotide_advantage", True),
            _aggregate("sequence_shuffle_null", "orbital_advantage", True),
        ),
        budgets=(8,),
    )
    assert not disagreed[0].promotion_gate_passed
    assert disagreed[0].consensus_direction == "unresolved"


def test_shuffle_composition_diagnostic_is_exact() -> None:
    diagnostics = sequence_shuffle_composition_diagnostics(
        SyntheticMotifConfig(
            sample_count=20,
            sequence_length=16,
            motif="CGCG",
            seeds=(3, 5),
            bootstrap_resamples=20,
        )
    )
    assert diagnostics.maximum_per_sequence_base_count_residual_after_shuffle == 0
    assert diagnostics.shuffle_composition_gate_passed
    assert "unresolved_design_limitation" in diagnostics.rejection_sampling_status


def test_small_end_to_end_analysis_and_bundle(tmp_path: Path) -> None:
    report = run_null_calibrated_analysis(
        sample_count=40,
        sequence_length=16,
        motif="CGCG",
        budgets=(8,),
        projection_master_seed=80801,
        projection_seed_count=2,
        data_master_seed=90901,
        data_seed_count=2,
        bootstrap_resamples=30,
        hierarchical_bootstrap_resamples=50,
        cluster_sign_flip_resamples=100,
    )
    assert len(report.seed_results) == 2 * 2
    assert len(report.aggregate_results) == 2
    assert len(report.decisions) == 1
    assert not report.biological_progression_allowed
    first = report.write_bundle(tmp_path / "result.json")
    second = report.write_bundle(tmp_path / "second" / "result.json")
    for left, right in zip(first, second):
        assert left.read_bytes() == right.read_bytes()
    summary = json.loads(first[2].read_text(encoding="utf-8"))
    assert "seed_results" not in summary


def test_repository_script_runs_without_editable_install(tmp_path: Path) -> None:
    repository_root = Path(__file__).resolve().parents[1]
    output = tmp_path / "result.json"
    completed = subprocess.run(
        [
            sys.executable,
            str(repository_root / "scripts" / "run_null_calibrated_contrast.py"),
            "--output",
            str(output),
            "--samples",
            "40",
            "--length",
            "16",
            "--budgets",
            "8",
            "--projection-seed-count",
            "2",
            "--data-seed-count",
            "2",
            "--bootstrap-resamples",
            "30",
            "--hierarchical-bootstrap-resamples",
            "50",
            "--cluster-sign-flip-resamples",
            "100",
        ],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert output.exists()
    assert output.with_name(output.name + ".gz").exists()
    assert output.with_name(output.stem + "_summary.json").exists()
