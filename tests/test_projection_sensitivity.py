import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from origins_dna.benchmark_types import SyntheticMotifConfig
from origins_dna.projection_sensitivity import (
    EXPLORATORY_PILOT_PROJECTION_SEEDS,
    hierarchical_bootstrap_confidence_interval,
    projection_seed_schedule,
    run_projection_seed_sensitivity,
)


def _small_config() -> SyntheticMotifConfig:
    return SyntheticMotifConfig(
        sample_count=40,
        sequence_length=16,
        motif="CGCG",
        seeds=(3, 5),
        bootstrap_resamples=50,
    )


def test_projection_seed_schedule_is_deterministic_unique_and_independent() -> None:
    first = projection_seed_schedule(master_seed=70_701, count=32)
    second = projection_seed_schedule(master_seed=70_701, count=32)
    assert first == second
    assert len(first) == len(set(first)) == 32
    assert set(first).isdisjoint(EXPLORATORY_PILOT_PROJECTION_SEEDS)


def test_projection_seed_schedule_validation() -> None:
    with pytest.raises(TypeError, match="master_seed"):
        projection_seed_schedule(master_seed=1.5, count=2)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="positive"):
        projection_seed_schedule(master_seed=1, count=0)


def test_hierarchical_bootstrap_is_deterministic_and_contains_central_mean() -> None:
    values = ((-0.2, -0.1, 0.0), (0.1, 0.2, 0.3), (0.0, 0.1, 0.2))
    first = hierarchical_bootstrap_confidence_interval(
        values,
        confidence=0.95,
        resamples=500,
        seed=17,
    )
    second = hierarchical_bootstrap_confidence_interval(
        values,
        confidence=0.95,
        resamples=500,
        seed=17,
    )
    central = sum(sum(group) for group in values) / 9
    assert first == second
    assert first[0] <= central <= first[1]


def test_confirmatory_report_structure_and_controls() -> None:
    report = run_projection_seed_sensitivity(
        _small_config(),
        budgets=(8, 12),
        projection_seeds=(101, 202, 303),
        hierarchical_bootstrap_resamples=100,
    )
    assert report.schema == "origins-dna.projection-seed-sensitivity.v1"
    assert report.projection_seeds == (101, 202, 303)
    assert report.budgets == (8, 12)
    assert len(report.seed_results) == 3 * 3 * 2
    assert len(report.aggregate_results) == 3 * 2
    assert {result.condition for result in report.seed_results} == {
        "observed_labels",
        "label_permutation_null",
        "sequence_shuffle_null",
    }
    assert all(result.projection_seed_count == 3 for result in report.aggregate_results)
    assert all(result.data_seed_count == 2 for result in report.aggregate_results)


def test_confirmatory_report_rejects_exploratory_seed_overlap() -> None:
    with pytest.raises(ValueError, match="overlap exploratory"):
        run_projection_seed_sensitivity(
            _small_config(),
            budgets=(8,),
            projection_seeds=(50_501,),
            hierarchical_bootstrap_resamples=50,
        )


def test_bundle_is_byte_reproducible(tmp_path: Path) -> None:
    report = run_projection_seed_sensitivity(
        _small_config(),
        budgets=(8,),
        projection_seeds=(101, 202),
        hierarchical_bootstrap_resamples=50,
    )
    first_dir = tmp_path / "first"
    second_dir = tmp_path / "second"
    first = report.write_bundle(first_dir / "result.json")
    second = report.write_bundle(second_dir / "result.json")
    for left, right in zip(first, second):
        assert left.read_bytes() == right.read_bytes()
    summary = json.loads(first[2].read_text(encoding="utf-8"))
    assert summary["schema"] == "origins-dna.projection-seed-sensitivity.v1"
    assert "seed_results" not in summary


def test_repository_script_runs_without_editable_install(tmp_path: Path) -> None:
    repository_root = Path(__file__).resolve().parents[1]
    output = tmp_path / "sensitivity.json"
    completed = subprocess.run(
        [
            sys.executable,
            str(repository_root / "scripts" / "run_projection_seed_sensitivity.py"),
            "--output",
            str(output),
            "--samples",
            "40",
            "--length",
            "16",
            "--data-seeds",
            "3,5",
            "--budgets",
            "8",
            "--projection-seeds",
            "101,202",
            "--bootstrap-resamples",
            "50",
            "--hierarchical-bootstrap-resamples",
            "50",
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


def test_default_confirmatory_result_regression(tmp_path: Path) -> None:
    report = run_projection_seed_sensitivity(
        SyntheticMotifConfig(),
        hierarchical_bootstrap_resamples=5_000,
    )
    observed = {
        (result.condition, result.budget): result
        for result in report.aggregate_results
    }
    expected = {
        ("observed_labels", 12): (0.02333984375, 0.004589843750000008, 0.04248046874999998, 0.6875),
        ("observed_labels", 24): (0.04440104166666666, 0.030078124999999997, 0.058984375, 0.84375),
        ("observed_labels", 48): (0.07747395833333334, 0.064453125, 0.09075520833333334, 1.0),
        ("label_permutation_null", 12): (-0.014388020833333341, -0.030371093750000005, 0.0014648437499999959, 0.3125),
        ("label_permutation_null", 24): (-0.024641927083333338, -0.039355468750000004, -0.010058593750000002, 0.21875),
        ("label_permutation_null", 48): (-0.015657552083333335, -0.02910156250000001, -0.0018880208333333338, 0.25),
        ("sequence_shuffle_null", 12): (0.01181640625, 0.00208333333333333, 0.021354166666666657, 0.71875),
        ("sequence_shuffle_null", 24): (0.009863281250000001, 0.0005208333333333367, 0.01936848958333334, 0.65625),
        ("sequence_shuffle_null", 48): (0.009570312499999997, 0.0009765625000000028, 0.018098958333333335, 0.71875),
    }
    assert set(observed) == set(expected)
    for key, (mean, low, high, positive_fraction) in expected.items():
        result = observed[key]
        assert result.mean_balanced_accuracy_difference == pytest.approx(mean, abs=1e-15)
        assert result.hierarchical_ci_low == pytest.approx(low, abs=1e-15)
        assert result.hierarchical_ci_high == pytest.approx(high, abs=1e-15)
        assert result.fraction_projection_means_positive == positive_fraction

    full_json, archive, summary = report.write_bundle(
        tmp_path / "v0.6_projection_seed_sensitivity.json"
    )
    assert hashlib.sha256(full_json.read_bytes()).hexdigest() == (
        "feeb5350584d9328474faa7886017b197268f6ccb0655d908ffa1f88f2f9802c"
    )
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == (
        "2e4710b9b2d05fe81157da40ad1a38dfe0a84ecca2d87543d50903defe5de4d4"
    )
    assert hashlib.sha256(summary.read_bytes()).hexdigest() == (
        "dc27c033028b17ed50f517715b1e74b226dd5914248cbbc57818d583c329599c"
    )
