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
