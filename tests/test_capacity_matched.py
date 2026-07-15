import json
from pathlib import Path
import subprocess
import sys

import pytest

from origins_dna.capacity_matched import (
    apply_projection,
    dimension_matched_comparison_summary,
    dimension_matched_representation_registry,
    make_signed_projection,
    paired_dimension_difference,
    run_dimension_matched_benchmark,
    write_dimension_matched_bundle,
)
from origins_dna.benchmark_types import SyntheticMotifConfig


def test_signed_projection_is_deterministic_and_dimensioned() -> None:
    first = make_signed_projection(4, 3, seed=7)
    second = make_signed_projection(4, 3, seed=7)
    assert first == second
    assert len(first) == 3
    assert {len(row) for row in first} == {4}
    assert len(apply_projection((1.0, 2.0, 3.0, 4.0), first)) == 3


def test_projection_validation() -> None:
    with pytest.raises(ValueError):
        make_signed_projection(0, 3, seed=1)
    with pytest.raises(TypeError):
        make_signed_projection(4, 3, seed=1.5)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="expected 2"):
        apply_projection((1.0,), ((1.0, 1.0),))
    with pytest.raises(ValueError, match="equal length"):
        apply_projection((1.0, 2.0), ((1.0, 1.0), (1.0,)))


def test_registry_has_exact_matched_dimensions() -> None:
    registry = dimension_matched_representation_registry(16, budgets=(8, 12), projection_seed=3)
    sequence = "ACGT" * 4
    assert len(registry["orbital_multilag_raw_48"](sequence)) == 48
    assert len(registry["tetranucleotide_raw_256"](sequence)) == 256
    for budget in (8, 12):
        assert len(registry[f"orbital_projected_{budget}"](sequence)) == budget
        assert len(registry[f"tetranucleotide_projected_{budget}"](sequence)) == budget


def test_registry_rejects_invalid_budgets() -> None:
    with pytest.raises(ValueError):
        dimension_matched_representation_registry(16, budgets=())
    with pytest.raises(ValueError):
        dimension_matched_representation_registry(16, budgets=(8, 8))
    with pytest.raises(ValueError):
        dimension_matched_representation_registry(3, budgets=(8,))


def _small_report():
    return run_dimension_matched_benchmark(
        SyntheticMotifConfig(
            sample_count=40,
            sequence_length=16,
            motif="CGCG",
            seeds=(3, 5),
            bootstrap_resamples=50,
        ),
        budgets=(8, 12),
        projection_seed=17,
    )


def test_dimension_matched_report_structure() -> None:
    report = _small_report()
    assert report.schema == "origins-dna.dimension-matched-benchmark.v1"
    assert report.representations == (
        "orbital_multilag_raw_48",
        "tetranucleotide_raw_256",
        "orbital_projected_8",
        "tetranucleotide_projected_8",
        "orbital_projected_12",
        "tetranucleotide_projected_12",
    )
    assert {result.condition for result in report.seed_results} == {
        "observed_labels",
        "label_permutation_null",
        "sequence_shuffle_null",
    }
    assert len(report.seed_results) == 2 * 3 * 6


def test_paired_summary_uses_identical_seed_pairs() -> None:
    report = _small_report()
    comparison = paired_dimension_difference(report, 8, seed=11)
    assert comparison["left"] == "tetranucleotide_projected_8"
    assert comparison["right"] == "orbital_projected_8"
    assert comparison["seed_count"] == 2
    summary = dimension_matched_comparison_summary(report)
    assert summary["budgets"] == [8, 12]
    assert len(summary["comparisons"]) == 2


def test_bundle_is_deterministic(tmp_path: Path) -> None:
    report = _small_report()
    first = write_dimension_matched_bundle(report, tmp_path / "first.json")
    second = write_dimension_matched_bundle(report, tmp_path / "second.json")
    for left, right in zip(first, second):
        assert left.read_bytes() == right.read_bytes()
    payload = json.loads(first[3].read_text(encoding="utf-8"))
    assert payload["schema"] == "origins-dna.dimension-matched-comparisons.v1"


def test_repository_script_runs_without_editable_install(tmp_path: Path) -> None:
    repository_root = Path(__file__).resolve().parents[1]
    output = tmp_path / "script.json"
    completed = subprocess.run(
        [
            sys.executable,
            str(repository_root / "scripts" / "run_dimension_matched_benchmark.py"),
            "--output",
            str(output),
            "--samples",
            "40",
            "--length",
            "16",
            "--seeds",
            "3,5",
            "--budgets",
            "8,12",
            "--bootstrap-resamples",
            "50",
        ],
        cwd=repository_root,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert output.exists()
    assert output.with_name(output.name + ".gz").exists()
    assert output.with_name(output.stem + "_summary.json").exists()
    assert output.with_name(output.stem + "_comparisons.json").exists()
