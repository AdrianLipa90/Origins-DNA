import json
from pathlib import Path
import subprocess
import sys

import pytest

from origins_dna import (
    SyntheticMotifConfig,
    kmer_frequency_vector,
    local_order_representation_registry,
    orbital_lag_tensor_vector,
    orbital_multilag_classes_vector,
    orbital_multilag_phase_vector,
    orbital_multilag_tensor_vector,
    orbital_pair_tensor_vector,
    run_controlled_benchmark,
    run_local_order_benchmark,
)


def _small_config() -> SyntheticMotifConfig:
    return SyntheticMotifConfig(
        sample_count=40,
        sequence_length=24,
        motif="CGCG",
        seeds=(3, 5),
        bootstrap_resamples=50,
    )


def test_kmer_frequency_vector_dimensions_and_mass() -> None:
    sequence = "ACGTACGT"
    for k, dimension in ((1, 4), (2, 16), (3, 64), (4, 256)):
        vector = kmer_frequency_vector(sequence, k)
        assert len(vector) == dimension
        assert sum(vector) == pytest.approx(1.0)


def test_kmer_frequency_vector_validation() -> None:
    with pytest.raises(TypeError, match="integer"):
        kmer_frequency_vector("ACGT", 2.0)
    with pytest.raises(ValueError, match="positive"):
        kmer_frequency_vector("ACGT", 0)


def test_orbital_tensor_dimensions_are_declared() -> None:
    sequence = "ACGTACGT"
    assert len(orbital_pair_tensor_vector(sequence)) == 16
    assert len(orbital_multilag_tensor_vector(sequence)) == 48
    assert len(orbital_multilag_phase_vector(sequence)) == 12
    assert len(orbital_multilag_classes_vector(sequence)) == 12


def test_orbital_lag_tensor_is_directional() -> None:
    assert orbital_pair_tensor_vector("CGCG") != orbital_pair_tensor_vector("GCGC")


def test_orbital_lag_tensor_validation() -> None:
    with pytest.raises(ValueError, match="must not be empty"):
        orbital_lag_tensor_vector("ACGT", lags=())
    with pytest.raises(ValueError, match="positive"):
        orbital_lag_tensor_vector("ACGT", lags=(0,))
    with pytest.raises(ValueError, match="unique"):
        orbital_lag_tensor_vector("ACGT", lags=(1, 1))
    with pytest.raises(ValueError, match="smaller"):
        orbital_lag_tensor_vector("ACGT", lags=(4,))
    with pytest.raises(ValueError, match="between 0 and 3"):
        orbital_lag_tensor_vector("ACGT", lags=(1,), channels=(4,))


def test_local_order_registry_is_stable() -> None:
    registry = local_order_representation_registry(24)
    assert tuple(registry) == (
        "dinucleotide_frequency",
        "trinucleotide_frequency",
        "tetranucleotide_frequency",
        "orbital_global_summary",
        "orbital_pair_tensor",
        "orbital_multilag_tensor",
        "orbital_multilag_phase_only",
        "orbital_multilag_classes_only",
    )
    assert {name: len(function("ACGT" * 6)) for name, function in registry.items()} == {
        "dinucleotide_frequency": 16,
        "trinucleotide_frequency": 64,
        "tetranucleotide_frequency": 256,
        "orbital_global_summary": 6,
        "orbital_pair_tensor": 16,
        "orbital_multilag_tensor": 48,
        "orbital_multilag_phase_only": 12,
        "orbital_multilag_classes_only": 12,
    }


def test_local_order_registry_rejects_too_short_sequences() -> None:
    with pytest.raises(ValueError, match="maximum lag"):
        local_order_representation_registry(3)


def test_local_order_benchmark_is_reproducible_and_controlled() -> None:
    first = run_local_order_benchmark(_small_config())
    second = run_local_order_benchmark(_small_config())
    assert first == second
    assert first.schema == "origins-dna.local-order-benchmark.v1"
    assert first.benchmark_status == "synthetic_local_order_diagnostic_result"
    assert {result.condition for result in first.seed_results} == {
        "observed_labels",
        "label_permutation_null",
        "sequence_shuffle_null",
    }
    assert len(first.seed_results) == 2 * 8 * 3


def test_v03_registry_and_schema_remain_frozen() -> None:
    report = run_controlled_benchmark(_small_config())
    assert report.schema == "origins-dna.benchmark.v1"
    assert report.representations == (
        "base_frequency",
        "dinucleotide_frequency",
        "positional_one_hot",
        "orbital_full",
        "orbital_phase_only",
        "orbital_classes_only",
    )


def test_local_order_report_bundle(tmp_path: Path) -> None:
    report = run_local_order_benchmark(_small_config())
    output = tmp_path / "local_order.json"
    full_json, archive, summary = report.write_bundle(output)
    assert json.loads(full_json.read_text(encoding="utf-8"))["schema"] == (
        "origins-dna.local-order-benchmark.v1"
    )
    assert archive.read_bytes()[:2] == b"\x1f\x8b"
    assert "seed_results" not in json.loads(summary.read_text(encoding="utf-8"))


def test_local_order_script_runs_without_editable_install(tmp_path: Path) -> None:
    repository_root = Path(__file__).resolve().parents[1]
    output = tmp_path / "local_order.json"
    completed = subprocess.run(
        [
            sys.executable,
            str(repository_root / "scripts" / "run_local_order_benchmark.py"),
            "--output",
            str(output),
            "--samples",
            "40",
            "--length",
            "24",
            "--seeds",
            "3,5",
            "--bootstrap-resamples",
            "50",
        ],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert output.exists()
    assert output.with_suffix(output.suffix + ".gz").exists()
    assert output.with_name(f"{output.stem}_summary.json").exists()
    assert output.with_name(f"{output.stem}_comparisons.json").exists()


def test_paired_comparison_summary_is_deterministic() -> None:
    from origins_dna import local_order_comparison_summary

    report = run_local_order_benchmark(_small_config())
    first = local_order_comparison_summary(report)
    second = local_order_comparison_summary(report)
    assert first == second
    assert first["schema"] == "origins-dna.local-order-comparisons.v1"
    assert len(first["comparisons"]) == 7
    assert all(item["seed_count"] == 2 for item in first["comparisons"])


def test_paired_difference_rejects_missing_representation() -> None:
    from origins_dna import paired_representation_difference

    report = run_local_order_benchmark(_small_config())
    with pytest.raises(ValueError, match="missing paired result"):
        paired_representation_difference(report, "missing", "dinucleotide_frequency")
