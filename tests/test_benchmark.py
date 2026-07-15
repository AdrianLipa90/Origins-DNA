import json
from pathlib import Path
import subprocess
import sys

import pytest

from origins_dna import (
    SyntheticMotifConfig,
    base_frequency_vector,
    bootstrap_mean_confidence_interval,
    confusion_matrix,
    dinucleotide_frequency_vector,
    fit_nearest_centroid,
    generate_synthetic_motif_dataset,
    orbital_classes_vector,
    orbital_full_vector,
    orbital_phase_vector,
    permute_labels,
    positional_one_hot_vector,
    representation_registry,
    run_controlled_benchmark,
    shuffle_sequence,
    stratified_split_indices,
)


def _small_config() -> SyntheticMotifConfig:
    return SyntheticMotifConfig(
        sample_count=40,
        sequence_length=24,
        motif="CGCG",
        seeds=(3, 5, 7),
        bootstrap_resamples=100,
    )


def test_config_rejects_odd_sample_count() -> None:
    with pytest.raises(ValueError, match="even"):
        SyntheticMotifConfig(sample_count=5)


def test_dataset_is_deterministic_balanced_and_composition_matched() -> None:
    config = _small_config()
    first = generate_synthetic_motif_dataset(config, seed=13)
    second = generate_synthetic_motif_dataset(config, seed=13)
    assert first == second
    assert first.labels.count(0) == first.labels.count(1) == 20
    assert sorted(first.motif) == sorted(first.decoy)
    assert all((first.motif in sequence) == bool(label) for sequence, label in zip(first.sequences, first.labels))


def test_shuffle_preserves_composition_and_is_deterministic() -> None:
    original = "AAAACCCCGGGGTTTT"
    shuffled = shuffle_sequence(original, seed=17)
    assert shuffled == shuffle_sequence(original, seed=17)
    assert sorted(shuffled) == sorted(original)


def test_label_permutation_preserves_class_counts() -> None:
    labels = (0, 0, 1, 1, 1)
    permuted = permute_labels(labels, seed=9)
    assert permuted == permute_labels(labels, seed=9)
    assert sorted(permuted) == sorted(labels)


def test_stratified_split_is_disjoint_and_class_complete() -> None:
    labels = (0, 0, 0, 0, 1, 1, 1, 1)
    train, test = stratified_split_indices(labels, train_fraction=0.5, seed=5)
    assert set(train).isdisjoint(test)
    assert set(train) | set(test) == set(range(len(labels)))
    assert {labels[index] for index in train} == {0, 1}
    assert {labels[index] for index in test} == {0, 1}


def test_representation_dimensions_are_stable() -> None:
    sequence = "ACGTACGT"
    assert len(base_frequency_vector(sequence)) == 4
    assert len(dinucleotide_frequency_vector(sequence)) == 16
    assert len(positional_one_hot_vector(sequence, 8)) == 32
    assert len(orbital_full_vector(sequence)) == 6
    assert len(orbital_phase_vector(sequence)) == 4
    assert len(orbital_classes_vector(sequence)) == 2
    assert tuple(representation_registry(8)) == (
        "base_frequency",
        "dinucleotide_frequency",
        "positional_one_hot",
        "orbital_full",
        "orbital_phase_only",
        "orbital_classes_only",
    )


def test_positional_one_hot_rejects_wrong_length() -> None:
    with pytest.raises(ValueError, match="expected_length"):
        positional_one_hot_vector("ACGT", 5)


def test_nearest_centroid_separates_simple_vectors() -> None:
    model = fit_nearest_centroid(((0.0,), (0.1,), (0.9,), (1.0,)), (0, 0, 1, 1))
    assert model.predict((0.05,)) == 0
    assert model.predict((0.95,)) == 1


def test_confusion_matrix_metrics() -> None:
    matrix = confusion_matrix((0, 0, 1, 1), (0, 1, 0, 1))
    assert matrix.accuracy == 0.5
    assert matrix.balanced_accuracy == 0.5
    assert matrix.matthews_correlation == 0.0


def test_bootstrap_interval_is_deterministic_and_contains_mean() -> None:
    values = (0.4, 0.5, 0.6, 0.7)
    first = bootstrap_mean_confidence_interval(values, confidence=0.95, resamples=200, seed=3)
    second = bootstrap_mean_confidence_interval(values, confidence=0.95, resamples=200, seed=3)
    assert first == second
    assert first[0] <= sum(values) / len(values) <= first[1]


def test_benchmark_is_reproducible_and_has_all_controls() -> None:
    config = _small_config()
    first = run_controlled_benchmark(config)
    second = run_controlled_benchmark(config)
    assert first == second
    conditions = {result.condition for result in first.seed_results}
    assert conditions == {
        "observed_labels",
        "label_permutation_null",
        "sequence_shuffle_null",
    }
    assert len(first.seed_results) == len(config.seeds) * len(first.representations) * 3
    assert all(0.0 <= result.balanced_accuracy <= 1.0 for result in first.seed_results)


def test_report_json_round_trip(tmp_path: Path) -> None:
    report = run_controlled_benchmark(_small_config())
    output = tmp_path / "report.json"
    full_json, archive, summary = report.write_bundle(output)
    payload = json.loads(full_json.read_text(encoding="utf-8"))
    summary_payload = json.loads(summary.read_text(encoding="utf-8"))
    assert payload["schema"] == "origins-dna.benchmark.v1"
    assert payload["benchmark_status"] == "synthetic_controlled_software_result"
    assert "biological mechanism" in payload["interpretation_boundary"]
    assert archive.read_bytes()[:2] == b"\x1f\x8b"
    assert "seed_results" not in summary_payload
    assert summary_payload["full_result_archive"] == archive.name


def test_bundle_archive_is_byte_reproducible(tmp_path: Path) -> None:
    report = run_controlled_benchmark(_small_config())
    _, first_archive, _ = report.write_bundle(tmp_path / "first.json")
    _, second_archive, _ = report.write_bundle(tmp_path / "second.json")
    assert first_archive.read_bytes() == second_archive.read_bytes()


def test_repository_script_runs_without_editable_install(tmp_path: Path) -> None:
    repository_root = Path(__file__).resolve().parents[1]
    output = tmp_path / "entrypoint.json"
    completed = subprocess.run(
        [
            sys.executable,
            str(repository_root / "scripts" / "run_controlled_benchmark.py"),
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
