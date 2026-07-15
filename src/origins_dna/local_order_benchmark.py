"""Local-order diagnostic benchmark introduced in Origins-DNA v0.4."""

from __future__ import annotations

import json
import math
from pathlib import Path

from .benchmark_data import derived_seed
from .benchmark_features import local_order_representation_registry
from .benchmark_runner import bootstrap_mean_confidence_interval, run_registered_benchmark
from .benchmark_types import BenchmarkReport, SyntheticMotifConfig

_COMPARISON_PAIRS = (
    ("orbital_multilag_tensor", "orbital_global_summary"),
    ("orbital_multilag_tensor", "orbital_pair_tensor"),
    ("orbital_multilag_tensor", "dinucleotide_frequency"),
    ("orbital_multilag_tensor", "trinucleotide_frequency"),
    ("tetranucleotide_frequency", "orbital_multilag_tensor"),
    ("orbital_multilag_tensor", "orbital_multilag_phase_only"),
    ("orbital_multilag_tensor", "orbital_multilag_classes_only"),
)


def run_local_order_benchmark(config: SyntheticMotifConfig) -> BenchmarkReport:
    """Compare local orbital tensors with explicit k-mer controls.

    This diagnostic uses the same deterministic synthetic task, held-out splits,
    null conditions, classifier and metrics as v0.3. It changes only the
    representation registry.
    """

    return run_registered_benchmark(
        config,
        local_order_representation_registry(config.sequence_length),
        schema="origins-dna.local-order-benchmark.v1",
        benchmark_status="synthetic_local_order_diagnostic_result",
        interpretation_boundary=(
            "Results diagnose information retention on one synthetic motif task only. "
            "They do not validate a biological mechanism and do not establish DNA decryption."
        ),
    )


def paired_representation_difference(
    report: BenchmarkReport,
    left: str,
    right: str,
    *,
    condition: str = "observed_labels",
    seed: int = 40_401,
) -> dict[str, object]:
    """Return a paired per-seed balanced-accuracy difference and bootstrap CI."""

    by_key = {
        (result.seed, result.representation, result.condition): result.balanced_accuracy
        for result in report.seed_results
    }
    differences: list[float] = []
    for benchmark_seed in report.config.seeds:
        try:
            left_score = by_key[(benchmark_seed, left, condition)]
            right_score = by_key[(benchmark_seed, right, condition)]
        except KeyError as error:
            raise ValueError(
                f"missing paired result for seed={benchmark_seed}, condition={condition}"
            ) from error
        differences.append(left_score - right_score)
    ci_low, ci_high = bootstrap_mean_confidence_interval(
        differences,
        confidence=report.config.confidence,
        resamples=report.config.bootstrap_resamples,
        seed=seed,
    )
    return {
        "left": left,
        "right": right,
        "condition": condition,
        "seed_count": len(differences),
        "mean_balanced_accuracy_difference": math.fsum(differences) / len(differences),
        "ci_low": ci_low,
        "ci_high": ci_high,
        "per_seed_differences": differences,
    }


def local_order_comparison_summary(report: BenchmarkReport) -> dict[str, object]:
    """Build the declared paired-comparison summary for the v0.4 report."""

    if report.schema != "origins-dna.local-order-benchmark.v1":
        raise ValueError("report must be a local-order benchmark")
    comparisons = [
        paired_representation_difference(
            report,
            left,
            right,
            seed=derived_seed(40_401, ordinal),
        )
        for ordinal, (left, right) in enumerate(_COMPARISON_PAIRS)
    ]
    return {
        "schema": "origins-dna.local-order-comparisons.v1",
        "source_report_schema": report.schema,
        "metric": "balanced_accuracy",
        "pairing": "identical held-out split and seed",
        "comparisons": comparisons,
        "interpretation_boundary": report.interpretation_boundary,
    }


def write_local_order_bundle(
    report: BenchmarkReport, path: str | Path
) -> tuple[Path, Path, Path, Path]:
    """Write the full report bundle plus paired-comparison summary."""

    full_json, archive, summary = report.write_bundle(path)
    summary_payload = json.loads(summary.read_text(encoding="utf-8"))
    summary.write_text(
        json.dumps(summary_payload, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    comparison_path = full_json.with_name(full_json.stem + "_comparisons.json")
    comparison_path.write_text(
        json.dumps(
            local_order_comparison_summary(report),
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n",
        encoding="utf-8",
    )
    return full_json, archive, summary, comparison_path
