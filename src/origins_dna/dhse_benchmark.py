"""Origins-DNA v0.8 DHSE projective-balance benchmark."""

from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path
from typing import Sequence

from .benchmark_data import derived_seed
from .benchmark_features import kmer_frequency_vector, orbital_multilag_tensor_vector
from .benchmark_runner import run_registered_benchmark
from .benchmark_types import SyntheticMotifConfig
from .capacity_matched import DEFAULT_DIMENSION_BUDGETS, projected_representation
from .dhse_exact import exact_stage_m_receipt
from .dhse_projective import (
    dhse_projective_vector,
    dhse_without_signed_log_vector,
    invariant_ws_control_vector,
    projective_only_vector,
    verify_complement_covariance,
)


def _vector_dimension(representation, sequence_length: int) -> int:
    probe = "ACGT" * ((sequence_length + 3) // 4)
    return len(representation(probe[:sequence_length]))


def dhse_representation_registry(
    sequence_length: int,
    *,
    budgets: Sequence[int] = DEFAULT_DIMENSION_BUDGETS,
    projection_seed: int = 100_801,
):
    if sequence_length < 24:
        raise ValueError("sequence_length must be at least 24 for the frozen v0.8 windows")

    @lru_cache(maxsize=None)
    def orbital(sequence: str):
        return orbital_multilag_tensor_vector(sequence)

    @lru_cache(maxsize=None)
    def tetranucleotide(sequence: str):
        return kmer_frequency_vector(sequence, 4)

    @lru_cache(maxsize=None)
    def dhse(sequence: str):
        return dhse_projective_vector(sequence)

    @lru_cache(maxsize=None)
    def projective(sequence: str):
        return projective_only_vector(sequence)

    @lru_cache(maxsize=None)
    def dhse_no_sign(sequence: str):
        return dhse_without_signed_log_vector(sequence)

    @lru_cache(maxsize=None)
    def ws_control(sequence: str):
        return invariant_ws_control_vector(sequence)

    families = {
        "orbital": orbital,
        "tetranucleotide": tetranucleotide,
        "dhse": dhse,
        "projective_only": projective,
        "dhse_without_signed_log": dhse_no_sign,
        "ws_invariant_control": ws_control,
    }
    dimensions = {
        name: _vector_dimension(representation, sequence_length)
        for name, representation in families.items()
    }
    registry = {}
    for ordinal, budget in enumerate(sorted(tuple(budgets))):
        family_seed = derived_seed(projection_seed, ordinal)
        for family_ordinal, (name, representation) in enumerate(families.items(), start=1):
            registry[f"{name}_projected_{budget}"] = projected_representation(
                representation,
                input_dimension=dimensions[name],
                output_dimension=budget,
                seed=derived_seed(family_seed, family_ordinal),
            )
    return registry


def run_dhse_projective_benchmark(
    *,
    output: str | Path,
    sample_count: int = 400,
    sequence_length: int = 64,
    motif: str = "CGCG",
    projection_seed: int = 100_801,
    bootstrap_resamples: int = 2_000,
    budgets: Sequence[int] = DEFAULT_DIMENSION_BUDGETS,
) -> tuple[Path, Path, Path, Path, Path]:
    """Execute exact Stage-M validation, covariance gates and the frozen task."""

    receipt = exact_stage_m_receipt()
    covariance_sequences = (
        "AAACCGTACGTT",
        "ACGTACGTACGT",
        "GGGAAATTTCCC",
        "CGCGATATCGCG",
    )
    covariance = {
        sequence: verify_complement_covariance(sequence)
        for sequence in covariance_sequences
    }
    covariance_passed = all(all(row.values()) for row in covariance.values())
    if receipt["technical_status"] != "PASS" or not covariance_passed:
        raise AssertionError("DHSE exact or complement-covariance gate failed")

    config = SyntheticMotifConfig(
        sample_count=sample_count,
        sequence_length=sequence_length,
        motif=motif,
        bootstrap_resamples=bootstrap_resamples,
    )
    report = run_registered_benchmark(
        config,
        dhse_representation_registry(
            sequence_length,
            budgets=budgets,
            projection_seed=projection_seed,
        ),
        schema="origins-dna.dhse-projective-benchmark.v1",
        benchmark_status="synthetic_exact_dhse_projective_balance_result",
        interpretation_boundary=(
            "This benchmark tests a frozen exact DHSE-derived representation on a synthetic motif task. "
            "Passing software and covariance gates does not establish a biological mechanism, abiogenesis "
            "pathway or DNA decryption. Predictive comparisons and null conditions must be interpreted separately."
        ),
    )
    output_path = Path(output)
    full_json, archive, summary = report.write_bundle(output_path)
    receipt_path = output_path.with_name(output_path.stem + "_dhse_receipt.json")
    receipt_path.write_text(
        json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    covariance_path = output_path.with_name(output_path.stem + "_covariance.json")
    covariance_path.write_text(
        json.dumps(
            {
                "schema": "origins-dna.dhse-covariance.v1",
                "sequences": covariance,
                "all_gates_passed": covariance_passed,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n",
        encoding="utf-8",
    )
    return full_json, archive, summary, receipt_path, covariance_path
