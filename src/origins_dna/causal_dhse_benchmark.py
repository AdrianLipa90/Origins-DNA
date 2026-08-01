"""Frozen Origins-DNA v0.9 causal DHSE benchmark."""
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
from .causal_dhse import (
    causal_dhse_vector,
    causal_projective_only_vector,
    conjugate_product_identity,
    double_reverse_complement_identity,
    order_is_noncommutative,
    reverse_complement_terminal_identity,
    terminal_only_vector,
)
from .dhse_projective import dhse_projective_vector
from .encoding import normalize_dna


def _dimension(representation, sequence_length: int) -> int:
    probe = ("ACGT" * ((sequence_length + 3) // 4))[:sequence_length]
    return len(representation(probe))


def count_matched_operator_vector(sequence: str) -> tuple[float, ...]:
    """Destroy order deterministically while preserving exact base counts."""
    return causal_dhse_vector("".join(sorted(normalize_dna(sequence))))


def reversed_causal_vector(sequence: str) -> tuple[float, ...]:
    return causal_dhse_vector(normalize_dna(sequence)[::-1])


def causal_representation_registry(
    sequence_length: int,
    *,
    budgets: Sequence[int] = DEFAULT_DIMENSION_BUDGETS,
    projection_seed: int = 109_001,
):
    if sequence_length < 24:
        raise ValueError("sequence_length must be at least 24")

    @lru_cache(maxsize=None)
    def orbital(sequence: str):
        return orbital_multilag_tensor_vector(sequence)

    @lru_cache(maxsize=None)
    def tetranucleotide(sequence: str):
        return kmer_frequency_vector(sequence, 4)

    @lru_cache(maxsize=None)
    def v08(sequence: str):
        return dhse_projective_vector(sequence)

    @lru_cache(maxsize=None)
    def causal(sequence: str):
        return causal_dhse_vector(sequence)

    @lru_cache(maxsize=None)
    def causal_projective(sequence: str):
        return causal_projective_only_vector(sequence)

    @lru_cache(maxsize=None)
    def terminal(sequence: str):
        return terminal_only_vector(sequence)

    @lru_cache(maxsize=None)
    def count_matched(sequence: str):
        return count_matched_operator_vector(sequence)

    @lru_cache(maxsize=None)
    def reversed_order(sequence: str):
        return reversed_causal_vector(sequence)

    families = {
        "tetranucleotide": tetranucleotide,
        "orbital": orbital,
        "composition_dhse_v08": v08,
        "causal_dhse": causal,
        "causal_projective_only": causal_projective,
        "causal_terminal_only": terminal,
        "causal_count_matched": count_matched,
        "causal_reversed_order": reversed_order,
    }
    dimensions = {
        name: _dimension(representation, sequence_length)
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


def exact_causal_gate_receipt() -> dict[str, object]:
    sequences = ("CGCG", "CCGG", "AACGTCGA", "ATGCCGTA")
    gates = {
        sequence: {
            "complement_operatorwise": conjugate_product_identity(sequence),
            "reverse_complement_terminal": reverse_complement_terminal_identity(sequence),
            "double_reverse_complement": double_reverse_complement_identity(sequence),
        }
        for sequence in sequences
    }
    noncommutative = order_is_noncommutative("CGCG", "CCGG")
    passed = noncommutative and all(all(row.values()) for row in gates.values())
    return {
        "schema": "origins-dna.causal-dhse-gates.v1",
        "technical_status": "PASS" if passed else "FAIL",
        "count_matched_noncommutativity_CGCG_vs_CCGG": noncommutative,
        "sequence_gates": gates,
    }


def run_causal_dhse_benchmark(
    *,
    output: str | Path,
    sample_count: int = 400,
    sequence_length: int = 64,
    motif: str = "CGCG",
    projection_seed: int = 109_001,
    bootstrap_resamples: int = 2_000,
    budgets: Sequence[int] = DEFAULT_DIMENSION_BUDGETS,
) -> tuple[Path, Path, Path, Path]:
    receipt = exact_causal_gate_receipt()
    if receipt["technical_status"] != "PASS":
        raise AssertionError("causal DHSE exact gates failed")
    config = SyntheticMotifConfig(
        sample_count=sample_count,
        sequence_length=sequence_length,
        motif=motif,
        bootstrap_resamples=bootstrap_resamples,
    )
    report = run_registered_benchmark(
        config,
        causal_representation_registry(
            sequence_length,
            budgets=budgets,
            projection_seed=projection_seed,
        ),
        schema="origins-dna.causal-dhse-benchmark.v1",
        benchmark_status="synthetic_preregistered_causal_dhse_result",
        interpretation_boundary=(
            "This benchmark tests whether ordered noncommutative DHSE trajectories carry synthetic "
            "sequence information beyond composition. A positive causal contrast is not proof that "
            "biological DNA literally executes the declared Möbius operators."
        ),
    )
    output_path = Path(output)
    full_json, archive, summary = report.write_bundle(output_path)
    receipt_path = output_path.with_name(output_path.stem + "_causal_gates.json")
    receipt_path.write_text(
        json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    return full_json, archive, summary, receipt_path
