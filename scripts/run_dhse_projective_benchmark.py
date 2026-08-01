#!/usr/bin/env python3
"""Run Origins-DNA v0.8 DHSE projective-balance validation and benchmark."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
_SOURCE_ROOT = _REPOSITORY_ROOT / "src"
if str(_SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(_SOURCE_ROOT))

from origins_dna.dhse_benchmark import run_dhse_projective_benchmark


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/v0.8_dhse_projective_benchmark.json"))
    parser.add_argument("--samples", type=int, default=400)
    parser.add_argument("--length", type=int, default=64)
    parser.add_argument("--motif", default="CGCG")
    parser.add_argument("--projection-seed", type=int, default=100_801)
    parser.add_argument("--bootstrap-resamples", type=int, default=2_000)
    args = parser.parse_args()
    paths = run_dhse_projective_benchmark(
        output=args.output,
        sample_count=args.samples,
        sequence_length=args.length,
        motif=args.motif,
        projection_seed=args.projection_seed,
        bootstrap_resamples=args.bootstrap_resamples,
    )
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
