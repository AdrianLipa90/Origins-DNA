#!/usr/bin/env python3
"""Run the preregistered Origins-DNA v0.9 causal DHSE benchmark."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from origins_dna.causal_dhse_benchmark import run_causal_dhse_benchmark


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/v0.9_causal_dhse_benchmark.json"))
    parser.add_argument("--samples", type=int, default=400)
    parser.add_argument("--length", type=int, default=64)
    parser.add_argument("--motif", default="CGCG")
    parser.add_argument("--projection-seed", type=int, default=109_001)
    parser.add_argument("--bootstrap-resamples", type=int, default=2_000)
    args = parser.parse_args()
    paths = run_causal_dhse_benchmark(
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
