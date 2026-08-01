from fractions import Fraction
import math

from origins_dna.dhse_exact import EXPECTED, exact_forcing_response
from origins_dna.dhse_projective import (
    DUAL_CHANNELS,
    INVARIANT_CHANNEL,
    projective_odds,
    projective_window_odds,
    reverse_complement_window_covariance,
    self_dual_residual,
    signed_log_odds,
    verify_complement_covariance,
)
from origins_dna.orbital_atomic import complement


def test_projective_odds_are_exact_positive_rationals() -> None:
    assert projective_odds("AAAT", "AT") == Fraction(4, 2)
    assert projective_odds("CCCG", "CG") == Fraction(4, 2)
    assert projective_odds("AGCT", "RY") == Fraction(3, 3)
    assert projective_odds("ACGT", INVARIANT_CHANNEL) == Fraction(3, 3)


def test_complement_inverts_every_dual_channel_exactly() -> None:
    sequence = "AAACCGTACGTT"
    transformed = complement(sequence)
    for channel in DUAL_CHANNELS:
        q = projective_odds(sequence, channel)
        assert projective_odds(transformed, channel) == 1 / q
        assert self_dual_residual(1 / q) == self_dual_residual(q)
        assert math.isclose(signed_log_odds(1 / q), -signed_log_odds(q), abs_tol=1e-12)


def test_ws_channel_is_complement_invariant() -> None:
    sequence = "AAAACCGT"
    assert projective_odds(complement(sequence), INVARIANT_CHANNEL) == projective_odds(
        sequence, INVARIANT_CHANNEL
    )


def test_exact_dhse_response_matches_stage_m_at_self_dual_point() -> None:
    response = exact_forcing_response(Fraction(1))
    assert response == tuple(EXPECTED[length]["q1_count"] for length in (1, 2, 3, 4))


def test_exact_dhse_response_is_reciprocal() -> None:
    for q in (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(4, 3), Fraction(2), Fraction(4)):
        assert exact_forcing_response(q) == exact_forcing_response(1 / q)


def test_window_rules_and_reverse_complement_covariance() -> None:
    sequence = "AACCGGTTAACCGGTT"
    assert projective_window_odds(sequence, 99) == ()
    assert reverse_complement_window_covariance(sequence, 3)


def test_full_covariance_gate() -> None:
    assert all(verify_complement_covariance("AAACCGTACGTT").values())
