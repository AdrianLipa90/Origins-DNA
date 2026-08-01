from __future__ import annotations

from fractions import Fraction

from origins_dna.causal_dhse import (
    BASES,
    base_operator_assignment,
    causal_dhse_vector,
    causal_projective_only_vector,
    causal_trajectory,
    conjugate_product_identity,
    determinant,
    double_reverse_complement_identity,
    operator_word,
    order_is_noncommutative,
    prefix_products,
    projective_image,
    reciprocal_conjugate,
    reciprocal_conjugation_identity,
    reverse_complement_expected_terminal,
    reverse_complement_terminal_identity,
    terminal_only_vector,
    terminal_product,
)
from origins_dna.orbital_atomic import complement, reverse_complement


def test_assignment_is_deterministic_distinct_and_complement_paired() -> None:
    assignment = base_operator_assignment()
    assert tuple(assignment) == tuple(BASES)
    assert len(set(assignment.values())) == 4
    assert assignment["T"] == reciprocal_conjugate(assignment["A"])
    assert assignment["G"] == reciprocal_conjugate(assignment["C"])
    assert all(determinant(matrix) > 0 for matrix in assignment.values())
    assert assignment == base_operator_assignment()


def test_prefix_products_are_ordered_and_exact() -> None:
    sequence = "ACGT"
    word = operator_word(sequence)
    products = prefix_products(sequence)
    assert len(products) == len(sequence)
    assert products[-1] == terminal_product(sequence)
    assert all(isinstance(value, int) for matrix in products for value in matrix)
    assert projective_image(products[0], Fraction(1)) > 0
    assert word != tuple(reversed(word))


def test_count_matched_words_are_noncommutative() -> None:
    assert sorted("CGCG") == sorted("CCGG")
    assert order_is_noncommutative("CGCG", "CCGG")
    assert terminal_product("CGCG") != terminal_product("CCGG")
    assert causal_trajectory("CGCG") != causal_trajectory("CCGG")


def test_reciprocal_conjugation_is_an_involutive_automorphism() -> None:
    matrices = tuple(base_operator_assignment().values())
    for matrix in matrices:
        assert reciprocal_conjugate(reciprocal_conjugate(matrix)) == matrix
        for right in matrices:
            assert reciprocal_conjugation_identity(matrix, right)


def test_complement_is_operatorwise_and_terminal_conjugation() -> None:
    for sequence in ("A", "CGCG", "AACGTCGA", "ATGCCGTA"):
        assert conjugate_product_identity(sequence)
        expected_word = tuple(
            reciprocal_conjugate(matrix) for matrix in operator_word(sequence)
        )
        assert operator_word(complement(sequence)) == expected_word
        assert terminal_product(complement(sequence)) == reciprocal_conjugate(
            terminal_product(sequence)
        )


def test_reverse_complement_terminal_identity() -> None:
    for sequence in ("A", "CGCG", "AACGTCGA", "ATGCCGTA"):
        assert reverse_complement_terminal_identity(sequence)
        assert terminal_product(
            reverse_complement(sequence)
        ) == reverse_complement_expected_terminal(sequence)


def test_double_reverse_complement_restores_word_and_trajectory() -> None:
    sequence = "ACGTCGAT"
    assert double_reverse_complement_identity(sequence)
    restored = reverse_complement(reverse_complement(sequence))
    assert operator_word(restored) == operator_word(sequence)
    assert causal_trajectory(restored) == causal_trajectory(sequence)


def test_feature_vectors_are_finite_fixed_and_order_sensitive() -> None:
    left = "CGCG" * 8
    right = "CCGG" * 8
    for representation in (
        causal_dhse_vector,
        causal_projective_only_vector,
        terminal_only_vector,
    ):
        left_vector = representation(left)
        right_vector = representation(right)
        assert left_vector
        assert len(left_vector) == len(right_vector)
        assert left_vector != right_vector
        assert all(
            value == value and value not in (float("inf"), float("-inf"))
            for value in left_vector
        )
