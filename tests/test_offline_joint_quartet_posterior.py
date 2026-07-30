import unittest

from diagnostics.offline_joint_quartet_posterior import (
    call_quartet,
    expected_fractions,
    unique_mixture_groupings,
)


class JointQuartetPosteriorTest(unittest.TestCase):
    def test_shared_allele_partitions_are_preserved(self):
        quartet = ("A*01:01", "A*01:01", "A*02:01", "A*03:01")
        partitions = unique_mixture_groupings(quartet)
        self.assertIn(
            (("A*01:01", "A*02:01"), ("A*01:01", "A*03:01")),
            partitions,
        )
        fractions = expected_fractions(partitions[0], 0.70)
        self.assertAlmostEqual(sum(fractions.values()), 1.0)

    def test_caller_recovers_repeated_shared_allele(self):
        counts = {
            "A*01:01": 500.0,
            "A*02:01": 350.0,
            "A*03:01": 150.0,
            "A*24:02": 5.0,
        }
        result = call_quartet(
            counts,
            ("A*01:01", "A*01:01", "A*02:01", "A*24:02"),
            major_fraction_prior=0.70,
            top_n=4,
        )
        self.assertEqual(
            result["quartet"],
            ("A*01:01", "A*01:01", "A*02:01", "A*03:01"),
        )

    def test_multiset_call_is_source_label_invariant(self):
        counts = {"A*01:01": 500.0, "A*02:01": 350.0, "A*03:01": 150.0}
        left = call_quartet(counts, tuple(counts), major_fraction_prior=0.70, top_n=3)
        right = call_quartet(
            counts, tuple(reversed(counts)), major_fraction_prior=0.70, top_n=3
        )
        self.assertEqual(left["quartet"], right["quartet"])

    def test_expands_candidates_for_supported_low_frequency_copy(self):
        counts = {
            "A*01:01": 500.0,
            "A*02:01": 300.0,
            "A*03:01": 190.0,
            "A*24:02": 4.0,
            "A*26:01": 3.0,
        }
        result = call_quartet(
            counts,
            ("A*01:01", "A*01:01", "A*02:01", "A*03:01"),
            major_fraction_prior=0.90,
            top_n=3,
            max_top_n=5,
            candidate_min_fraction=0.002,
        )
        self.assertGreaterEqual(result["candidate_count"], 5)


if __name__ == "__main__":
    unittest.main()