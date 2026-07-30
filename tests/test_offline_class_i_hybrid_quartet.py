import unittest
from argparse import Namespace
from pathlib import Path
from tempfile import TemporaryDirectory

from diagnostics.offline_class_i_hybrid_quartet import (
    compare_quartets,
    pair_evidence,
    proposal_is_supported,
    read_tsv,
    write_tsv,
)


class ClassIHybridQuartetTest(unittest.TestCase):
    def test_paired_end_concordance_increases_same_allele_evidence(self):
        owners = {"A" * 31: ("A*01:01",), "C" * 31: ("A*02:01",)}
        evidence, private = pair_evidence(
            "A" * 31, "I" * 31, "A" * 31, "I" * 31, owners, 31, 0.5
        )
        self.assertGreater(evidence["A*01:01"], 2.0)
        self.assertEqual(private, {"A*01:01"})

    def test_read_evidence_prefers_supported_repeated_copy(self):
        baseline = ("A*01:01", "A*02:01", "A*03:01", "A*24:02")
        proposal = ("A*01:01", "A*01:01", "A*02:01", "A*03:01")
        rows = [
            {"A*01:01": 8.0},
            {"A*01:01": 7.0},
            {"A*02:01": 6.0},
            {"A*03:01": 5.0},
        ]
        result = compare_quartets(rows, baseline, proposal)
        self.assertGreater(result["log_bayes_factor"], 0.0)

    def test_quartet_comparison_is_source_label_free(self):
        quartet = ("A*01:01", "A*01:01", "A*02:01", "A*03:01")
        rows = [{"A*01:01": 5.0}, {"A*02:01": 4.0}, {"A*03:01": 3.0}]
        left = compare_quartets(rows, quartet, tuple(reversed(quartet)))
        self.assertAlmostEqual(left["log_bayes_factor"], 0.0)

    def test_write_tsv_accepts_branch_specific_fields(self):
        rows = [{"decision": "same"}, {"decision": "proposal", "total_pairs": 12}]
        with TemporaryDirectory() as directory:
            path = Path(directory) / "hybrid.tsv"
            write_tsv(path, rows)
            written = read_tsv(path)
        self.assertEqual(len(written), 2)
        self.assertEqual(written[0]["total_pairs"], "")
        self.assertEqual(written[1]["total_pairs"], "12")

    def test_depth_normalized_gate_is_invariant_to_pair_count_scaling(self):
        args = Namespace(
            min_log_bayes_factor=-float("inf"),
            min_log_bayes_factor_per_informative_pair=0.1,
            min_discriminating_pairs=3,
            min_proposal_private_pairs=10,
            min_private_pair_ratio=5.0,
            private_pair_slack=0,
        )
        shallow = {
            "log_bayes_factor": 120.0,
            "informative_pairs": 1000,
            "discriminating_pairs": 50,
        }
        deep = {
            "log_bayes_factor": 360.0,
            "informative_pairs": 3000,
            "discriminating_pairs": 150,
        }
        self.assertTrue(proposal_is_supported(shallow, 2, 20, args))
        self.assertTrue(proposal_is_supported(deep, 6, 60, args))

    def test_depth_normalized_gate_rejects_weak_private_ratio(self):
        args = Namespace(
            min_log_bayes_factor=-float("inf"),
            min_log_bayes_factor_per_informative_pair=0.1,
            min_discriminating_pairs=3,
            min_proposal_private_pairs=10,
            min_private_pair_ratio=5.0,
            private_pair_slack=0,
        )
        comparison = {
            "log_bayes_factor": 200.0,
            "informative_pairs": 1000,
            "discriminating_pairs": 50,
        }
        self.assertFalse(proposal_is_supported(comparison, 11, 46, args))


if __name__ == "__main__":
    unittest.main()