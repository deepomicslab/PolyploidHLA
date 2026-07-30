import tempfile
import unittest
from collections import deque
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import apply_quartet_optimization as optimization


class QuartetOptimizationTest(unittest.TestCase):
    def test_slot_quartet_maps_major_group_to_recipient(self):
        baseline = ("A*01:01", "A*02:01", "A*03:01", "A*04:01")

        result = optimization.slot_quartet(
            ("A*02:01", "A*05:01"),
            ("A*03:01", "A*06:01"),
            0.8,
            baseline,
        )

        self.assertEqual(result, ("A*05:01", "A*02:01", "A*03:01", "A*06:01"))

    def test_slot_quartet_maps_major_group_to_donor_when_recipient_is_minor(self):
        baseline = ("A*01:01", "A*02:01", "A*03:01", "A*04:01")

        result = optimization.slot_quartet(
            ("A*03:01", "A*05:01"),
            ("A*01:01", "A*06:01"),
            0.2,
            baseline,
        )

        self.assertEqual(result, ("A*01:01", "A*06:01", "A*03:01", "A*05:01"))

    def test_lift_quartet_requires_mapping_for_every_allele(self):
        mapping = {"A*01:01": deque(["A*01:01:01"])}

        self.assertIsNone(optimization.lift_quartet(("A*01:01", "A*02:01"), mapping))

    def test_class_i_gate_falls_back_when_gene_fastqs_are_missing(self):
        with tempfile.TemporaryDirectory() as directory:
            args = SimpleNamespace(
                spechla_root=Path(directory),
                sample="sample",
                imgt=Path(directory) / "missing.fa",
            )
            accepted, audit = optimization.class_i_read_gate(
                args,
                "HLA-A",
                ("A*01:01", "A*02:01", "A*03:01", "A*04:01"),
                ("A*01:01", "A*02:01", "A*03:01", "A*05:01"),
            )

        self.assertFalse(accepted)
        self.assertEqual(audit["decision"], "fallback")

    def test_unexpected_gene_error_is_isolated_as_fallback(self):
        args = SimpleNamespace(sample="sample", profile="shadow")

        with patch.object(optimization, "optimize_gene", side_effect=RuntimeError("test")):
            audit = optimization.safely_optimize_gene(args, "HLA-A", 0.8)

        self.assertEqual(audit["decision"], "fallback")
        self.assertEqual(audit["reason"], "unexpected_error:RuntimeError")
        self.assertEqual(audit["applied"], "0")


if __name__ == "__main__":
    unittest.main()