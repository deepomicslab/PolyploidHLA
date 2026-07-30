import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import apply_class2_joint_rescue as rescue


class Class2RescueGeneScopeTest(unittest.TestCase):
    def test_drb1_only_scope_skips_dqb1_and_dpb1(self):
        rows = [
            {"gene": "HLA-DRB1"},
            {"gene": "HLA-DQB1"},
            {"gene": "HLA-DPB1"},
        ]
        args = SimpleNamespace(rescue_genes={"HLA-DRB1"})
        drb1_proposal = {
            "current_2field": ["DRB1*01:01"] * 4,
            "new_2field": ["DRB1*03:01"] * 4,
            "rule": "test",
            "reason": "test",
        }

        with (
            patch.object(rescue, "read_tsv", return_value=rows),
            patch.object(rescue, "propose_class1_target90", return_value=None),
            patch.object(rescue, "propose_drb1", return_value=drb1_proposal) as propose_drb1,
            patch.object(rescue, "propose_dqb1_high_copy") as propose_dqb1,
            patch.object(rescue, "propose_dpb1") as propose_dpb1,
        ):
            _rows, proposals = rescue.proposals_for_sample(
                Path("/asm"), Path("/spechla"), "SIM0001", args, {}, {}
            )

        self.assertEqual(proposals, [drb1_proposal])
        propose_drb1.assert_called_once()
        propose_dqb1.assert_not_called()
        propose_dpb1.assert_not_called()


if __name__ == "__main__":
    unittest.main()