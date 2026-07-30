import tempfile
import unittest
from pathlib import Path

from benchmark.scripts.simulate_wgsim_benchmark import cleanup_caller_intermediates


class SimulationCleanupTest(unittest.TestCase):
    def test_removes_regenerable_files_and_preserves_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sample_root = root / "work" / "spechla_out" / "SIM0001"
            em_root = sample_root / "em_refine"
            drb345_root = sample_root / "drb345"
            em_root.mkdir(parents=True)
            drb345_root.mkdir()
            removable = [
                sample_root / "SIM0001.map_database.bam",
                sample_root / "A.R1.fq.gz",
                em_root / "HLA-A.aug.fa.bwt",
                drb345_root / "HLA-DRB345.aug.fa.sa",
                root / "SIM0001.R1.fastq.gz",
                root / "SIM0001.R2.fastq.gz",
            ]
            evidence = em_root / "HLA-A.tf_counts.tsv"
            for path in removable + [evidence]:
                path.touch()

            removed = cleanup_caller_intermediates(
                root / "work", "SIM0001", removable[-2], removable[-1]
            )

            self.assertEqual(removed, len(removable))
            self.assertTrue(evidence.exists())
            self.assertFalse(any(path.exists() for path in removable))


if __name__ == "__main__":
    unittest.main()