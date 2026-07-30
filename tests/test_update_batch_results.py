import csv
import tempfile
import unittest
from pathlib import Path

import update_batch_results as batch


class BatchResultsTest(unittest.TestCase):
    def test_parse_metadata_rejects_unknown_keys(self):
        with self.assertRaisesRegex(ValueError, "unsupported metadata key"):
            batch.parse_metadata(["unknown=value"])

    def test_atomic_write_and_replace_key(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "results.tsv"
            base = {field: "" for field in batch.FIELDS}
            first = {
                **base,
                "experiment": "v4",
                "condition": "graft10_cov300x",
                "sample": "SIM0001",
                "gene": "HLA-A",
                "allele_multiset": "A*01:01,A*02:01,A*03:01,A*04:01",
            }
            batch.write_atomic(output, [first])
            rows = batch.read_tsv(output)
            self.assertEqual(rows, [first])

            replacement = {**first, "allele_multiset": "A*11:01,A*12:01,A*13:01,A*14:01"}
            key = ("v4", "graft10_cov300x", "SIM0001", "HLA-A")
            combined = [
                row for row in rows
                if (row["experiment"], row["condition"], row["sample"], row["gene"]) != key
            ]
            combined.append(replacement)
            batch.write_atomic(output, combined)

            with output.open() as handle:
                written = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(written), 1)
            self.assertEqual(written[0]["allele_multiset"], replacement["allele_multiset"])


if __name__ == "__main__":
    unittest.main()