import unittest

from iterative_remap_em import rescue_class_i_distinct


class ClassIDistinctRescueTest(unittest.TestCase):
    def setUp(self):
        self.counts = {
            "A*01:01": 450.0,
            "A*02:01": 430.0,
            "A*03:01": 55.0,
            "A*24:02": 45.0,
            "A*26:01": 20.0,
        }
        self.winners = ("A*01:01", "A*02:01", "A*03:01", "A*03:01")

    def test_replaces_one_duplicate_when_support_is_separated(self):
        rescued, detail, residual = rescue_class_i_distinct(
            self.counts, self.winners, 0.10, min_gap=1.5
        )
        self.assertEqual(set(rescued), {"A*01:01", "A*02:01", "A*03:01", "A*24:02"})
        self.assertIsNotNone(detail)
        self.assertIsNotNone(residual)

    def test_keeps_quartet_when_fourth_and_fifth_are_not_separated(self):
        counts = dict(self.counts, **{"A*24:02": 36.0, "A*26:01": 30.0})
        rescued, detail, residual = rescue_class_i_distinct(
            counts, self.winners, 0.10, min_gap=1.5
        )
        self.assertEqual(rescued, self.winners)
        self.assertIsNone(detail)
        self.assertIsNone(residual)

    def test_keeps_quartet_with_two_distinct_alleles(self):
        winners = ("A*01:01", "A*01:01", "A*02:01", "A*02:01")
        rescued, detail, residual = rescue_class_i_distinct(
            self.counts, winners, 0.10, min_gap=1.5
        )
        self.assertEqual(rescued, winners)
        self.assertIsNone(detail)
        self.assertIsNone(residual)


if __name__ == "__main__":
    unittest.main()