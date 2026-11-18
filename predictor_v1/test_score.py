import json
from pathlib import Path
import unittest

from .score import Predictor


class FrozenScoringTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scorer = Predictor()
        cls.examples = Path(__file__).parent / "examples"

    def test_archive_replay(self):
        expected = json.loads((self.examples / "expected.json").read_bytes())
        for route, item in expected.items():
            with self.subTest(route=route):
                request = json.loads((self.examples / (item["example"] + ".json")).read_bytes())
                actual = self.scorer.score(request)
                self.assertEqual(actual["variant_id"], item["variant"])
                self.assertEqual(actual["route"], route)
                self.assertAlmostEqual(actual["score"], item["score"], places=12)
                self.assertEqual(actual["screen_positive"], item["score"] >= actual["cutoff"])

    def test_reference_mismatch(self):
        req = json.loads((self.examples / "aligned_protein.json").read_bytes())
        req["ref"] = "G" if req["ref"] != "G" and req["alt"] != "G" else "A"
        with self.assertRaises(ValueError):
            self.scorer.score(req)

    def test_incomplete_alignment(self):
        req = json.loads((self.examples / "aligned_protein.json").read_bytes())
        del req["protein"]["alignment"]["log_ref_alt"]
        with self.assertRaises(ValueError):
            self.scorer.score(req)

    def test_full_context_rejects_unqualified_fold(self):
        req = json.loads((self.examples / "full_context.json").read_bytes())
        req["full_context"]["held_fold"] = -1
        with self.assertRaises(ValueError):
            self.scorer.score(req)


if __name__ == "__main__":
    unittest.main()
