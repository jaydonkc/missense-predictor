import hashlib
import json
import unittest

from . import report as assay


class ReportTest(unittest.TestCase):
    def test_frozen_prediction_and_roles(self):
        first = assay.report("T279A")
        self.assertEqual(first["protein"], "ENSP00000377192.3")
        self.assertAlmostEqual(first["turnover_fraction"], 0.37456000654616783)
        self.assertFalse(first["training_residue"])
        self.assertEqual(first["predicted_assay_profile"], "low_activity__preserved_abundance")
        self.assertTrue(assay.report("F501C")["training_residue"])
        with self.assertRaisesRegex(ValueError, "No qualified prediction"):
            assay.report("T279*")

    def test_protocol_and_fit_binding(self):
        pin = json.loads((assay.HERE / "score-lock.json").read_bytes())
        for name, field in [("PROTOCOL.md", "protocol_sha256"), ("predict.py", "predictor_sha256"),
                            ("train_lock.py", "training_code_sha256")]:
            self.assertEqual(hashlib.sha256((assay.HERE / name).read_bytes()).hexdigest(), pin[field])


if __name__ == "__main__":
    unittest.main()
