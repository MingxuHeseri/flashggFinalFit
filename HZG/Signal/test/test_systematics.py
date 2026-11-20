"""Run with Python in the HZG signal environment."""
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

HZG = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(HZG / "tools"), str(HZG / "Signal" / "tools")]
from finalModel import FinalModel


class SystematicsTest(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.path = Path(directory.name) / "cat.pkl"

    def load(self, rows):
        pd.DataFrame(rows).to_pickle(self.path)
        coefficients = {}
        def record(name, mean, sigma, rate, kind, opts):
            coefficients[name] = (mean, sigma, rate)
        model = SimpleNamespace(procSyst="ggH", cat="cat", flav="mu", sqrts="13TeV",
                                scales="energy", scalesCorr="", scalesGlobal="", smears="",
                                NuisanceMap={}, makeNuisance=record)
        FinalModel.buildNuisanceMap(model, self.path)
        return coefficients["energy"]

    def test_flavor_and_qualified_columns(self):
        ele = dict(proc="ggH", cat="cat", flav="ele", energy_13TeVscale_mean=0.1,
                   energy_13TeVscale_sigma=0.2, energy_13TeVscale_rate=0.3)
        mu = dict(ele, flav="mu", energy_13TeVscale_mean=0.4)
        self.assertEqual(self.load([ele, dict(mu, cat="other"), mu]), (0.4, 0.2, 0.3))

    def test_shared_row_and_invalid_matches(self):
        row = dict(proc="ggH", energy_mean=0.1, energy_sigma=0.2, energy_rate=0.3)
        self.assertEqual(self.load([row]), (0.1, 0.2, 0.3))
        for rows in ([dict(row, proc="VH")], [row, row], [dict(row, energy_rate=float("nan"))]):
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                self.load(rows)


if __name__ == "__main__":
    unittest.main()
