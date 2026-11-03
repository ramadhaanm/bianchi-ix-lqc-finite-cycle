"""Scientific regression checks for the revised statistical definitions."""
import unittest
import sys
import numpy as np
from reanalyze import interval, wilson, Run
from pathlib import Path


class ReanalysisTests(unittest.TestCase):
    def test_single_pair_does_not_claim_zero_sampling_uncertainty(self):
        result=interval([4.3])
        self.assertEqual(result['median'],4.3)
        self.assertIsNone(result['lo'])
        self.assertIsNone(result['hi'])

    def test_zero_counts_are_upper_limits_not_pseudo_probabilities(self):
        lo,hi=wilson(0,400)
        self.assertAlmostEqual(lo,0.)
        self.assertGreater(hi,0.)
        self.assertLess(hi,.01)

    def test_median_of_total_is_not_sum_of_medians(self):
        # Three physical pairs whose phase gains telescope exactly.
        phases=np.array([[0.,0.,9.],[0.,9.,0.]])
        direct=interval(phases.sum(axis=0))['median']
        wrong=sum(interval(p)['median'] for p in phases)
        self.assertEqual(direct,9.)
        self.assertEqual(wrong,0.)

    def test_small_signal_survives_hilo_difference(self):
        # The hi components coincide, but the separation is real at Double64 precision.
        a=np.array([20.,2e-25,30.,3e-25,40.,4e-25])
        b=np.array([20.,1e-25,30.,1e-25,40.,1e-25])
        diff=(a[:6:2]-b[:6:2])+(a[1:6:2]-b[1:6:2])
        np.testing.assert_allclose(diff,[1e-25,2e-25,3e-25],atol=1e-40)

    def test_physical_exit_never_produces_later_events(self):
        path=Path(__file__).resolve().parents[1]/'results/revised/phase_L20_K_N4_on0.5_b300000.csv'
        run=Run(path)
        for row in run.rows.values():
            reached=np.isfinite(row['a'][:,0])
            if not reached.all():
                first_missing=np.flatnonzero(~reached)[0]
                self.assertFalse(reached[first_missing:].any())


if __name__=='__main__':
    unittest.main(testRunner=unittest.TextTestRunner(stream=sys.stdout))
