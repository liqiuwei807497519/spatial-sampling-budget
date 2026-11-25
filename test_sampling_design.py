import unittest
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd

from sampling_design import POLICIES,draws,recover,repeat_plan,compare_blocks,calibration

class PlanningChecks(unittest.TestCase):
    def test_recovery(self):
        pool=list("aabbxxx")
        for k in range(1,8):
            draws_=list(combinations(pool,k))
            expected=sum("a" in d and "b" in d for d in draws_)/len(draws_)
            self.assertAlmostEqual(recover(7,2,2,k),expected)
    def test_repeats(self):
        self.assertIsNone(repeat_plan(0,.8))
        self.assertEqual(repeat_plan(1,.95),1)
        self.assertEqual(repeat_plan(.044763785,.8),36)
        self.assertEqual(repeat_plan(.044763785,.8,.10),443)
    def test_designs(self):
        frame=pd.read_csv(Path(__file__).parent/"example_sites.csv")
        for policy in POLICIES:
            ids,w,ok=draws(frame,12,50,12,policy)
            self.assertTrue((np.diff(np.sort(ids,axis=1))>0).all())
            np.testing.assert_allclose(w[ok].sum(axis=1),0,atol=1e-10)
    def test_overlap(self):
        frame=pd.read_csv(Path(__file__).parent/"example_sites.csv")
        blocks=pd.DataFrame({"block":["a","b"],"zone":["a","b"],"sites":["s0|s1","s1|s2"]})
        r=compare_blocks(frame,blocks,["Acan","Col2a1"],50)
        self.assertEqual((len(r),r.fit.sum(),r.strict.sum()),(1,0,0))

if __name__=="__main__":unittest.main()
