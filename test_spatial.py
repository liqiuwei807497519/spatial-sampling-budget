from pathlib import Path
import unittest
import numpy as np
import pandas as pd
from design_atlas import policies,evaluate
from sampling_design import benchmark

class SpatialChecks(unittest.TestCase):
    def setUp(self):
        y,x=np.mgrid[:6,:8]
        self.f=pd.DataFrame(dict(x=x.ravel(),y=y.ravel(),site=np.arange(x.size),stratum=np.where(y.ravel()<3,'dark','light'),total=100,Acan=1,Col2a1=1))
    def test_spatial_baselines(self):
        for method in ['maximin','inhibition']:
            ids,w,ok=policies(self.f,12,32,7,method)
            self.assertTrue(np.all(np.diff(np.sort(ids,axis=1))>0))
            np.testing.assert_array_equal(ids,policies(self.f,12,32,7,method)[0])
            np.testing.assert_allclose(w[ok].sum(1),0,atol=1e-14)
    def test_dual_vs_pooled_events(self):
        # One molecule per site, with complementary single-target donor sites.
        f=pd.DataFrame(dict(total=[2,2],Acan=[1,0],Col2a1=[0,1],stratum=['dark','light']))
        vals=evaluate(f,np.array([[0,1],[0,1]]),np.array([[1.,-1.],[1.,-1.]]),np.array([True,True]),1,{'x':('Acan','Col2a1')})[0]
        self.assertAlmostEqual(vals['dual_probability'],.25)
        self.assertAlmostEqual(vals['rmse_per10k'],10000*np.sqrt(.5))
    def test_synthetic_frame(self):
        base=Path(__file__).parent
        frame=pd.read_csv(base/'synthetic_pilot.csv').query("section=='Synthetic-A' and spatial_cell>=0")
        result=benchmark(frame,['Acan','Col2a1'],12,50,120,42,['Acan','Col2a1'])
        self.assertEqual(len(result),15)
        self.assertTrue(np.isfinite(result.rmse).all())

if __name__=='__main__':unittest.main()
