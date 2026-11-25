import unittest,itertools
import numpy as np,pandas as pd
from resource_analysis import nondominated,ht_weights
from sampling_design import partition
class ResourceChecks(unittest.TestCase):
    def test_known_population_totals(self):
        f=pd.DataFrame(dict(x=np.arange(6),y=[0]*6,stratum=['dark','light','light','dark','dark','light']))
        val=np.array([0,1,2,5,1,3],float)
        bins=partition(np.arange(6),f[['x','y']].to_numpy(),2)
        draws=np.array(list(itertools.product(*bins)))
        est=(ht_weights(f,draws)*val[draws]).sum(1)
        truth=val[f.stratum.eq('dark')].mean()-val[f.stratum.eq('light')].mean()
        self.assertAlmostEqual(est.mean(),truth)
    def test_resource_dominance(self):
        f=nondominated([[3,-.1,1],[3,-.2,1],[2,-.3,3],[3,-.2,1],[4,-.1,1]])
        self.assertEqual(f.tolist(),[False,True,True,True,False])
