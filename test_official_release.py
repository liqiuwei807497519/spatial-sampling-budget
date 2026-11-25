"""Full archive integration test. Set OMIX_ARCHIVES to the official five-file directory."""
import os,tempfile,unittest
from pathlib import Path
import pandas as pd
from omix_readout import score_inputs
from ngdc_files import mapped_inputs
from check_public import verify
class OfficialRelease(unittest.TestCase):
    @unittest.skipUnless(os.environ.get('OMIX_ARCHIVES'),'set OMIX_ARCHIVES for full official-file integration test')
    def test_five_archives_match_frozen_endpoints(self):
        base=Path(__file__).parent
        mapped=mapped_inputs(Path(os.environ['OMIX_ARCHIVES']))
        self.assertEqual(len({p for _,p in mapped}),5)
        self.assertEqual(len({s for s,_ in mapped}),5)
        axis=pd.read_csv(base/'reference_21479_to_20512.csv').query('eligible==1').set_index('stable_gene_id').symbol
        with tempfile.TemporaryDirectory() as dr:
            dest=Path(dr);score_inputs(mapped,axis,dest);verify(dest)
if __name__=='__main__':unittest.main()
