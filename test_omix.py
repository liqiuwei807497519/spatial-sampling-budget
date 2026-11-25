import gzip,math,tempfile,unittest
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from ngdc_files import records,verify
from omix_readout import score_inputs,no_event_mean
from check_public import elementary

class OmixTests(unittest.TestCase):
    def test_patient_map(self):
        self.assertEqual([r['library'] for r in records()],['IVD4','IVD1','IVD5','IVD2','IVD3'])
        self.assertEqual([int(r['matrix_members']) for r in records()],[1,5,1,1,1])
    def test_bad_digest(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'f';p.write_bytes(b'x')
            with self.assertRaises(ValueError):verify(p,dict(bytes=1,sha256='invalid'))
    def test_elementary(self):
        from itertools import combinations
        a=[.1,.2,.7,.85,.99]
        brute=np.mean([np.prod(q) for q in combinations(a,3)])
        self.assertAlmostEqual(elementary(a,3),brute)
        self.assertAlmostEqual(no_event_mean(a,3),brute)
    def test_assigned_directory_only_and_separate_targets(self):
        # Each of twelve complete barcodes contains exactly one of two targets.
        # Extra directories in the container must not enter the selected sample.
        with tempfile.TemporaryDirectory() as d:
            base=Path(d);p=base/'fixture.zip'
            feat='ENSG00000157766\tACAN\tGene Expression\nENSG00000139219\tCOL2A1\tGene Expression\nx\tX\tGene Expression\n'
            cols=''.join(f'BC{i}\n' for i in range(1,13))
            mtx='%%MatrixMarket matrix coordinate integer general\n3 12 24\n'+''.join(f'{1+i%2} {i} 1\n3 {i} 521\n' for i in range(1,13))
            with ZipFile(p,'w') as z:
                for sample in ['IVD1','extra']:
                    for name,val in [('features.tsv.gz',feat),('barcodes.tsv.gz',cols),('matrix.mtx.gz',mtx)]:z.writestr(f'{sample}/matrix/{name}',gzip.compress(val.encode()))
            axis=pd.Series(['ACAN','COL2A1','X'],index=['ENSG00000157766','ENSG00000139219','x'])
            score_inputs([('IVD1',p)],axis,base/'out')
            df=pd.read_csv(base/'out/OMIX017533_suspension.csv')
            self.assertEqual(list(df.library),['IVD1']);self.assertEqual(df.released.item(),12)
            self.assertAlmostEqual(df.within_barcode_strict.item(),0)
            self.assertAlmostEqual(df.strict_dual.item(),(1-(1/522)**6)**2)
            with self.assertRaises(ValueError):score_inputs([('IVD1',p),('IVD1',p)],axis,base/'out')
if __name__=='__main__':unittest.main()
