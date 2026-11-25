"""Retrieve or reuse all five official archives, score their designated samples and check both endpoints."""
import argparse
from pathlib import Path
import pandas as pd
from ngdc_files import mapped_inputs
from omix_readout import score_inputs
from check_public import verify

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--archives',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--offline',action='store_true',help='Require all five verified files to be available locally')
    a=ap.parse_args();base=Path(__file__).parent
    ref=pd.read_csv(base/'reference_21479_to_20512.csv').query('eligible==1').set_index('stable_gene_id').symbol
    assert len(ref)==20512 and ref.index.is_unique
    if a.archives.is_file():ap.error('--archives is the directory of the five official files')
    score_inputs(mapped_inputs(a.archives,not a.offline),ref,a.out)
    verify(a.out)
if __name__=='__main__':main()
