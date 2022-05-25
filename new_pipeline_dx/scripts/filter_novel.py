#!/usr/bin/env python3

import pandas as pd
import gtfparse

def main():
    parser = argparse.ArgumentParser("IO file locations for make pacbio cds gtf")
    parser.add_argument("--gffcompare_gtf", action="store", dest = "gffcompare_gtf")
    results = parser.parse_args()
    combined=gtfparse.read_gtf(results.gffcompare_gtf)
    exclude=combined[combined['class_code']=='=']['transcript_id'].unique()
    combined=combined[~combined['transcript_id'].isin(exclude)]
    df=pd.read_table(results.gffcompare_gtf, names=['seqname','source','feature','start','end','score','strand','frame','extra'])
    df.merge(combined[['seqname','source','feature','start','end']]).to_csv('gffcmp.combined.filtered.gtf',sep='\t',index=False,header=False)

if __name__ == "__main__":
    main()