#!/usr/bin/env python3

import pandas as pd
import gtfparse
import sys

def main():
    combined=gtfparse.read_gtf(sys.argv[1])
    exclude=combined[combined['class_code']=='=']['transcript_id'].unique()
    combined=combined[~combined['transcript_id'].isin(exclude)]
    df=pd.read_table(sys.argv[1], names=['seqname','source','feature','start','end','score','strand','frame','extra'])
    df.merge(combined[['seqname','source','feature','start','end']]).to_csv('gffcmp.combined.filtered.gtf',sep='\t',index=False,header=False)

if __name__ == "__main__":
    main()