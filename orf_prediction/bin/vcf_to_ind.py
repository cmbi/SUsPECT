#!/usr/bin/env python3

import pandas as pd
import argparse
import io

def combined_vcf_to_mapping(path):
    '''take the combined vcf file and make a 2 column mapping file of variants to patient IDs who have the variant
    can use this as alternative to '--individual all' on the VEP command line, which prints one line for every patient.
    '''
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    df=pd.read_csv(io.StringIO(''.join(lines)),sep='\t').rename(columns={'#CHROM': 'CHROM'}).drop(columns=['ID','QUAL','FILTER','INFO','FORMAT'])
    df=df.set_index(['CHROM','POS','REF','ALT']).replace(':.*','',regex=True).replace({'./.':None}).apply(lambda x : x.name + ':' + x if x is not None else x)
    df['patients']=df.apply(lambda x: ';'.join(x.dropna().astype(str)),axis=1)
    df[['patients']].reset_index().drop_duplicates().to_csv('variants_to_patients.tsv',sep='\t',index=False)
    #pd.DataFrame(df.notna().dot(df.columns+',').str.rstrip(','), index=df.index, columns=['patients'])

def main():
    parser = argparse.ArgumentParser(description='VEP parser for polyphen-2')
    parser.add_argument('vcf', help='vep output')
    args=parser.parse_args()
    combined_vcf_to_mapping(args.vcf)

if __name__ == '__main__':
    main()