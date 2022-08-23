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
    df=pd.read_csv(io.StringIO(''.join(lines)),sep='\t').rename(columns={'#CHROM': 'CHROM'})
    df=df.drop(columns=['ID','QUAL','FILTER','INFO','FORMAT']).set_index(['CHROM','POS','REF','ALT']).replace({'./.:.:.:.:.':None})
    pd.DataFrame(df.notna().dot(df.columns+',').str.rstrip(','), index=df.index, columns=['patients']).reset_index().drop_duplicates().to_csv('variants_to_patients.tsv',sep='\t',index=False)

def main():
    parser = argparse.ArgumentParser(description='VEP parser for polyphen-2')
    parser.add_argument('vcf', help='vep output')
    args=parser.parse_args()
    combined_vcf_to_mapping(args.vcf)

if __name__ == '__main__':
    main()