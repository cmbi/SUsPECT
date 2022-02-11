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
    df['Uploaded_variation']=df['CHROM']+'_'+df['POS'].astype(str)+'_'+df['REF']+'/'+df['ALT'].apply(lambda x: x.replace(',','/') if ',' in x else x)
    df=df.drop(columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']).set_index('Uploaded_variation').replace({'./.:.:.:.:.':None})
    pd.DataFrame(df.notna().dot(df.columns+',').str.rstrip(','), index=df.index, columns=['patients']).reset_index().to_csv('variants_to_patients.csv',index=False)

def main():
    parser = argparse.ArgumentParser(description='VEP parser for polyphen-2')
    parser.add_argument('--vcf', help='vep output', required=True)
    args=vars(parser.parse_args())
    combined_vcf_to_mapping(args['vcf'])

if __name__ == '__main__':
    main()