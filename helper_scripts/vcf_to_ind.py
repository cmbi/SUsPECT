#!/usr/bin/env python3

import pandas as pd
import argparse
import io

def parse(file):
    df=pd.DataFrame()
    with open(file,'r') as f:
        for line in f:
            if not line.startswith('#'):
                start=line.split('PASS')[0]
                start=start.split('\t')
                ident='_'.join([start[0],start[1],start[3]])+'/'+start[4]
                # count=total_patients-line.count("./.:.:.:.:.")
                count_hom=line.count("1/1:")
                count_het=line.count("0/1:")
                df = df.append({'Uploaded_variation': ident, 'patients_hom':count_hom,'patients_het':count_het}, ignore_index=True)
    return(df)

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
    # df
 

def main():
    parse()