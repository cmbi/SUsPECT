#!/usr/bin/env python3

import pandas as pd
import argparse

'''This script only needs to be run once, after the creation of an abundance file by talon
this will then be input (--ab option) for hg38_vs_custom
'''

def convert_to_ratio(df):
    '''get the sum for that gene
    input: slice of the abundance dataframe that corresponds to one gene
    output: same slice with extra columns containing gene totals
    '''
    df['calb_total']=df['calb'].sum()
    df['rpmi_total']=df['rpmi'].sum()
    df['saur_total']=df['saur'].sum()
    return(df)

def main():
    '''
    input: talon abundance file
    output: talon abundance file with a few extra columns corresponding to totals
    '''
    parser = argparse.ArgumentParser(description='mapping of variants to allele frequency')
    parser.add_argument('--ab', help='path talon abundances', required=True)
    parser.add_argument('--out', help='output file name', required=True)
    args=vars(parser.parse_args())
    df=pd.read_csv(args['ab'],sep='\t')
    genes=df['annot_gene_name'].unique()
    out_df=pd.DataFrame()
    for g in genes:
        out_df=pd.concat([out_df,convert_to_ratio(df[df['annot_gene_name']==g])])
    out_df.to_csv(args['out'],sep='\t',index=False)

if __name__ == '__main__':
    main()