#!/usr/bin/env python3

import pandas as pd
import argparse
import io

'''
This script parses the tab output of VEP and creates input file(s) for external pathogenicity prediction
'''

def read_vep(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    df=pd.read_csv(io.StringIO(''.join(lines)),sep='\t').rename(columns={'#Uploaded_variation': 'Uploaded_variation'})
    df=df[(df['Feature'].str.contains('novel'))&(df['Consequence']=='missense_variant')]
    df['Amino_acids']=df['Amino_acids'].str.split('/')
    df[['old','new']]=pd.DataFrame(df['Amino_acids'].tolist(),index=df.index)
    df['rare']=df['Extra'].apply(af_filter)
    rare,common=df[df['rare']],df[~df['rare']]
    return(rare[['Feature','Protein_position','old','new']],common[['Feature','Protein_position','old','new']])

def af_filter(info):
    '''given the information, separate rare/novel from common'''
    if 'MAX_AF' in info:
        if float(info.split('MAX_AF=')[1].split(';',1)[0])<0.01:
            return(True)
        else:
            return(False)
    else:
        return(True)

def main():
    parser = argparse.ArgumentParser(description='VEP parser for polyphen-2')
    parser.add_argument('--vep', help='vep output', required=True)
    args=vars(parser.parse_args())
    rare,common=read_vep(args['vep'])
    rare.drop_duplicates().to_csv('polyphen_in_rare.input',sep='\t',header=False,index=False)
    common.drop_duplicates().to_csv('polyphen_in_common.input',sep='\t',header=False,index=False)
    
if __name__ == '__main__':
    main()