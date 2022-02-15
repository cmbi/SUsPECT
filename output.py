#!/usr/bin/env python3


import pandas as pd
import argparse
import re, io

'''
integrates multiple outputs with each other to create a candidate variant file to be analyzed further in a jupyter notebook.
Adds transcript/exon abundance, pph scores, patients to VEP output for final variant prioritization
'''

def read_vep(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    df=pd.read_csv(io.StringIO(''.join(lines)),sep='\t').rename(columns={'#Uploaded_variation': 'Uploaded_variation'})
    # df=df[((df.Extra.str.contains('IMPACT=MODERATE'))|(df.Extra.str.contains('IMPACT=HIGH')))&(df.Feature.str.contains('novelT'))]# filter for more severe in new annotation
    df['tmp']=df.Extra.apply(mine_info)
    df[['Consequence_level','max_af','conservation']] = pd.DataFrame(df.tmp.tolist(), index= df.index)
    return(df.drop(columns=['tmp']).drop_duplicates(subset=['Uploaded_variation','Feature']))

def simple_vep(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    df=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',usecols=[0,1,2,6,13]).rename(columns={'#Uploaded_variation': 'Uploaded_variation'})
    df['Consequence_level']=df.Extra.apply(classify_consequence)
    df['pp_hg38']=df.Extra.apply(extract_polyphen)
    return(df.drop(columns=['Extra']).sort_values('Consequence_level',ascending=False).groupby('Uploaded_variation').first())
    
def extract_af(s):
    for p in s.split(';'):
        if 'MAX_AF=' in p:
            return(re.sub('MAX_AF=','',p))

def extract_conservation(s):
    for p in s.split(';'):
        if 'Conservation=' in p:
            return(re.sub('Conservation=','',p))

def extract_polyphen(extra):
    if 'PolyPhen=' in extra:
        extra=extra.split('PolyPhen=')[1]
        return(extra.split('(')[0])
    return(None)

def classify_consequence(extra):
    '''information about vep consequence classification can be found at https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
    '''
    if 'IMPACT=' in extra:
        extra=extra.split('IMPACT=')[1]
        extra=extra.split(';',1)[0]
        if 'MODERATE' in extra:
            return(2)
        elif 'HIGH' in extra:
            return(3)
        elif 'LOW' in extra:
            return(1)
    return(0) #modifier

def mine_info(extra):
    info=[]
    maxaf=''
    conservation=''
    for p in extra.split(';'):
        if 'IMPACT=' in p:
            if 'MODERATE' in p:
                info.append(2)
            elif 'HIGH' in p:
                info.append(3)
            elif 'LOW' in p:
                info.append(1)
            else:
                info.append(0)
        elif 'MAX_AF=' in p:
            maxaf=re.sub('MAX_AF=','',p)
        elif 'Conservation=' in p:
            conservation=re.sub('Conservation=','',p)
    info.extend([maxaf,conservation])
    assert len(info)==3, 'Info not collected properly'
    return(info)

def polyphen_read(pph):
    df=pd.read_csv(pph,sep='\t').drop_duplicates()
    df.columns=df.columns.str.strip()
    df['Amino_acids']=df['o_aa1'].str.strip()+'/'+df['o_aa2'].str.strip()
    df['Feature']=df['#o_acc'].str.strip()
    df['Protein_position']=df.o_pos.astype(str)
    return(df[['Amino_acids','Feature','Protein_position','prediction', 'based_on', 'effect', 'site', 'region']])

def main():
    parser = argparse.ArgumentParser(description='Cross reference VEP output file with transcript abundance and other relevant information')
    parser.add_argument('--vep', help='hg38+pacbio vep annotated file (only most severe)', required=True)
    parser.add_argument('--ab', help='abundance file made by transcript_dominance.py', required=True)
    parser.add_argument('--pp', help='polyphen output', required=False)
    parser.add_argument('--hg38', help='VEP run with hg38 only, to compare annotations', required=False)
    parser.add_argument('--xn', help='exon-level abundances made by variants_in_novel_edges.py', required=True)
    parser.add_argument('--go', help='go terms', required=False)
    parser.add_argument('--pat', help='patient mapping created by function in vcf_to_ind.py', required=False)
    parser.add_argument('--out', help='Outfile name', required=True)
    args=vars(parser.parse_args())
    #read in the files
    if args['pat']:
        combi=read_vep(args['vep']).merge(pd.read_csv(args['pat']),on='Uploaded_variation')
    else:
        combi=read_vep(args['vep'])
    combi['transcript']=combi['Feature'].apply(lambda x: x.split('.')[0] if '.' in x else x)
    abund=pd.read_csv(args['ab'],sep='\t')
    abund['transcript']=abund['annot_transcript_id'].apply(lambda x: x.split('.')[0] if '.' in x else x)
    # be aware the next step removes variants that occur in transcripts that are not found in the transcriptome
    combi=combi.merge(abund,on='transcript').merge(pd.read_csv(args['xn']).rename({'calb':'calb_exon','rpmi':'rpmi_exon','saur':'saur_exon'},axis=1),how='left')
    if args['pp']:
        combi=combi.merge(polyphen_read(args['pp']),on=['Feature','Amino_acids','Protein_position'],how='left')
    if args['go']:
        combi['Ensembl_gene_identifier']=combi.annot_gene_id.apply(lambda x: x.split('.')[0])
        go=pd.read_csv(args['go'],sep='\t',usecols=['Ensembl_gene_identifier','GO_ID','GO_term','Category'])
        go=go[go.Category=='Process'].groupby('Ensembl_gene_identifier').agg({'GO_ID':lambda x: ','.join(list(x)),'GO_term':lambda x: ','.join(list(x))}).reset_index()
        combi=combi.merge(go,how='left')
    if args['hg38']:
        combi=combi.merge(simple_vep(args['hg38']),on=['Uploaded_variation','Location','Allele'],suffixes=('','_hg38'))   
    combi.to_csv(args['out'],index=False)

def twofile():
    '''optional filterings'''
    l=['benign','unknown', None]
    df[(df.prediction.str.contains('damaging'))&(df.pp_hg38.isin(l))&(df.Consequence_hg38=='missense_variant')].to_csv('twofile_more_severe_missense.csv',index=False)
    df[(df.Consequence_level!=0)|(df.Consequence_level_hg38!=0)].to_csv('twofile_analysis_nonzero.csv',index=False)
    
if __name__ == '__main__':
    main()