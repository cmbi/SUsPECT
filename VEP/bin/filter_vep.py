#!/usr/bin/env python3

# A series of filters to reduce the number of variants 
# that should be run through SIFT/Polyphen2.
# 
# input: gff file used as input to VEP and VEP output file (vcf).
#   
# ouput: A filtered VEP vcf output file, an intermediate file listing the 
# variants that had a more severe annotation in the custom annot. 

#%%
import pandas as pd
import numpy as np
import argparse
import io, os
from collections import ChainMap
from vcf_parser import VCFParser

#%%
def get_transcript_ids(gff_file):
    gff=pd.read_table(gff_file,names=['seqname', 'source', 'feature', 'start','end', 'score', 'strand','frame','extra'])
    gff=gff[gff.feature=='mRNA'][['extra']]
    gff=gff[gff.extra.str.contains('_ORF_')] # make sure CPAT orfs
    gff['Transcripts']=gff.extra.str.split(';Parent').apply(lambda x: x[0])
    gff['Transcripts']=gff['Transcripts'].str.strip('ID=')
    return(gff[['Transcripts']])

def get_target_variants(path,transcripts):
    '''return variants of interest given novel transcripts'''
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    df=pd.read_csv(io.StringIO(''.join(lines)),dtype={'#CHROM': str, 'POS': str, 'ID': str, 'REF': str, 'ALT': str,'QUAL': str, 'FILTER': str, 'INFO': str},sep='\t').rename(columns={'#CHROM': 'CHROM'})
    df['Transcripts']=df.INFO.str.findall('Transcript\|(\w+)\|')
    df['variant_id']=df['CHROM']+'_'+df['POS']+'_'+df['REF']+'_'+df['ALT']
    df=df[['variant_id','Transcripts']].explode('Transcripts').merge(transcripts).drop_duplicates('variant_id')
    return(set(df.variant_id.unique()))

def read_vcf_vep(vcfparserobj,target_variants):
    '''to parse VCF formatted VEP output'''
    variant_info=pd.DataFrame()
    for variant in vcfparserobj:
        if variant['variant_id'] in target_variants:
            #if variant['MAX_AF']=='' or 'e' in variant['MAX_AF'] or float(variant['MAX_AF'])<0.01:
            l=[*variant['vep_info'].values()]
            l=[item for sublist in l for item in sublist]
            subdf=pd.DataFrame(l)
            subdf['Chrom']=variant['CHROM']
            subdf['Pos']=variant['POS']
            subdf['Ref']=variant['REF']
            subdf['Alt']=variant['ALT']
            subdf['Uploaded_variation']=variant['variant_id']
            variant_info=pd.concat([variant_info,subdf])
    variant_info['Impact_level']=variant_info['IMPACT'].apply(classify_impact)
    variant_info['MAX_AF']=variant_info['MAX_AF'].str.replace('','0') # novel variants given maxaf of 0
    variant_info['MAX_AF']=variant_info['MAX_AF'].apply(lambda x: 0 if 'e' in x else x) # this python package has a bug! cannot parse AF with exp.
    return(variant_info)

def classify_impact(impact):
    '''information about vep consequence classification can be found at https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
    '''
    if impact=='HIGH':
        return(3)
    elif impact=='MODERATE':
        return(2)
    elif impact=='LOW':
        return(1)
    elif impact=='MODIFIER':
        return(0)
    return()

def black_list(novel_dict,ref_list):
    black=[]
    for acc,aa in novel_dict.items():
        if aa in ref_list:
            black.append(acc)
    return(black)

def case_both_missense(df):
    missense_only=df[df.Amino_acids!='-']
    missense_only['aa_dict']=missense_only.apply(lambda x: {x['Feature']:x['Amino_acids']},axis=1)
    missense_only=missense_only.groupby(['Uploaded_variation','nov'])['aa_dict'].apply(list).unstack().dropna()
    missense_only['Novel']=missense_only['Novel'].apply(lambda a: dict(ChainMap(*a)))
    missense_only['Reference']=missense_only['Reference'].apply(lambda a: dict(ChainMap(*a)))
    missense_only['blacklist']=missense_only.apply(lambda x: black_list(x['Novel'],[*x['Reference'].values()]),axis=1)
    return(missense_only['blacklist'].explode().dropna().unique())

def print_filtered_vep(parser,good_dict,outfile):
    handle=open(outfile,'w')
    for line in parser.metadata.print_header():
        handle.write(line+'\n')
    for variant in parser:
        key=(variant['CHROM'],variant['POS'])
        if key in good_dict and list(variant['vep_info'].keys())[0] in list(good_dict[key]):
            handle.write('\t'.join([variant[head] for head in parser.header])+'\n')
    handle.close()

#%%
def main():
    aparser=argparse.ArgumentParser(description='Filters VEP VCF output for potentially interesting variants')
    aparser.add_argument('input_gff_file', metavar='INPUT_GFF', type=str,help='Input file gff file used as VEP input') #gff file used for vep - to get the names of the novel sequences
    aparser.add_argument('input_vep_file', metavar='INPUT_VCF', type=str,help='Input file in the VEP VCF output file format') #VEP output file (when using both cache and gff file as annotation sources)
    aparser.add_argument('output_file', metavar='OUTPUT', type=str,help='The output VCF file')
    args = aparser.parse_args()
    # read in data
    tids=get_transcript_ids(args.input_gff_file)
    good_vars=get_target_variants(args.input_vep_file,tids)
    parser= VCFParser(infile=args.input_vep_file, split_variants=True, check_info=True)
    parse_out=parser
    vep=read_vcf_vep(parser,good_vars)
    # remove variants that are only in reference sequences
    vep['nov']=vep.Feature.isin(tids) # True = novel seq
    variants_in_novel=vep[vep.nov]['Uploaded_variation'].unique() # get variants that fall into at least 1 novel transcript
    vep['nov'] = np.where(vep['nov'], 'Novel', 'Reference')
    vep_nov=vep[vep.Uploaded_variation.isin(variants_in_novel)] 
    # remove variants where (more) severe in reference set
    severity=vep_nov.groupby(['Uploaded_variation','nov'])['Impact_level'].apply(max).unstack().reset_index()
    severity[(severity['Novel']>severity['Reference'])&(severity['Reference']<2)][['Uploaded_variation']].to_csv('Variants_more_severe_in_new_annotation.csv',index=False,header=False) # intermediate list
    filter_severity=severity[(severity['Novel']>=severity['Reference'])&(severity['Reference']<3)]['Uploaded_variation'].unique() # less strict filtering
    vep_nov=vep_nov[vep_nov.Uploaded_variation.isin(filter_severity)] 
    # filter allele frequency - only interested in rare/novel variants
    vep_nov=vep_nov[vep_nov.MAX_AF.astype(float)<0.01]
    # when both ref and custom are missense,
    # remove variant effects in custom seqs 
    # where aa subs are same in ref & custom
    blacklist=case_both_missense(vep_nov) 
    vep_nov=vep_nov[~vep_nov.Feature.isin(blacklist)]
    greendict=vep_nov.groupby(['Chrom','Pos'])['Alt'].apply(set).to_dict()
    #write to output file
    print_filtered_vep(parse_out,greendict,args.output_file)
    # vep[vep.Feature.isin(greenlist)].drop(columns=['nov']).to_csv('filtered_vep.tab',header=False,index=False,sep='\t') # tabular format

def test():
    tids=get_transcript_ids('orf_prediction/testing/test/final_gtf/complete_sorted.gff3') # make sure to sed the first line first
    parser= VCFParser(infile='orf_prediction/testing/novelonly.vepout.vcf', split_variants=True, check_info=True)
    vep['Uploaded_variation']=vep['Chrom']+'_'+vep['Pos'].str+'_'+vep['Ref']+'/'+vep['Allele']

if __name__ == '__main__':
    main()

    