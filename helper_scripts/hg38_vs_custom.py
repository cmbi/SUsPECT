#!/usr/bin/env python3


import pandas as pd
import argparse
import re, io

'''
integrates multiple outputs with each other to create a candidate variant file to be analyzed further in a jupyter notebook.
Adds transcript/exon abundance, pph scores, patients to VEP output
**this version is designed to take input from '--individuals all' on the command line.

example usage:
python hg38_vs_custom.py --vep all_cond_annot.conservation.tab --pat ../../hg38/liftover_vcfs/variants_to_patients.csv --ab /mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/talon/all_talon_abundance_filtered_genetotals.tsv --go /mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/ensembl2go.human.tab --xn variant_exon_presence.csv --pp pph_output_all.features --out onefile_hg38vcustom.csv
python hg38_vs_custom.py --vep twofile_all_bare.csv --pat ../../hg38/liftover_vcfs/variants_to_patients.csv --ab /mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/talon/all_talon_abundance_filtered_genetotals.tsv --go /mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/ensembl2go.human.tab --xn variant_exon_presence.csv --pp pph_output_all.features --out twofile_all_filled.csv
'''

def read_vep(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    df=pd.read_csv(io.StringIO(''.join(lines)),sep='\t').rename(columns={'#Uploaded_variation': 'Uploaded_variation'})
    # df=df[((df.Extra.str.contains('IMPACT=MODERATE'))|(df.Extra.str.contains('IMPACT=HIGH')))&(df.Feature.str.contains('novelT'))]# filter for more severe in new annotation
    df['tmp']=df.Extra.apply(mine_info)
    df[['Consequence_level','max_af','conservation']] = pd.DataFrame(df.tmp.tolist(), index= df.index)
    # patients=df.groupby(['Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature']).agg({'patient':lambda x: ','.join(list(x))}).rename({'patient':'patients'},axis=1).reset_index()
    return(df.drop(columns=['tmp']).drop_duplicates(subset=['Uploaded_variation','Feature']))

def extract_af(s):
    for p in s.split(';'):
        if 'MAX_AF=' in p:
            return(re.sub('MAX_AF=','',p))

def extract_conservation(s):
    for p in s.split(';'):
        if 'Conservation=' in p:
            return(re.sub('Conservation=','',p))


def classify_consequence(consequence):
    '''information about vep consequence classification can be found at https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
    '''
    if consequence in ['transcript_ablation','splice_acceptor_variant','splice_donor_variant','stop_gained','frameshift_variant','stop_lost','start_lost','transcript_amplification']:
        return(3) #high
    elif consequence in ['inframe_insertion','inframe_deletion','missense_variant','protein_altering_variant','	regulatory_region_ablation']:
        return(2) #moderate
    elif consequence in ['splice_region_variant','incomplete_terminal_codon_variant','start_retained_variant','stop_retained_variant','synonymous_variant']:
        return(1) #low
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
    # parser.add_argument('--np', help='num patients per variant', required=True)
    parser.add_argument('--pp', help='polyphen output', required=False)
    parser.add_argument('--xn', help='exon-level abundances made by variants_in_novel_edges.py', required=True)
    parser.add_argument('--go', help='go terms', required=False)
    parser.add_argument('--pat', help='patient mapping created by function in vcf_to_ind.py (hg38/liftover_vcfs/variants_to_patients.csv)', required=True)
    parser.add_argument('--out', help='Outfile name', required=True)
    args=vars(parser.parse_args())
    #read in the files
    # combi=read_vep(args['vep']).merge(pd.read_csv(args['pat']),on='Uploaded_variation')
    combi=pd.read_csv(args['vep']).merge(pd.read_csv(args['pat']),on='Uploaded_variation') # use alternative() function first
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
    combi.to_csv(args['out'],index=False)

def alternative():
    '''this is to compare the annotations in a different way: 
    1 file for output annotation with hg38 and 1 file for transcriptome-only annotation
    this way you can see not only what was more severe but what was not more severe and get other information
    '''
    def simple_vep(path):
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        df=pd.read_csv(io.StringIO(''.join(lines)),sep='\t').rename(columns={'#Uploaded_variation': 'Uploaded_variation'})
        df['Consequence_level']=df.Extra.apply(classify_consequence)
        return(df.sort_values('Consequence_level',ascending=False).groupby('Uploaded_variation').first())
    def classify_consequence(consequence):
        '''information about vep consequence classification can be found at https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
        '''
        if 'IMPACT=HIGH' in consequence:
            return(3) #high
        elif 'IMPACT=MODERATE' in consequence:
            return(2) #moderate
        elif 'IMPACT=LOW' in consequence:
            return(1) #low
        return(0) #modifier
    def extract_polyphen(extra):
        if 'PolyPhen=' in extra:
            extra=extra.split('PolyPhen=')[1]
            return(extra.split('(')[0])
        return(None)
    combi=read_vep('/mnt/xomics/renees/data/host_pathogen_PID/unsolved_cases/hg38/vep_builtin/all_patients_148.annotated.extended.conservation.vcf').merge(simple_vep('/mnt/xomics/renees/data/host_pathogen_PID/unsolved_cases/pacbio_anno/talon/observed_transcriptome_3cond_annot.tab'),on=['Uploaded_variation','Location','Allele'],suffixes=('_hg38',''))
    combi['pp_hg38']=combi.Extra_hg38.apply(extract_polyphen)
    combi[combi.Feature!='-'].to_csv('twofile_all_bare.csv',index=False)
    ### then run this file through the above code
    l=['benign','unknown', None]
    df[(df.prediction.str.contains('damaging'))&(df.pp_hg38.isin(l))&(df.Consequence_hg38=='missense_variant')].to_csv('twofile_more_severe_missense.csv',index=False)
    df[(df.Consequence_level!=0)|(df.Consequence_level_hg38!=0)].to_csv('twofile_analysis_nonzero.csv',index=False)
    # combi=combi.rename({'Feature_custom':'Feature','Amino_acids_custom':'Amino_acids','Protein_position_custom':'Protein_position'},axis=1).merge(polyphen_read('/mnt/xomics/renees/data/host_pathogen_PID/unsolved_cases/comparison/hg38_vs_custom/pph_output_all.features'),on=['Feature','Amino_acids','Protein_position'],how='left')
    # combi['prediction']=combi.prediction.str.strip()
    # combi['change']=combi.apply(lambda x: consequence_change(x['Consequence_level_hg38'],x['Consequence_level_custom']),axis=1)
    # combi=combi.merge(pd.read_csv('/mnt/xomics/renees/data/host_pathogen_PID/unsolved_cases/comparison/hg38_vs_custom/onefile_hg38vcustom.csv')[['Uploaded_variation','Location','Allele','max_af', 'conservation', 'patients','transcript', 'gene_ID', 'transcript_ID', 'annot_gene_id','annot_transcript_id', 'annot_gene_name', 'annot_transcript_name','n_exons', 'length', 'gene_novelty', 'transcript_novelty','ISM_subtype', 'calb', 'rpmi', 'saur', 'calb_total', 'rpmi_total','saur_total', 'annotation', 'calb_exon', 'rpmi_exon', 'saur_exon','prediction', 'based_on', 'effect', 'site', 'region','Ensembl_gene_identifier', 'GO_ID', 'GO_term']],on=['Uploaded_variation','Location','Allele'])
    # combi[(combi.Consequence_level_hg38!=0)&(combi.Consequence_level_custom!=0)].to_csv('/mnt/xomics/renees/data/host_pathogen_PID/unsolved_cases/comparison/hg38_vs_custom/twofile_hg38vcustom.csv',index=False)
    # test=combi[(combi.Consequence_hg38=='missense_variant')&(combi.Consequence_custom=='missense_variant')][['pp_hg38','prediction']].dropna(subset=['prediction'])
    # combi[(combi.Consequence_hg38=='missense_variant')&(combi.Consequence_custom=='missense_variant')&(combi.Feature_custom.str.contains('ENST'))]
    # test['change']=test.pp_hg38.astype('str')+'->'+test.prediction

if __name__ == '__main__':
    main()