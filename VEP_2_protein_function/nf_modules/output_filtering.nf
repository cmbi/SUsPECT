process filter_high_severity {
  container 'ensemblorg/ensembl-vep:latest'

  input:
    file annotated_vcf

  output:
    tuple path('benign_to_pathogenic_dup.vcf'),path('benign_ref_dup.vcf')
    
  """
  filter_vep -i $annotated_vcf --only_matched --filter "Feature matches _ORF_"  | filter_vep -o benign_to_pathogenic_dup.vcf --only_matched --filter "PPH2_score matches damaging or IMPACT is HIGH"
  filter_vep -i $annotated_vcf -o benign_ref_dup.vcf --only_matched --filter "Feature not matches _ORF_"
  """
}

process dedup_output {
  container 'biocontainers/bcftools:v1.5_cv3'

  input:
    tuple path(pathogenic),path(benign)

  output:
    path 'benign_to_pathogenic.vcf'
    path 'benign_ref.vcf'

  """
  bcftools norm -d all -o benign_to_pathogenic.vcf $pathogenic
  bcftools norm -d all -o benign_ref.vcf $benign
  """
}

process map_individuals {
  container 'rlsalz/biopj:0.1.1'

  input:
  path final_vcf
  
  output:
  path "variants_to_patients.tsv"
  
  script:
  """
  python $workflow.projectDir/bin/vcf_to_ind.py $final_vcf
  """
}

process add_annotations_to_indiv {
  container 'griffithlab/vatools:5.0.1'

  input:
    file annotated_vcf
    file mapped_indiv

  output:
    path 'candidates.tsv'

  """
  vep-annotation-reporter -t $mapped_indiv -o candidates.tsv $annotated_vcf Allele Consequence IMPACT Feature CDS_position Protein_position Amino_acids PPH2_score
  """
}

process get_ref_annotations {
  container 'griffithlab/vatools:5.0.1'

  input:
    file annotated_vcf

  output:
    path 'candidates_ref.tsv'

  """
  vep-annotation-reporter -o candidates_ref.tsv $annotated_vcf Consequence SYMBOL Amino_acids Gene BIOTYPE MAX_AF
  """
}

process combine_custom_ref_candidates {
  container 'rlsalz/biopj:0.1.1'
  publishDir "${params.outdir}/"

  input:
    file candidates
    file candidates_ref

  output:
    path 'candidates_info.tsv'

  """
  #!/usr/bin/env python3
  import pandas as pd

  def remove_dups(commasep):
    return(','.join(list(set(str(commasep).split(',')))))

  def checkaa(sub,sublist,cons):
    return(sub in sublist and cons=='missense_variant')
  
  df=pd.read_table($candidates).merge(pd.read_table($candidates_ref),on=['CHROM','POS','REF','ALT'],suffixes=('','_ref')).fillna('-').applymap(remove_dups).applymap(lambda x: str(x).strip(','))
  df['same_aa']=df.apply(lambda x: checkaa(x['Amino_acids'],x['Amino_acids_ref'].split(','),x['Consequence']),axis=1)
  df[~df['same_aa']].drop(columns=['same_aa']).to_csv('candidates_info.tsv',sep='\t',index=False)
  
  
  """
}