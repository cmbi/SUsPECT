process filter_high_severity {
  container 'ensemblorg/ensembl-vep:latest'

  input:
    file annotated_vcf

  output:
    path 'benign_to_pathogenic.vcf'

  """
  filter_vep -i $annotated_vcf --only_matched --filter "Feature matches _ORF_"  | filter_vep -o benign_to_pathogenic_dup.vcf --only_matched --filter "PPH2_score matches damaging or IMPACT is HIGH"
  bcftools norm -d all -o benign_to_pathogenic.vcf benign_to_pathogenic_dup.vcf
  """
}

process filter_ref_only {
  container 'ensemblorg/ensembl-vep:latest'

  input:
    file annotated_vcf

  output:
    path 'benign_ref.vcf'

  """
  filter_vep -i $annotated_vcf -o benign_ref_dup.vcf --only_matched --filter "Feature matches _ORF_"
  bcftools norm -d all -o benign_ref.vcf benign_ref_dup.vcf
  """
}

process map_individuals {
  container 'rlsalz/biopj:0.1.1'

  publishDir "${params.outdir}/${params.name}/cds/", mode: 'copy'

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
  vep-annotation-reporter -t $mapped_indiv -o candidates.tsv $annotated_vcf Allele Consequence IMPACT SYMBOL CDS_position Protein_position Amino_acids PPH2_score
  """
}

process get_ref_annotations {
  container 'griffithlab/vatools:5.0.1'

  input:
    file annotated_vcf

  output:
    path 'candidates_ref.tsv'

  """
  vep-annotation-reporter -o candidates_ref.tsv $annotated_vcf Consequence SYMBOL Gene BIOTYPE MAX_AF
  """
}

process combine_custom_ref_candidates {
  container 'rlsalz/biopj:0.1.1'

  input:
    file candidates
    file candidates_ref

  output:
    path 'candidates_info.tsv'

  """
  #!/usr/bin/env python3
  import pandas as pd
  pd.read_table($candidates).merge(pd.read_table($candidates_ref),on=['CHROM','POS','REF','ALT'],suffixes=('','_ref')).to_csv('candidates_info.tsv',sep='\t')
  
  """
}