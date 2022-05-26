process create_subs {
  /*
  Generate amino acid substitutions
  */

  tag "$vcf"
  // container "biocontainers/bcftools:v1.9-1-deb_cv1"
  memory '4 GB'
  // errorStrategy 'ignore'

  input:
    path vcf

  output:
    path '*.subs'

  """
  export BCFTOOLS_PLUGINS=${params.bcftools_plugins}
  bcftools +split-vep $vcf -d -A tab -s :missense \
           -f '%CHROM-%POS %Feature %Protein_position %Amino_acids\n' | \
           sed 's|/|\\t|g' > protein.subs
  awk '{print \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 > ( \$1 "_" \$2 ".subs" ) }' protein.subs  
  rm protein.subs
  """
}
