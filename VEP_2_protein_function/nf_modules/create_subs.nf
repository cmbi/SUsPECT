process create_subs {
  /*
  Generate amino acid substitutions
  */

  tag "$vcf"
  container "quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0"
  memory '4 GB'
  // errorStrategy 'ignore'

  input:
    path vcf

  output:
    path '*.subs'

  """
  bcftools +split-vep $vcf -d -A tab -s :missense \
           -f '%CHROM-%POS-%REF-%ALT %Feature %Protein_position %Amino_acids\n' | \
           sed 's|/|\\t|g' > protein.subs
  awk '{ file=\$1 "-" \$2 ".subs"; print \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 > file; close(file) }' protein.subs
  rm protein.subs
  """
}
