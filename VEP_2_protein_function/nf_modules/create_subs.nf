process create_subs {
  /*
  Generate amino acid substitutions
  */

  tag "$vcf"
  container "quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0"

  input:
    path vcf

  output:
    path 'protein.subs'

  """
  echo "chrom\tpos\tref\talt\tfeature\tprotein_pos\taa_ref\taa_alt" > protein.subs
  bcftools +split-vep $vcf -d -A tab -s :missense \
           -f '%CHROM\t%POS\t%REF\t%ALT\t%Feature\t%Protein_position\t%Amino_acids\n' |\
           sed 's|/|\\t|g' >> protein.subs
  """
}
