process linearise_fasta {
  /*
  Linearise FASTA file
  */

  tag "${fasta.baseName}"

  input:
    path fasta

  output:
    path 'linear.fa'

  script:
  """
  sed -e 's/\\(^>.*\$\\)/#\\1#/' $fasta | \
    tr -d "\\r" | \
    tr -d "\\n" | \
    sed -e 's/\$/#/' | \
    tr "#" "\\n" | sed -e '/^\$/d' > linear.fa
  """
}

process filter_peptide {
  /*
  Filters linearlised FASTA and substitions file based on peptide identifier
  */

  tag "${fasta.baseName}"

  input:
    val peptide_id
    path fasta
    path subs

  output:
    tuple val(peptide_id), path('filtered.fa'), path('filtered.subs')

  """
  grep ${peptide_id} $subs > filtered.subs
  grep ${peptide_id} -A1 $fasta > filtered.fa
  """
}
