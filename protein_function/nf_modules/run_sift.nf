#!/usr/bin/env nextflow

/* 
 * Predict protein function using SIFT
 */

blastdb_name = file(params.blastdb).name

process alignProteins {
  /*
  Run multiple alignment for protein sequence
  */

  tag "$fasta"
  container "nunoagostinho/sift:6.2.1"
  memory '4 GB'
  errorStrategy 'ignore'

  input:
    path fasta
    path blastdb_dir

  output:
    path '*.alignedfasta'

  """
  #!/bin/csh
  setenv tmpdir "."
  setenv NCBI "/opt/blast/bin/"
  seqs_chosen_via_median_info.csh $fasta \
                                  $blastdb_dir/$blastdb_name \
                                  $params.median_cutoff
  """
}

process sift {
  /*
  Run SIFT

  Returns
  -------
  Returns 1 file:
      1) Output 'protein.SIFTprediction'
  */

  tag "${aln}"
  container "nunoagostinho/sift:6.2.1"
  memory '4 GB'
  errorStrategy 'ignore'
  publishDir "${params.outdir}/sift"

  input:
    path aln
    path subs

  output:
    path '*.SIFTprediction'

  """
  info_on_seqs $aln $subs protein.SIFTprediction
  """
}
