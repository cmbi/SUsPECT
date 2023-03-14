#!/usr/bin/env nextflow

/* 
 * Predict protein function using SIFT
 */

process align_peptide {
  /*
  Run multiple alignment for protein sequence
  */

  tag "${peptide_id}"
  label 'sift'
  errorStrategy 'ignore'

  input:
    tuple val(peptide_id), path(fasta), path(subs)
    path blastdb_dir
    val blastdb_name

  output:
    tuple val(peptide_id), path('*.alignedfasta'), path(subs), optional: true

  """
  #!/bin/csh
  setenv tmpdir "."
  setenv NCBI "/opt/blast/bin/"
  seqs_chosen_via_median_info.csh ${fasta} \
                                  ${blastdb_dir}/${blastdb_name} \
                                  ${params.median_cutoff}
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

  tag "${peptide_id}"
  label 'sift'
  errorStrategy 'ignore'

  input:
    tuple val(peptide_id), path(aln), path(subs)

  output:
    path 'sift.out'

  """
  awk '{print \$7\$6\$8}' ${subs} > var.subs
  info_on_seqs ${aln} var.subs protein.SIFTprediction

  # Clean results and append variant location and transcript identifier
  echo "#Chrom,Pos,Ref,Alt,Transcript,Prediction,Score" |\
    sed 's/,/\t/g' > sift.out
  paste <(cut -f 1-5 ${subs}) <(cut -f 2-3 protein.SIFTprediction) >> sift.out
  """
}
