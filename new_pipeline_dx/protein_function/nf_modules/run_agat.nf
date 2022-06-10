#!/usr/bin/env nextflow

/*
 * AGAT: Another GTF/GFF Analysis Toolkit
 */

process getTranslation {
  /*
  Translate nucleotide FASTA sequences based on GTF features

  Returns
  -------
  Returns 1 file:
      1) Protein FASTA sequence 'translated.fa'
  */

  tag "${gtf}"
  container "quay.io/biocontainers/agat:0.9.0--pl5321hdfd78af_0"
  memory '20 GB'
  publishDir "${params.outdir}"

  input:
    path gtf
    path fasta

  output:
    path '*_translated.fa'

  script:
    """
    agat_sp_extract_sequences.pl -g $gtf -f $fasta --protein \
                                 -o ${gtf.baseName}_translated.fa
    """
}
