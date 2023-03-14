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
  label 'agat'
  publishDir "${params.outdir}"
  memory = { ["4 GB", "8 GB", "20 GB", "40 GB"][task.attempt - 1] }

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
