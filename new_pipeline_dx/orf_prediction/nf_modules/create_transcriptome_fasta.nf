/*--------------------------------------------------
Create transcript fasta
 * Takes novel GTF input and creates a fasta file of cDNA
---------------------------------------------------*/
process create_transcriptome_fasta {
  cpus 1
  conda 'bioconda::transdecoder'
  publishDir "${params.outdir}/${params.name}/transcriptome_fasta/", mode: 'copy'

  input:
  path genome_fasta
  path sample_gtf
  path transdecoder_dir
  
  output:
  path "novel_transcripts.fasta"
  
  script:
  """
  $transdecoder_dir/util/gtf_genome_to_cdna_fasta.pl $sample_gtf $genome_fasta > novel_transcripts.fasta
  """
}