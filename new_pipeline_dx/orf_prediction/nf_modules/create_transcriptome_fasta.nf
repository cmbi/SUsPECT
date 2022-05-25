/*--------------------------------------------------
Create transcript fasta
 * Takes novel GTF input and creates a fasta file of cDNA
---------------------------------------------------*/
process create_transcriptome_fasta {
  cpus 1
  container 'quay.io/biocontainers/agat:0.9.0--pl5321hdfd78af_0'
  publishDir "${params.outdir}/${params.name}/transcriptome_fasta/", mode: 'copy'

  input:
  path genome_fasta
  path sample_gtf
  
  output:
  path "novel_transcripts.fasta"
  
  script:
  """
  agat_sp_extract_sequences.pl -g $sample_gtf -f $genome_fasta -o novel_transcripts.fasta -t exon --merge
  """
}