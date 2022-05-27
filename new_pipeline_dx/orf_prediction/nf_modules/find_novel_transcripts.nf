/*--------------------------------------------------
Find novel transcripts
For those that want to submit a full GTF file
 * Takes GTF input and creates a fasta file of cDNA
---------------------------------------------------*/
process identify_novel {
  publishDir "${params.outdir}/${params.name}/novel_gtf/", mode: 'copy'
  tag "${params.name} ${reference_gtf} ${sample_gtf}"
  cpus 1
  container 'quay.io/biocontainers/gffcompare:0.11.2--h6bb024c_0'

  input:
      path sample_gtf
      path reference_gtf
  output:
      path "gffcmp.combined.gtf"
      path "*"
      

  script:
      """
      gffcompare -r $reference_gtf -R -A -T $sample_gtf 
      """
}

process filter_novel {
  publishDir "${params.outdir}/${params.name}/novel_gtf/", mode: 'copy'
  tag "${params.name} ${reference_gtf} ${sample_gtf}"
  cpus 1
  container 'quay.io/biocontainers/gtfparse:1.2.1--pyh864c0ab_0'

  input:
      path novel_marked_gtf
  output:
      path "gffcmp.combined.filtered.gtf"
      

  script:
      """
      python -r $workflow.projectDir/scripts/filter_novel.py --gffcompare_gtf $novel_marked_gtf 
      """
}