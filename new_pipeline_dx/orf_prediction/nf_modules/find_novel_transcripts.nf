/*--------------------------------------------------
Find novel transcripts
For those that want to submit a full GTF file
 * Takes GTF input and creates a fasta file of cDNA
---------------------------------------------------*/
process identify_novel {
  publishDir "${params.outdir}/${params.name}/novel_gtf/", mode: 'copy'
  tag "${params.name} ${reference_gtf} ${sample_gtf}"
  cpus 1
  conda 'bioconda::gffcompare'

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
  container ""

  input:
      path novel_marked_gtf
  output:
      path "gffcmp.combined.filtered.gtf"
      

  script:
      """
      python -r $workflow.projectDir/cripts/filter_novel.py --gffcompare_gtf $novel_marked_gtf 
      """
}