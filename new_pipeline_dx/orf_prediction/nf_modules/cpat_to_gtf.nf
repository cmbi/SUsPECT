
/*--------------------------------------------------
CDS GTF 
 * create GTF that also contains CDS attributes
---------------------------------------------------*/
process make_cds_gtf {
  cpus 1
  container 'quay.io/biocontainers/gtfparse:1.2.1--pyh864c0ab_0'

  publishDir "${params.outdir}/${params.name}/cds/", mode: 'copy'

  input:
  path sample_gtf
  path cpat_orfs
  
  output:
  path "novel_cds.gtf"
  
  script:
  """
  python $workflow.projectDir/bin/cpat_to_gtf.py --sample_gtf $sample_gtf --cpat_orfs $cpat_orfs
  """
}

/*--------------------------------------------------
Combine/organize the final gtf of novel sequences
 * Create a final GTF file to be given to VEP
---------------------------------------------------*/
process create_final_gtf{
  publishDir "${params.outdir}/${params.name}/final_gtf/", mode: 'copy'
  cpus 1
  container "quay.io/biocontainers/agat:0.9.0--pl5321hdfd78af_0"

  input:
      path sample_gtf
      path sample_cds
  output:
      path "novel_complete.gtf"
      

  script:
      """
      agat_sp_complement_annotations.pl --ref $sample_gtf --add $sample_cds --out novel_complete.gtf 
      """
}