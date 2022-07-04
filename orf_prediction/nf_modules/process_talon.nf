process fetch_novel {
  publishDir "${params.outdir}/${params.name}/novel_gtf/", mode: 'copy'
  cpus 1
  
  input:
      path novel_gtf_raw
  output:
      path "talon_novelonly.gtf"
      

  script:
      """
      grep -e '\tgene\t' -e 'transcript_status\s\"NOVEL\"' $novel_gtf_raw | sed '/\tgene\t/{\$!N;/\\n.*\tgene\t/!P;D}' > talon_novelonly.gtf
      """
}
