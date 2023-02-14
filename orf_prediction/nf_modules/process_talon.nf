process fetch_novel_talon {
  publishDir "${params.outdir}/${params.name}/novel_gtf/", mode: 'copy'
  cpus 1
  
  input:
      path novel_gtf_raw
  output:
      path "talon_novelonly.gtf"
      

  script:
      """
      grep -e '\\tgene\\t' -e 'transcript_status\s\"NOVEL\"' $novel_gtf_raw | sed '/\\tgene\\t/{\$!N;/\\n.*\\tgene\\t/!P;D}' > talon_novelonly.gtf
      """
}

process fetch_novel_isoquant {
  publishDir "${params.outdir}/${params.name}/novel_gtf/", mode: 'copy'
  cpus 1
  
  input:
      path novel_gtf_raw
  output:
      path "isoquant_novelonly.gtf"
      

  script:
      """
      grep -e '\\sgene\\s' -e '\\.nic' -e '\\.nnic' $novel_gtf_raw | awk '! ( \$1 ~ /_/ )' | awk '! ( \$7=="." )' | sed '/\\sgene\\s/{\$!N;/\\n.*\\sgene\\s/!P;D}' > isoquant_novelonly.gtf
      """
}
