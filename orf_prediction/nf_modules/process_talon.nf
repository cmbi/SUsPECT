process fetch_novel {
  publishDir "${params.outdir}/${params.name}/novel_gtf/", mode: 'copy'
  cpus 1
  
  input:
      path novel_gtf_raw
      val idprefix
  output:
      path "talon_novelonly.gtf"
      

  script:
      """
      grep $idprefix $novel_gtf_raw > talon_novelonly.gtf
      """
}

process uppercase {
  publishDir "${params.outdir}/${params.name}/novel_gtf/", mode: 'copy'
  cpus 1
  container 'quay.io/biocontainers/gtfparse:1.2.1--pyh864c0ab_0'

  input:
      path novel_gtf
    
  output:
      path "talon_novelonly_UP.gtf"
    
  script:
      """
      #!/usr/bin/env python3

      import pandas as pd
      import csv
      def upper_tid(extra):
        marker='transcript_id "'
        if marker in extra:
            begin=extra.split(marker,1)
            tid=begin[1]
            begin=begin[0]
            end=tid.split('"',1)
            return(begin+marker+end[0].upper()+'"'+end[1])
        else:
            return(extra)
      df=pd.read_table('$novel_gtf', names=['seqname','source','feature','start','end','score','strand','frame','extra'])
      df['extra']=df['extra'].apply(upper_tid) #make everything upper case to prepare for CPAT
      df.to_csv('talon_novelonly_UP.gtf',sep='\t',index=False,header=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar="")
      """
}