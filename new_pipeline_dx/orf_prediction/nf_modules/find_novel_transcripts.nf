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
      path sample_gtf
  output:
      path "gffcmp.combined.filtered.gtf"
      

  script:
      """
      #!/usr/bin/env python3

      import pandas as pd
      import gtfparse
      import csv
      
      #determine the novel seqs and save them
      combined=gtfparse.read_gtf('$novel_marked_gtf')
      exclude=combined[combined['class_code']=='=']['transcript_id'].unique()
      combined=combined[~combined['transcript_id'].isin(exclude)]
      combined=combined[combined['oId'].astype(bool)]
      acceptable_genes=combined.gene_name.unique()
      acceptable_transcripts=combined.oId.unique()
      #filter the gtf
      old=gtfparse.read_gtf('$sample_gtf')
      old=old[((old.gene_name.isin(acceptable_genes))&(old.feature=='gene'))|(old.transcript_id.isin(acceptable_transcripts))]
      #print the filtered gtf
      df=pd.read_table('$sample_gtf', names=['seqname','source','feature','start','end','score','strand','frame','extra'])
      df.merge(old[['seqname','feature','start','end']]).drop_duplicates().to_csv('gffcmp.combined.filtered.gtf',sep='\t',index=False,header=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar="")
      """
}