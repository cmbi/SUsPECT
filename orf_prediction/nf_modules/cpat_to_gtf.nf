
/*--------------------------------------------------
CDS GTF 
 * create GTF that also contains CDS attributes
---------------------------------------------------*/
process cpat_to_bed {
  cpus 1
  container 'rlsalz/biopj:0.1.1'

  storeDir "${params.outdir}/${params.name}/cds/"

  input:
  path transcript_bed
  path cpat_orfs
  
  output:
  path "orfs.bed"
  
  script:
  """
  python $workflow.projectDir/bin/cpat_to_bed.py $cpat_orfs $transcript_bed orfs.bed
  """
}

process combine_bed {
  storeDir "${params.outdir}/${params.name}/final_gtf/"
  cpus 1
  container "quay.io/biocontainers/gtfparse:1.2.1--pyh864c0ab_0"

  input:
      path transcript_bed
      path orf_bed
  output:
      path "transcripts_and_orfs.bed"
      

  script:
      """
      #!/usr/bin/env python3

      import pandas as pd

      transcripts=pd.read_table("${transcript_bed}",names = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])
      orfs=pd.read_table("${orf_bed}", header=0, names = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])
      orfs['transcript']=orfs['name'].str.split('_',1).apply(lambda x: x[0])
      combined=transcripts.merge(orfs[['transcript','name','thickStart','thickEnd']],left_on='name',right_on='transcript',suffixes=('','_orf'))
      combined[['chrom', 'chromStart', 'chromEnd', 'name_orf', 'score', 'strand', 'thickStart_orf', 'thickEnd_orf', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']].to_csv("transcripts_and_orfs.bed",sep='\t',header=False,index=False)
      """
}

process bed_to_genepred {
  storeDir "${params.outdir}/${params.name}/final_gtf/"
  cpus 1
  container "quay.io/biocontainers/ucsc-bedtogenepred:377--ha8a8165_3"

  input:
      path whole_bed
  output:
      path "transcripts_and_orfs.genepred"
      

  script:
      """
      bedToGenePred $whole_bed transcripts_and_orfs.genepred
      """
}

process genepred_to_gtf {
  storeDir "${params.outdir}/${params.name}/final_gtf/"
  cpus 1
  container "quay.io/biocontainers/ucsc-genepredtogtf:377--ha8a8165_5"

  input:
      path sample_genepred
  output:
      path "transcripts_and_orfs.gtf"
      

  script:
      """
      genePredToGtf -utr -honorCdsStat file $sample_genepred transcripts_and_orfs.gtf
      """
}

process add_genes {
  storeDir "${params.outdir}/${params.name}/final_gtf/"
  cpus 1

  input:
      path final_gtf
  output:
      path "transcripts_and_orfs_with_gene.gtf"
      

  script:
      """
      sed '/\ttranscript\t/p' $final_gtf | sed -e 's/\ttranscript\t/\tgene\t/g;n' | sed '/\tgene\t/{s/transcript_id[^.]*\$//}' | sed '/gene_id\s".*"/{s/_ORF//}' > transcripts_and_orfs_with_gene.gtf
      """
}

// process merge_cds_with_rest {
//   publishDir "${params.outdir}/${params.name}/final_gtf/", mode: 'copy'
//   cpus 1
//   container 'quay.io/biocontainers/gtfparse:1.2.1--pyh864c0ab_0'
//   memory '16 GB'

//   input:
//       path novel_cds_gff
//       path whole_gtf
//   output:
//       path "novel_cds_formatted.gtf"
      

//   script:
//       """
//       #!/usr/bin/env python3

//       import pandas as pd
//       import gtfparse
//       import csv

//       def get_parent(extra):
//         if 'Parent=' in extra:
//           return(extra.split('Parent=')[1])
//         else:
//           p=extra.split('ID=')[1]
//           return(p.split(';',1)[0])
      
//       def make_last_column(type_feature,gene_id,gene_name,transcript_id,exon_id):
//         if type_feature=='gene':
//           return('ID='+gene_id+';Name='+gene_name)
//         elif type_feature=='transcript':
//           return('ID='+transcript_id+';Parent='+gene_id)
//         elif type_feature=='exon':
//           return('ID='+exon_id+';Parent='+transcript_id)
      
//       #read in CPAT orf bed file
//       gff_orfs=pd.read_table("${novel_cds_gff}",names=['seqname', 'source', 'feature', 'start','end', 'score', 'strand','frame','extra'])
//       gff_orfs['parent']=gff_orfs.extra.apply(get_parent)
//       gff_orfs['ccdsid']=gff_orfs.extra.apply(lambda x: x.split(';Parent=')[0] if ';Parent=' in x else 'ID=')
//       gff_orfs['ccdsid']=gff_orfs['ccdsid'].str.strip('ID=')
//       #create mapping file
//       mapper=gff_orfs[gff_orfs.extra.str.contains(';Name=')][['extra']]
//       mapper['extra']=mapper.extra.str.split(';blockCount').apply(lambda x: x[0])
//       mapper['extra']=mapper.extra.str.split(';Name=')
//       mapper[['parent','ORF_id']]=pd.DataFrame(mapper.extra.tolist(),index=mapper.index)
//       mapper['parent']=mapper.parent.str.strip('ID=')
//       mapper['transcript_id']= mapper.ORF_id.str.split('_',1).apply(lambda x: x[0])
//       #read in full gtf
//       gtf_whole=gtfparse.read_gtf("${whole_gtf}").fillna('.')
//       gtf_whole['transcript_id']=gtf_whole['transcript_id'].str.upper()
//       mapper=mapper.drop(columns='extra').merge(gtf_whole[['transcript_id','gene_id']].drop_duplicates())
//       acceptable_genes=mapper.gene_id.unique()
//       acceptable_transcripts=mapper.transcript_id.unique()
//       gtf_whole=gtf_whole[((gtf_whole.gene_id.isin(acceptable_genes))&(gtf_whole.feature=='gene'))|(gtf_whole.transcript_id.isin(acceptable_transcripts))]
//       gtf_whole_dup=gtf_whole.merge(mapper.drop(columns='gene_id'),on='transcript_id',how='left')
//       gtf_whole_dup['extra']=gtf_whole_dup.apply(lambda x: make_last_column(x['feature'],x['gene_id'],x['gene_name'],x['ORF_id'],x['exon_id']),axis=1)
//       gtf_whole_dup['feature']=gtf_whole_dup['feature'].str.replace('transcript','mRNA',regex=False)
//       out_orfs=gff_orfs[gff_orfs.feature=='CDS'].merge(mapper,on='parent')
//       out_orfs['extra']='ID='+out_orfs['ccdsid']+';Parent='+out_orfs['ORF_id']
//       #print out final gff3 file
//       pd.concat([gtf_whole_dup[['seqname', 'source', 'feature', 'start','end', 'score', 'strand','frame','extra']],out_orfs[['seqname', 'source', 'feature', 'start','end', 'score', 'strand','frame','extra']]]).fillna('.').to_csv('novel_cds_formatted.gtf',sep='\t',index=False,header=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar="")
//       """
// }

// /*--------------------------------------------------
// Combine/organize the final gtf of novel sequences
//  * Create a final GFF file to be given to VEP
// ---------------------------------------------------*/
// process create_final_gff{
//   publishDir "${params.outdir}/${params.name}/final_gtf/", mode: 'copy'
//   cpus 1
//   container "quay.io/biocontainers/agat:0.9.0--pl5321hdfd78af_0"
//   memory '8 GB'

//   input:
//       path unsorted_gff
//   output:
//       path "complete_sorted.gff3"
      

//   script:
//       """
//       agat_convert_sp_gxf2gxf.pl -g $unsorted_gff --tabix -o complete_sorted.gff3 
//       """
// }
