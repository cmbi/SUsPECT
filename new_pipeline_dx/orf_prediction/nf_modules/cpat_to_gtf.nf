
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
    path "novel_with_cds.gtf"
  
  script:
  """
  import pandas as pd
  from collections import defaultdict
  import copy
  import gtfparse

  def get_first_block_index(orf_coord, cblens):
      # get the index corresponding to the first block containing the orf start
      # return index, and the dela (spacing upstream of end)
      for i, cblen in enumerate(cblens):
          if orf_coord <= cblen:
              delta = cblen - orf_coord
              return i, delta

  def make_coords_trimmed_to_orf_range(i1, delta1, i2, delta2, coords):
      orf_coords = copy.deepcopy(coords)
      orf_coords = orf_coords[i1: i2+1]
      # trim ends to orf start/end
      orf_coords[0][0] = orf_coords[0][1] - delta1
      orf_coords[-1][1] = orf_coords[-1][1] - delta2
      return(orf_coords)

  def make_coords_trimmed_to_orf_range_neg_strand(i1, delta1, i2, delta2, coords):
      orf_coords = copy.deepcopy(coords)
      orf_coords = orf_coords[i1: i2+1] 
      # trim ends to orf start/end
      orf_coords[0][1] = orf_coords[0][0] + delta1
      orf_coords[-1][0] = orf_coords[-1][0] + delta2
      return(orf_coords)

  def make_pacbio_cds_gtf(sample_gtf, cpat_orfs):
      # import gtf, only exon info.
      # only move forward with representative pb isoform (for same-protein groups)
      gtf = gtfparse.read_gtf(sample_gtf)
      gtf['transcript_id']=gtf['transcript_id'].str.upper()
      talon_gene = pd.Series(gtf.gene_id.values, index=gtf.transcript_id).to_dict()


      gtf = gtf[['seqname', 'feature', 'start', 'end', 'strand', 'transcript_id']]
      gtf = gtf[gtf['feature'] == 'exon']
      gtf.columns = ['chr', 'feat', 'start', 'end', 'strand', 'acc']


      # pb coords into dict
      pbs = defaultdict(lambda: ['chr', 'strand', [], [], []]) # pb -> [chr, strand, [start, end], [block lengths], [cum. block lengths]]
      # PB.1.1 -> ['chr1', '+', [[100,150], [200,270]], [50, 70], [50, 120], [150-200]]
      for i, row in gtf.iterrows():
          chr, feat, start, end, strand, acc = row
          pbs[acc][0] = chr
          pbs[acc][1] = strand
          pbs[acc][2].append([int(start), int(end)])
      # sort all coords, calc blocks
      def make_cumulative_blens(blocks):
          cblocks = []
          cbl = 0 # cumulative block length
          for b in blocks:
              cbl += b
              cblocks.append(cbl)
          return cblocks
      for acc, infos in pbs.items():
          strand = infos[1]
          if strand == '+':
              infos[2] = sorted(infos[2])
          elif strand == '-':
              infos[2] = sorted(infos[2], reverse=True)
          infos[3] = [end-start+1 for [start, end] in infos[2]]
          infos[4] = make_cumulative_blens(infos[3])


      # read in the ranges of orf on pb transcripts
      ranges = pd.read_table(cpat_orfs)
      ranges= ranges[ranges['Coding_prob']>0.364][['ID','ORF_start','ORF_end']] # human cutoff
      ranges['acc']= ranges['ID'].apply(lambda x: x.split('_')[0])


      with open('novel_with_cds.gtf', 'w') as ofile:
          for i, row in ranges.iterrows():
              oid, orf_start, orf_end, acc = row
              if acc in pbs:
                  gene = talon_gene[acc]
                  infos = pbs[acc]
                  chr, strand, coords, blens, cblens = infos
                  # NOTE - uncomment to do chr22-oly test 
                  # if chr != 'chr22': continue
                  i1, delta1 = get_first_block_index(orf_start, cblens)
                  i2, delta2 = get_first_block_index(orf_end, cblens)
                  if strand == '+':
                      orf_coords = make_coords_trimmed_to_orf_range(i1, delta1, i2, delta2, coords)
                  elif strand == '-':
                      orf_coords = make_coords_trimmed_to_orf_range_neg_strand(i1, delta1,i2, delta2, coords)
                  # write out the coordinates
                  acc_w_gene_w_cpm = gene + '|' + oid
                  out_acc = 'gene_id "{}"; transcript_id "{}";'.format(gene, acc_w_gene_w_cpm)
                  for [start, end] in coords:
                      ofile.write('\t'.join([chr, 'hg38_canon', 'exon', str(start), str(end), '.', strand, '.', out_acc]) + '\n')
                  for [start, end] in orf_coords:
                      ofile.write('\t'.join([chr, 'hg38_canon', 'CDS', str(start), str(end), '.', strand, '.', out_acc]) + '\n')


  make_pacbio_cds_gtf('$sample_gtf', '$cpat_orfs')
  """
}

/*--------------------------------------------------
Combine/organize the final gtf of novel sequences
 * Create a final GTF file to be given to VEP
---------------------------------------------------*/
process create_final_gtf{
    publishDir "${params.outdir}/${params.name}/final_gtf/", mode: 'copy'
    tag "${params.name} ${reference_gtf} ${sample_gtf}"
    cpus 1
    container "quay.io/biocontainers/agat:0.9.0--pl5321hdfd78af_0"

    input:
        path sample_gtf
        path sample_cds
    output:
        path "${params.name}_complete.gtf"
        

    script:
        """
        agat_sp_complement_annotations.pl --ref $sample_gtf --add $sample_cds --out ${params.name}_complete.gtf 
        """
}