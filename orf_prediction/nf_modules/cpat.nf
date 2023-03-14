process convert_to_bed {
  cpus 1
  label 'agat'
  storeDir "${params.outdir}/${params.name}/transcriptome_fasta/"
  memory '4 GB'

  input:
  path sample_gtf_cleaned
  
  output:
  path "novel_transcripts.bed"
  
  script:
  """
  agat_convert_sp_gff2bed.pl --gff $sample_gtf_cleaned -o novel_transcripts.bed
  """
}

/*--------------------------------------------------
CPAT
 * CPAT is a bioinformatics tool to predict an RNA’s coding probability 
 * based on the RNA sequence characteristics. 
 * To achieve this goal, CPAT calculates scores of sequence-based features 
 * from a set of known protein-coding genes and background set of non-coding genes.
 *     ORF size
 *     ORF coverage
 *     Fickett score
 *     Hexamer usage bias
 * 
 * CPAT will then builds a logistic regression model using these 4 features as 
 * predictor variables and the “protein-coding status” as the response variable. 
 * After evaluating the performance and determining the probability cutoff, 
 * the model can be used to predict new RNA sequences.
 *
 * https://cpat.readthedocs.io/en/latest/
---------------------------------------------------*/
process cpat {
  cpus 1
  label 'cpat'
  
  storeDir "${params.outdir}/${params.name}/cpat/"

  input:
  path hexamer
  path logit_model
  path sample_bed
  path reference_fasta

  output:
  path "novel_seqs.ORF_prob.tsv"
  path "novel_seqs.ORF_prob.best.tsv"
  path "novel_seqs.ORF_seqs.fa"

  script:
  """
  cpat.py -x $hexamer -d $logit_model -g $sample_bed -r $reference_fasta -o novel_seqs 
  """
}

process cpat_orf_to_protein {
  cpus 1
  label 'emboss'
  
  storeDir "${params.outdir}/${params.name}/cpat/"

  input:
  path orf_seqs

  output:
  path 'novel_seqs.prot_seqs.fa'

  script:
  """
  transeq $orf_seqs novel_seqs.prot_seqs.fa -trim
  sed -i 's/_[^_]*//3g' novel_seqs.prot_seqs.fa
  """

}
