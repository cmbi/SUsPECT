process convert_to_bed {
  cpus 1
  container 'quay.io/biocontainers/agat:0.9.0--pl5321hdfd78af_0'
  publishDir "${params.outdir}/${params.name}/transcriptome_fasta/", mode: 'copy'
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
  container "quay.io/biocontainers/cpat:3.0.4--py38h17adfb0_1"
  
  publishDir "${params.outdir}/${params.name}/cpat/", mode: 'copy'

  input:
  path hexamer
  path logit_model
  path sample_bed
  path reference_fasta

  output:
  path "novel_seqs.ORF_prob.tsv"
  path "novel_seqs.ORF_prob.best.tsv"
  path "novel_seqs.ORF_seqs.fa"
  path "*"

  script:
  """
  cpat.py -x $hexamer -d $logit_model -g $sample_bed -r $reference_fasta -o novel_seqs 
  """
}
