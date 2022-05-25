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
  tag "${hexamer}, ${logit_model}, ${sample_fasta}"

  publishDir "${params.outdir}/${params.name}/cpat/", mode: 'copy'

  input:
  path hexamer
  path logit_model
  path sample_fasta

  output:
  path "${params.name}.ORF_prob.tsv"
  path "${params.name}.ORF_prob.best.tsv"
  path "${params.name}.ORF_seqs.fa"
  path "*"

  script:
  """
  cpat.py -x $hexamer -d $logit_model -g $sample_fasta -o ${params.name} 
  """
}