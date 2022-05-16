#!/usr/bin/env nextflow
/*
 * 
 *   This script was adapted from the "long read proteogenomics" pipeline written by Sheynkman Lab.
 *
 *
 * @authors
 * Renee Salz
 * Nuno Saraiva-Agostinho
 */

def helpMessage() {
    log.info logHeader()
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
      
      #TODO

    Other:

      --max_cpus                    Maximum number of CPUs (int)
      --max_memory                  Maximum memory (memory unit)
      --max_time                    Maximum time (time unit)

    See here for more info: 
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

log.info "ORF prediction novel sequences - N F  ~  version 0.1"
log.info "====================================="
// Header log info
log.info "\nPARAMETERS SUMMARY"
log.info "genome_fasta                          : ${params.genome_fasta}"
log.info "hexamer                               : ${params.hexamer}"
log.info "logit_model                           : ${params.logit_model}"
log.info "transdecoder_dir                      : ${params.transdecoder_dir}"
// log.info "max_cpus                              : ${params.max_cpus}"
log.info "name                                  : ${params.name}"
log.info "vcf                                   : ${params.vcf}"
log.info "outdir                                : ${params.outdir}"
log.info "sample_gtf                            : ${params.sample_gtf}"
log.info ""


if (!params.sample_gtf) exit 1, "Cannot find gtf file for parameter --sample_gtf: ${params.sample_gtf}"
ch_sample_gtf = Channel.value(file(params.sample_gtf))

if (!params.hexamer) exit 1, "Cannot find headmer file for parameter --hexamer: ${params.hexamer}"
ch_hexamer = Channel.value(file(params.hexamer))

if (!params.transdecoder_dir) exit 1, "Cannot find headmer file for parameter --transdecoder_dir: ${params.transdecoder_dir}"
ch_transdecoder_dir = Channel.value(file(params.transdecoder_dir))

if (!params.logit_model) exit 1, "Cannot find any logit model file for parameter --logit_model: ${params.logit_model}"

if (!params.vcf) exit 1, "Cannot find any vcf file for parameter --vcf: ${params.vcf}"
ch_normalized_ribo_kallisto = Channel.value(file(params.vcf))

if (params.vcf.endsWith('.gz')){
   ch_vcf = Channel.value(file(params.vcf))
} else {
   ch_vcf_uncompressed = Channel.value(file(params.vcf))
}

if (params.genome_fasta.endsWith('.gz')){
   ch_genome_fasta = Channel.value(file(params.genome_fasta))
} else {
   ch_genome_fasta_uncompressed = Channel.value(file(params.genome_fasta))
}

if (params.logit_model.endsWith('.gz')) {
   ch_logit_model = Channel.value(file(params.logit_model))
} else {
   ch_logit_model_uncompressed = Channel.value(file(params.logit_model))
}


// Implements logic for cloud compatibility, NO_TOML_FILE as variable only works for envs with local file system
projectDir = workflow.projectDir


/*--------------------------------------------------
Decompress Logit Model
---------------------------------------------------*/
if (params.logit_model.endsWith('.gz')) {
   process gunzip_logit_model {
      tag "decompress logit model"
      cpus 1

      input:
      file(logit_model) from ch_logit_model

      output:
      file("*.RData") into ch_logit_model_uncompressed

      script:
      """
      gunzip -f ${logit_model}
      """
   }
}

/*--------------------------------------------------
Decompress VCF file
---------------------------------------------------*/
if (params.vcf.endsWith('.gz')) {
   process gunzip_vcf {
      tag "decompress VCF file"
      cpus 1

      input:
      file(vcf) from ch_vcf

      output:
      file("*.vcf") into ch_vcf_uncompressed

      script:
      """
      gunzip -f ${vcf}
      """
   }
}

/*--------------------------------------------------
Decompress genome fasta file
---------------------------------------------------*/
if (params.genome_fasta.endsWith('.gz')) {
   process gunzip_gencome_fasta {
   tag "decompress gzipped genome fasta"
   cpus 1

   input:
   file(genome_fasta) from ch_genome_fasta

   output:
   file("*.{fa,fasta}") into ch_genome_fasta_uncompressed

   script:
   """
   gunzip -f ${genome_fasta}
   """
   }
}



/*--------------------------------------------------
Create transcript fasta
 * Takes novel GTF input and creates a fasta file of cDNA
---------------------------------------------------*/
process create_transcriptome_fasta {
  cpus 1
  conda 'bioconda::transdecoder'
  publishDir "${params.outdir}/${params.name}/transcriptome_fasta/", mode: 'copy'

  input:
  file(genome_fasta) from ch_genome_fasta_uncompressed
  file(sample_gtf) from ch_sample_gtf
  file(transdecoder_dir) from ch_transdecoder_dir
  
  output:
  file("novel_transcripts.fasta") into ch_sample_fasta
  
  script:
  """
  $transdecoder_dir/util/gtf_genome_to_cdna_fasta.pl $sample_gtf $genome_fasta > novel_transcripts.fasta
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
  conda 'bioconda::cpat'
  tag "${hexamer}, ${logit_model}, ${sample_fasta}"

  publishDir "${params.outdir}/${params.name}/cpat/", mode: 'copy'

  input:
  file(hexamer) from ch_hexamer
  file(logit_model) from ch_logit_model_uncompressed
  file(sample_fasta) from ch_sample_fasta_cpat

  output:
  file("${params.name}.ORF_prob.tsv") into ch_cpat_all_orfs
  file("${params.name}.ORF_prob.best.tsv") into ch_cpat_best_orf
  file("${params.name}.ORF_seqs.fa") into ch_cpat_protein_fasta
  file("*")

  script:
  """
  cpat.py \
  -x $hexamer \
  -d $logit_model \
  -g $sample_fasta \
  -o ${params.name} \
  1> ${params.name}_cpat.output \
  2> ${params.name}_cpat.error
  """
}


/*--------------------------------------------------
CDS GTF 
 * create GTF that also contains CDS attributes
---------------------------------------------------*/
process make_cds_gtf {
  cpus 1
  conda 'conda_environments/make_cds.yml'

  publishDir "${params.outdir}/${params.name}/cds/", mode: 'copy'

  input:
    file(sample_gtf) from ch_sample_gtf_cds
    file(cpat_orfs) from ch_cpat_all_orfs //should I make this best orfs instead of all orfs?
  
  output:
    file("${params.name}_with_cds.gtf") into ch_sample_cds
    file("*")
  
  script:
  """
  python cpat_to_gtf.py \
  --sample_gtf $sample_gtf \
  --cpat_orfs $cpat_orfs \
  --output_cds ${params.name}_with_cds.gtf 
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
    conda 'bioconda::agat'

    input:
        file(sample_gtf) from ch_sample_gtf
        file(sample_cds) from ch_sample_cds
    output:
        // file("*")
        file("${params.name}_complete.gtf") into ch_lr_annotation_unformatted
        

    script:
        """
        agat_sp_complement_annotations.pl \
        --ref $sample_gtf \
        --add $sample_cds \
        --out "${params.name}_complete.gtf 
        """
}

/*--------------------------------------------------
VEP pre-formatting
---------------------------------------------------*/
process gtf_for_vep{
    publishDir "${params.outdir}/${params.name}/final_gtf/", mode: 'copy'
    cpus 1
    conda 'bioconda::ensembl-vep'

    input:
        file(complete_gtf) from ch_lr_annotation_unformatted
    output:
      //   file("${params.name}_complete.gtf") into lr_annotation
        file("*")

    script:
        """
        grep -v "#" $complete_gtf | sort -k1,1 -k4,4n -k5,5n -t\$'\t' | bgzip -c > ${complete_gtf}.gz
        tabix -p gff ${complete_gtf}.gz
        """
}

/*--------------------------------------------------
Run VEP to get missense substitutions
 * TODO @ Nuno
---------------------------------------------------*/

// include { run_vep } from '../VEP/workflows/run_vep.nf'

/*--------------------------------------------------
Format missense substitutions to polyphen/sift input
 * TODO @ Nuno @ Renee
 * Example of how Renee normally does this:
   * vep --cache --merged --assembly GRCh38 --cache_version 104 \
   --max_af --sift b --polyphen b --protein --domains \
   --gff LR_transdecoder.gff.gz --pick_allele --pick_order rank,biotype,length \
   --fasta $genomefile --no_stats --offline --fork 4 --force -i $variantfile -o LR_all_patients.tab
   * https://github.com/cmbi/VEP_custom_annotations/blob/main/helper_scripts/parse_vep.py
   * Pros: using pick_allele makes sure I only gather missense substitutions for more severe effect variants, quicker
   * Cons: can't tell if the missense substitutions are the same in hg38 (e.g. A->D & A->D which is not interesting)
   * Ideally would like to run vep with no picking, re-write parse_vep.py to process effects of each variant 
   *  on multiple transcripts and decide whether missense is interesting enough to put through sift/polyphen
---------------------------------------------------*/
// process format_missense{


//        input:
//         file(patient_vcf) from ch_vcf_uncompressed
//        output:
//         file() into inital_vep_output

//     script:
// }

/*--------------------------------------------------
Run polyphen/sift
 * TODO @ Nuno
---------------------------------------------------*/
// include { pph2 } from '../protein_function/nf_modules/run_polyphen2.nf'

/*--------------------------------------------------
Format polyphen/sift output to VEP matrix
 * TODO @ Nuno
 * can potentially re-format Mark's perl script for this, contact Mark for help?
---------------------------------------------------*/
// process parse_patho_pred{

//        input:
//         file(patient_vcf) from ch_vcf_uncompressed
//        output:
//         file() into inital_vep_output

//     script:

// }

/*--------------------------------------------------
Run VEP, adding polyphen scores 
 * TODO @ Nuno
---------------------------------------------------*/
// process final_vep{

//        input:
//         file(patient_vcf) from ch_vcf_uncompressed
//        output:
//         file() into inital_vep_output

//     script:

// }

/*--------------------------------------------------
Last processing step
 * Output 2 files, one with all results and one with filtered results
 * TODO @ Renee
---------------------------------------------------*/
// process gtf_for_vep{

// }

workflow {
   // create transcript fasta 
   create_transcriptome_fasta()
   // do ORF prediction
   cpat(create_transcriptome_fasta.out)
   // make cds gtf out of output
   make_cds_gtf(cpat.out[0])
   // combine input gtf with predicted CDS using agat
   create_final_gtf()
   // format the GTF for VEP
   gtf_for_vep()
   // ...
}

def logHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return 
}





