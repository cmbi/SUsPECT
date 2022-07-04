#!/usr/bin/env nextflow
/*
 * @authors
 * Renee Salz 
 *  
 */

nextflow.enable.dsl=2

params.help = null
params.vep_config = "../VEP/nf_config/vep.ini"
params.vcf = null

def helpMessage() {
    log.info logHeader()
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:

     nextflow run_all.nf \
         --name run_name \
         --outdir outDir \
         --hexamer Human_Hexamer.tsv \
         --logit_model Human_logitModel.RData \
         --talon_gtf filtered_talon_observedOnly.gtf \
         --talon_idprefix talon \
         --genome_fasta hg38.fa.gz \
         --vcf patient.vcf.gz \
         --vep_dir_cache vep_cache \
         --vep_config input/vep.ini \
         --polyphen2_data /path/to/pph2/data

    Other parameters:

      --max_cpus                    Maximum number of CPUs (int)
      --max_memory                  Maximum memory (memory unit)
      --max_time                    Maximum time (time unit)

    See here for more info: https://github.com/cmbi/VEP_custom_annotations
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
log.info "reference_gtf                          : ${params.reference_gtf}"
log.info "hexamer                               : ${params.hexamer}"
log.info "logit_model                           : ${params.logit_model}"
// log.info "max_cpus                              : ${params.max_cpus}"
log.info "name                                  : ${params.name}"
log.info "vcf                                   : ${params.vcf}"
log.info "outdir                                : ${params.outdir}"
log.info "sample_gtf                            : ${params.sample_gtf}"
log.info "talon_gtf                             : ${params.talon_gtf}"
log.info ""


if (!params.sample_gtf && !params.talon_gtf) exit 1, "Must submit sample gtf"


if (!params.reference_gtf && params.sample_gtf) exit 1, "A reference gtf must be provided to determine novelty"

if (params.talon_gtf) {
   ch_talon_gtf=file(params.talon_gtf)
}

if (params.sample_gtf) {
   ch_sample_gtf=file(params.sample_gtf)
   ch_reference_gtf=file(params.reference_gtf)
}

if (!params.hexamer) exit 1, "Cannot find headmer file for parameter --hexamer: ${params.hexamer}"
ch_hexamer = file(params.hexamer)

if (!params.logit_model) exit 1, "Cannot find any logit model file for parameter --logit_model: ${params.logit_model}"

if (!params.vcf) exit 1, "Cannot find any vcf file for parameter --vcf: ${params.vcf}"
ch_vcf = file(params.vcf)

//import all the stuff
include { gunzip_genome_fasta; gunzip_logit_model } from './nf_modules/decompression.nf'
include { identify_novel; filter_novel; clean_gxf } from './nf_modules/find_novel_transcripts.nf'
include { fetch_novel } from './nf_modules/process_talon.nf'
include { convert_to_bed; cpat } from './nf_modules/cpat.nf'
include { cpat_to_bed; combine_bed; bed_to_genepred; genepred_to_gtf; add_genes } from './nf_modules/cpat_to_gtf.nf'
include { gtf_for_vep } from './nf_modules/prepare_for_vep.nf'
include { predict_protein_function } from '../VEP_2_protein_function/main.nf'

workflow full_gtf_input {
   take: 
    full_gtf
    reference_gtf
   main:
    identify_novel(full_gtf,reference_gtf)
    filter_novel(identify_novel.out[0])
    clean_gxf(filter_novel.out)
   emit:
    clean_gxf.out
}

workflow talon_gtf_input {
   take:
    talon_gtf
   main:
    fetch_novel(file(talon_gtf))
   emit:
    fetch_novel.out
}


workflow {
   if (params.talon_gtf) {
      talon_gtf_input(ch_talon_gtf)
      ch_novel_gtf=talon_gtf_input.out
   } else {
      full_gtf_input(ch_sample_gtf,ch_reference_gtf)
      ch_novel_gtf=full_gtf_input.out
   }
   // do ORF prediction
   if (params.logit_model.endsWith('.gz')) 
      ch_logit_model = gunzip_logit_model(file(params.logit_model))
   else
      ch_logit_model = file(params.logit_model)

   convert_to_bed(ch_novel_gtf)
   cpat(ch_hexamer, ch_logit_model, convert_to_bed.out, file(params.genome_fasta))
   // make cds bed out of output
   cpat_to_bed(convert_to_bed.out,cpat.out[0])
   // create whole bed of transcriptome with cds
   combine_bed(convert_to_bed.out,cpat_to_bed.out)
   // convert to gtf for use with vep
   bed_to_genepred(combine_bed.out)
   genepred_to_gtf(bed_to_genepred.out)
   add_genes(genepred_to_gtf.out)
   // format the GTF for VEP
   gtf_for_vep(add_genes.out)
   // run VEP for single aa subs
   predict_protein_function( gtf_for_vep.out,
                            file(params.genome_fasta),
                            file(params.vep_config),
                            cpat.out[2])
}
