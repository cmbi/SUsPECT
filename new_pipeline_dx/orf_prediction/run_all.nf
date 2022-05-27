#!/usr/bin/env nextflow
/*
 * 
 *   This script was adapted from the "long read proteogenomics" pipeline written by Sheynkman Lab.
 *
 *
 * @authors
 * Renee Salz
 */

nextflow.enable.dsl=2

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
log.info "reference_gtf                          : ${params.reference_gtf}"
log.info "hexamer                               : ${params.hexamer}"
log.info "logit_model                           : ${params.logit_model}"
// log.info "max_cpus                              : ${params.max_cpus}"
log.info "name                                  : ${params.name}"
log.info "vcf                                   : ${params.vcf}"
log.info "outdir                                : ${params.outdir}"
log.info "sample_gtf                            : ${params.sample_gtf}"
log.info ""


if (!params.sample_gtf) exit 1, "Must submit sample gtf"
ch_sample_gtf = Channel.value(params.sample_gtf)

if (!params.reference_gtf) exit 1, "A reference gtf must be provided to determine novelty"
ch_reference_gtf= Channel.value(params.reference_gtf)

if (!params.hexamer) exit 1, "Cannot find headmer file for parameter --hexamer: ${params.hexamer}"
ch_hexamer = Channel.value(params.hexamer)

if (!params.logit_model) exit 1, "Cannot find any logit model file for parameter --logit_model: ${params.logit_model}"

if (!params.vcf) exit 1, "Cannot find any vcf file for parameter --vcf: ${params.vcf}"
ch_normalized_ribo_kallisto = Channel.value(params.vcf)

//import all the stuff
include { gunzip_genome_fasta; gunzip_logit_model } from './nf_modules/decompression.nf'
include { identify_novel; filter_novel } from './nf_modules/find_novel_transcripts.nf'
include { cpat } from './nf_modules/cpat.nf'
include { create_transcriptome_fasta } from './nf_modules/create_transcriptome_fasta.nf'
include { make_cds_gtf; create_final_gtf } from './nf_modules/cpat_to_gtf.nf'
include { gtf_for_vep } from './nf_modules/prepare_for_vep.nf'

workflow {
   // extract only novel sequences
   identify_novel(ch_sample_gtf,ch_reference_gtf)
   filter_novel(identify_novel.out[0])
   // create transcript fasta 
   if (params.genome_fasta.endsWith('.gz')) 
      ch_genome_fasta = Channel.value(gunzip_genome_fasta(params.genome_fasta))
   else
      ch_genome_fasta = Channel.value(params.genome_fasta)
   create_transcriptome_fasta(ch_genome_fasta,filter_novel.out)
   // do ORF prediction
   if (params.logit_model.endsWith('.gz')) 
      ch_logit_model = Channel.value(gunzip_logit_model(params.logit_model))
   else
      ch_logit_model = Channel.value(params.logit_model)
   cpat(ch_hexamer, ch_logit_model, create_transcriptome_fasta.out)
   // make cds gtf out of output
   make_cds_gtf(ch_sample_gtf,cpat.out[0])
   // combine input gtf with predicted CDS using agat
   create_final_gtf(ch_sample_gtf,make_cds_gtf.out)
   // format the GTF for VEP
   gtf_for_vep(create_final_gtf.out)
   // run VEP for single aa subs
   // if (params.vcf.endsWith('.gz')) 
   //    gunzip_vcf(params.vcf)
   // else
   //    params.vcf
   // ...
}