/*
 * From VCF to protein prediction
 */

nextflow.enable.dsl=2

params.help = false
params.outdir = "outdir"
params.gtf = null
params.vep_config = "../VEP/nf_config/vep.ini"
params.model = "HumDiv.UniRef100.NBd.f11.model"

// print usage
if (params.help) {
  log.info """
  Pipeline to run VEP and SIFT/PolyPhen2 based on a VCF file
  ----------------------------------------------------------

  Usage: 
    nextflow run run_all.nf --fasta $fasta --gtf $gtf --vcf $vcf -profile lsf -resume

  Mandatory arguments:
    --vcf VCF     VCF containing missense variants for protein function prediction
    --gtf
    --fasta

  Optional arguments:
    --vep_config  Custom VEP configuration INI file
    --outdir
  """
  exit 1
}

include { run_vep; run_vep as run_vep_plugin } from '../VEP/workflows/run_vep.nf'
include { append_fasta_gtf_to_config; prepare_vep_transcript_annotation; create_exclusion_variants; exclude_pathogenic; filter_common_variants } from '../VEP/nf_modules/utils.nf'
include { pph2; weka } from '../protein_function/nf_modules/run_polyphen2.nf'
include { create_subs } from './nf_modules/create_subs.nf'
include { linearise_fasta; get_fasta } from './nf_modules/fasta.nf'
include { filter_high_severity; dedup_output; map_individuals; add_annotations_to_indiv; get_ref_annotations; combine_custom_ref_candidates } from './nf_modules/output_filtering.nf'

log.info """\
SIFT/PPH2 prediction       v 0.1
================================
fasta    : $params.fasta
gtf      : $params.gtf
vcf      : $params.vcf
outdir   : $params.outdir
"""

workflow predict_protein_function {
  take:
    gtf
    fasta
    vep_config
    protein_fasta
  main:
    // Get translated FASTA
    linearise_fasta( protein_fasta )

    // Filter out common variants
    vep_config_complete = append_fasta_gtf_to_config(vep_config, fasta, gtf)
    run_vep( file( params.vcf ), file( "${params.vcf}.tbi" ), vep_config_complete )
    create_exclusion_variants ( run_vep.out.vcfFile )
    exclude_pathogenic ( run_vep.out.vcfFile, create_exclusion_variants.out )
    filter_common_variants( exclude_pathogenic.out )
    // filter_common_variants( run_vep.out.vcfFile )

    // Get substitutions
    create_subs( filter_common_variants.out.vep_filtered_vcf )
    subs = create_subs.out.flatten()
    get_fasta ( linearise_fasta.out, subs )

    // For each transcript, predict protein function
    pph2(get_fasta.out)
    weka(params.model, pph2.out.results)
    res = weka.out.processed.collectFile(name: "weka_results.out")

    // Incorporate PolyPhen-2 scores into VEP
    prepare_vep_transcript_annotation( res, vep_config_complete, file("../VEP_2_protein_function/VEP_plugins") )
    run_vep_plugin( filter_common_variants.out.originally_benign_af_vcf,
                    filter_common_variants.out.originally_benign_af_vcf_index,
                    prepare_vep_transcript_annotation.out )

    // create output files
    filter_high_severity( run_vep_plugin.out.vcfFile )
    dedup_output( filter_high_severity.out )
    map_individuals(dedup_output.out[0])
    add_annotations_to_indiv(dedup_output.out[0], map_individuals.out)
    get_ref_annotations(dedup_output.out[1])
    combine_custom_ref_candidates(add_annotations_to_indiv.out, get_ref_annotations.out)

}

workflow {
  predict_protein_function( file(params.gtf), file(params.fasta),
                            file(params.vep_config) )
}
