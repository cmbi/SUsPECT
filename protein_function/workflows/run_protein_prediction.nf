/* 
 * Nextflow pipeline to predict protein function using SIFT and PolyPhen-2
 */

nextflow.enable.dsl=2

// default params
params.help = false
params.outdir = "outdir"
params.singularity_dir = "singularity-images"

// params.gtf="input/annotation/sample.gtf"
params.gtf="input/annotation/sample_human.gtf"
// params.fasta="input/fasta/Canis_lupus_familiarisbasenji.Basenji_breed-1.1.dna_sm.toplevel.fa"
params.fasta="input/fasta/Homo_sapiens.GRCh38.dna_sm.toplevel.fa"
params.program="polyphen2"

// SIFT-specific parameters
params.median_cutoff = 2.75 // as per SIFT's README
params.blastdb="input/sift/uniref100/uniref100_seg.sift.fa"
blastdb_name = file(params.blastdb).name
blastdb_dir = file(params.blastdb).parent

// PolyPhen-2-specific parameters
params.model = "HumDiv.UniRef100.NBd.f11.model"

// module imports
include { getTranslation }            from '../nf_modules/run_agat.nf'
include { sift; alignProteins }       from '../nf_modules/run_sift.nf'
include { pph2; weka }                from '../nf_modules/run_polyphen2.nf'
include { getAminoacidSubstitutions } from '../nf_modules/create_substitutions.nf'

// print usage
if (params.help) {
  log.info """
  Pipeline to predict protein function using SIFT and PolyPhen-2
  --------------------------------------------------------------

  Usage:
    nextflow -C nf_config/nextflow.config run \
             workflows/run_protein_prediction.nf \
             --gtf [path/to/gtf] --fasta [path/to/fasta] \
             -profile lsf -resume

  Options:
    --gtf FILE             Annotation GTF file
    --fasta FILE           Genomic sequence FASTA file
    --program VAL          "sift" (default) or "polyphen2"
    --outdir DIRNAME       Name of output dir. Default: outdir

  SIFT options:
    --blastdb DIR          Path to SIFT-formatted BLAST database, e.g. uniref100
    --median_cutoff VALUE  Protein alignment's median cutoff. Default: 2.75

  PolyPhen-2 options:
    --model VALUE          WEKA model used to calculate PolyPhen-2 scores; can
                           be either "HumDiv.UniRef100.NBd.f11.model" (default)
                           or "HumVar.UniRef100.NBd.f11.model"
  """
  exit 1
}

log.info """
  Starting workflow...
    GTF     : ${params.gtf}
    FASTA   : ${params.fasta}
    program : ${params.program}
"""

workflow run_protein_prediction {
  // Translate transcripts from GTF and FASTA
  getTranslation(file(params.gtf), file(params.fasta))
  translated = getTranslation.out.splitFasta(by: 1, file: true)

  // Generate aminoacid substitutions
  translated_rec = translated.splitFasta(record: [id: true, seqString: true ])
  getAminoacidSubstitutions(translated_rec, params.program)
  subs = getAminoacidSubstitutions.out

  if ( params.program == "sift" ) {
    // Align translated sequences against BLAST database
    alignProteins(translated, blastdb_dir)
    aligned = alignProteins.out

    // Run SIFT
    sift(aligned, subs)
  } else if (params.program == "polyphen2" ) {
    // Run PolyPhen-2
    pph2(translated, subs)
    weka(params.model, pph2.out.results)
  }
}

workflow {
  run_protein_prediction()
}

// Print summary
workflow.onComplete {
    println ( workflow.success ? """
        Workflow summary
        ----------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}
