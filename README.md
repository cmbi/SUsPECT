# VEP Custom
Variant effect prediction based on custom long-read transcriptomes improves clinical variant annotation

## Description:

 Long-read transcriptomes often contain lots of novelty- especially in samples containing less well-studied tissues/conditions. This has implications for variant effect predictions. Variant effects may change with new or updated transcript annotations. This Nextflow pipeline runs VEP using your transcript file and appends missense effect predictions. The output may help expose potentially disease-causing variants that were previously overlooked.

## Installation:

Running the pipeline requires:
 - [Nextflow](https://www.nextflow.io)
 - [Polyphen-2 data](http://genetics.bwh.harvard.edu/pph2/dokuwiki/downloads)
 - [VEP cache](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html)
 - [Singularity](https://sylabs.io/singularity/)

## Input:

The input for this pipeline is a patient (combined) VCF and a long-read transcriptome. 

The files needed include:
- TALON output (more details below)
- a genome FASTA file (bgzipped)
- a GTF of long read transcripts (bgzipped)
- a VCF file of patient(s)
- Human_logitModel.RData & Human_Hexamer.tsv from https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/

If you have multiple patients over multiple VCF files, please combine them into one VCF file with [vcftools](https://vcftools.github.io/perl_module.html#vcf-merge).


## How to run:

### Perform TALON analysis

The long-read transcriptome should preferably be processed with [TALON](https://github.com/mortazavilab/TALON) before running this pipeline. The idprefix argument given in ```talon_initialize_database``` and the gtf file output from the ```talon_create_GTF``` step are both used as input for this pipeline.

We are currently adapting the pipeline to accept TAMA/SQANTI3 output. Stay tuned.

### Run the Nextflow pipeline: ORF prediction, VEP annotation and PolyPhen scores

To start running the Nextflow pipeline, run the following code with your input:

```
cd new_pipeline_dx/orf_prediction
nextflow run_all.nf \
         --name run_name \
         --outdir outDir \
         --hexamer Human_Hexamer.tsv \
         --logit_model Human_logitModel.RData \
         --talon_gtf filtered_talon_observedOnly.gtf \
         --talon_idprefix talon \
         --genome_fasta hg38.fa.gz \
         --vcf patient.vcf.gz \
         --vep_dir_cache vep_cache
```

The pipeline automatically downloads the required Singularity containers.

A final VCF output will contain VEP annotation with PolyPhen-2 scores.
