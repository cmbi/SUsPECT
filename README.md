# VEP custom annotations
Variant effect prediction based on custom long-read transcriptomes

## Description:

 Variant effects may change with new or updated transcript annotation. This Nextflow pipeline runs VEP and appends Polyphen-2 predictions to your variant annotation file using Singularity containers.

## Installation:

Running the pipeline requires:
 - [Nextflow](https://www.nextflow.io)
 - [Polyphen-2 data](http://genetics.bwh.harvard.edu/pph2/dokuwiki/downloads)
 - [VEP cache](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html)
 - [Singularity](https://sylabs.io/singularity/)

We recommend separately installing the conservation file (e.g. gerp_conservation_scores.homo_sapiens.GRCh38.bw) from ftp://ftp.ensembl.org/pub/current_compara/conservation_scores/. In case you do not want to do this, you can edit the vep command directly in ```run_orfpred_and_vep.sh```.

## Input:

The input for this pipeline is a patient (combined) VCF and a long-read transcriptome. 

The files needed include:
- TALON output (more details below)
- a genome FASTA file (bgzipped)
- a GTF of long read transcripts (bgzipped)
- a VCF file of patient(s)

If you have multiple patients over multiple VCF files, please combine them into one VCF file with [vcftools](https://vcftools.github.io/perl_module.html#vcf-merge).

With your combined VCF file (containing multiple patients), please activate the conda environment and run:
```python helper_scripts/vcf_to_ind.py --vcf your_combined.vcf```
this will make a file called 'variants_to_patients.csv' in your current working directory which will be input for the final script in the pipeline.

## How to run:

### Perform TALON analysis

The long-read transcriptome should preferably be processed with [TALON](https://github.com/mortazavilab/TALON) before running this pipeline. There is an example pre-processing script in the preprocessing folder. You may also choose to run ```transcript_dominance.py``` if you would like your final output to include transcript dominance information.

### Run the Nextflow pipeline: ORF prediction, VEP annotation and PolyPhen scores

To start running the Nextflow pipeline, run the following code with your input:

```
cd new_pipeline_dx/orf_prediction
nextflow run_all.nf -resume \
         --name test \
         --outdir testing/ \
         --hexamer input/Human_Hexamer.tsv \
         --logit_model input/Human_logitModel.RData \
         --talon_gtf input/filtered_talon_observedOnly.gtf \
         --talon_idprefix novel \
         --genome_fasta input/hg38.fa.gz \
         --vcf input/homo_sapiens_GRCh38.vcf.gz \
         --vep_dir_cache input/vep_cache
```

The pipeline automatically downloads the required Singularity containers.

A final VCF output will contain VEP annotation with PolyPhen-2 scores.
