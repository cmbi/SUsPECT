# VEP_custom_annotations
Variant effect prediction based on custom long-read transcriptomes

## Description:

 Variant effects may change with new or updated transcript annotation. This pipeline is an addition to VEP that adds new Polyphen-2 predictions to your variant annotation file.

## Installation:

Running the pipeline requires:
 - conda 
 - a working Polyphen-2 installation (http://genetics.bwh.harvard.edu/pph2/dokuwiki/downloads)
 - pfam hmm file (http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz)
 - a genome fasta file

Please install the conda environment in the 'environments' folder using

```conda env create -f vep_pipeline.yml```

VEP is included in the conda environment. To install the cache, please activate the environment and use the command ```vep_install```.

We recommend separately installing the conservation file (e.g. gerp_conservation_scores.homo_sapiens.GRCh38.bw) from ftp://ftp.ensembl.org/pub/current_compara/conservation_scores/. In case you do not want to do this, you can edit the vep command directly in ```run_orfpred_and_vep.sh```.

If you have multiple patients over multiple vcf files, please combine them into 1 VCF file with vcftools (https://vcftools.github.io/perl_module.html#vcf-merge).

With your combined VCF file (containing multiple patients), please activate the conda environment and run:
```python helper_scripts/vcf_to_ind.py --vcf your_combined.vcf```
this will make a file called 'variants_to_patients.csv' in your current working directory which will be input for the final script in the pipeline.

## Input:

The input for this pipeline is a patient (combined) VCF and a long-read transcriptome. 

The files needed include:
- a gtf of long read transcripts
- a VCF file of patient(s)

## How to run:

### Step 0: Perform TALON analysis

The long-read transcriptome should preferably be processed with TALON (https://github.com/mortazavilab/TALON) before running this pipeline. There is an example pre-processing script in the preprocessing folder. You may also choose to run ```transcript_dominance.py``` if you would like your final output to include transcript dominance information.

### Step 1: ORF prediction & VEP

Activate your vep_pipeline environment by running ```conda activate vep_pipeline```

Do ORF prediction and run VEP ```run_orfpred_and_vep.sh``` script. Run it in the path where you would like your results. Attention: It assumes all of your novel sequence ids are prefeixed with the word 'novel'. 

### Step 2: Polyphen

Run polyphen using the ```polyphen_wrapper.sh``` script:

```bash polyphen_wrapper.sh -s orfs.faa polyphen_in_rare.input```

### Step 3: Generate output

At the moment, the output generation is a simple python script that adds polyphen and GO annotation to the VEP output and outputs CSV format. In future implementation, the output would be formatted and VEP would be run again. (TODO)

```
Usage: python output.py [options]

Options:
    --vep   VEP output tab file with custom annotation enabled
    --pp    polyphen output
    --ab    Abundance file from talon_abundance
    --pat   patient mapping created by function in vcf_to_ind.py (optional)
    --hg38  VEP output tab file with cache hg38 only (optional)
    --go    ensembl to GO mapping file (optional)
    --xn    exon-level abundances made by variants_in_novel_edges.py (optional)
    --out   output file name
```
The output file is a csv file that will depend on what (optional) inputs were provided.
The file will contain all the standard vep columns, as well as information from sources provided to the script such as transcript/exon dominance information, hg38 consequence prediction, polyphen output, GO annotation, etc.
