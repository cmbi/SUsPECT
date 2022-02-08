#!/usr/bin/env bash

conda activate VEP
tabix -p gff five_conditions_transdecoder.gff.gz
vep --cache --merged --assembly GRCh38 --cache_version 104 --max_af --sift b --polyphen b --protein --domains --plugin Conservation,/home/renees/.vep/gerp_conservation_scores.homo_sapiens.GRCh38.bw --gff five_conditions_transdecoder.gff.gz --pick_allele --pick_order rank,biotype,length --fasta /mnt/xomics/renees/data/ref_grch38/ensembl/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --no_stats --offline --fork 4 --force -i /mnt/xomics/renees/data/host_pathogen_PID/unsolved_cases/hg38/liftover_vcfs/all_patients_148.vcf -o five_conditions_all_patients.tab
python /mnt/xomics/renees/data/host_pathogen_PID/host_pathogen_interactions/reannotation/parse_vep.py --vep five_conditions_all_patients.tab