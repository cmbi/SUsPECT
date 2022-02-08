#!/usr/bin/env bash

conda activate db_creation
cd /mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/transdecoder
grep 'novel' ../talon/filtered_talon_observedOnly.gtf | grep 'protein_coding' > filtered_observed_novel_pc.gtf
perl /mnt/xomics/renees/tools/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl filtered_observed_novel_pc.gtf /mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/data/hg38.fa > filtered_observed_novel_pc.fasta
perl /mnt/xomics/renees/tools/TransDecoder-TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl filtered_observed_novel_pc.gtf > filtered_observed_novel_pc.gff3
rm -r filtered_observed_novel_pc.fasta.transdecoder_dir*
TransDecoder.LongOrfs -t filtered_observed_novel_pc.fasta
cd filtered_observed_novel_pc.fasta.transdecoder_dir
pyfasta split -n 40 longest_orfs.pep
rm longest_orfs.pep.flat
rm longest_orfs.pep.gdx
parallel --jobs 20 -I% hmmsearch --cpu 2 --domtblout %.domtblout /mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/Pfam-A.hmm % ::: longest_orfs.pep.*
cat *.domtblout > longest_orfs.domtblout
cd ..
TransDecoder.Predict -t filtered_observed_novel_pc.fasta --retain_pfam_hits filtered_observed_novel_pc.fasta.transdecoder_dir/longest_orfs.domtblout
perl /mnt/xomics/renees/tools/TransDecoder-TransDecoder-v5.5.0/util/cdna_alignment_orf_to_genome_orf.pl filtered_observed_novel_pc.fasta.transdecoder.gff3 filtered_observed_novel_pc.gff3 filtered_observed_novel_pc.fasta > filtered_observed_novel_pc.fasta.transdecoder.genome.gff3
perl /mnt/xomics/renees/tools/TransDecoder-TransDecoder-v5.5.0/util/gff3_file_to_proteins.pl --gff3 filtered_observed_novel_pc.fasta.transdecoder.genome.gff3 --fasta /mnt/xomics/renees/data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/data/hg38.fa > filtered_observed_novel_pc.orfs.faa
sed -i 's/*//g' filtered_observed_novel_pc.orfs.faa
grep -v "#" filtered_observed_novel_pc.fasta.transdecoder.genome.gff3 | awk NF | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > /mnt/xomics/renees/data/host_pathogen_PID/unsolved_cases/comparison/hg38_vs_custom/five_conditions_transdecoder.gff.gz
