#!/usr/bin/env bash

for bam in rpmi.flnc.bam saur.flnc.bam calb.flnc.bam; do bedtools bamtofastq -i $bam -fq "${bam%.*}".fastq && minimap2 -t 30 -ax splice:hq -uf --MD /data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/data/hg38.fa "${bam%.*}".fastq > "${bam%.*}".aln.sam; done
for db in rpmi saur calb lps polyic; do python /tools/TranscriptClean/TranscriptClean.py --sam "$db".flnc.aln.sam --genome /data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/data/hg38.fa --spliceJns /data/host_pathogen_PID/pbmc_isoseq/gencode.v29.SJs.tsv --threads 30 --canonOnly --outprefix /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/cleaned/$db; done
for db in rpmi saur calb lps polyic; do mkdir "$db"_tmp && talon_label_reads --f /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/cleaned/"$db"_clean.sam --g /data/host_pathogen_PID/pbmc_isoseq/analysis/swanvis/data/hg38.fa --t 30 --tmpDir "$db"_tmp --deleteTmp --o /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/labeled/"$db"; done
rm /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/gv29.db
talon_initialize_database --f /data/ref_grch38/gencode/gencode.v29.annotation.gtf --g GRCh38 --a GENCODEv29 --idprefix novel --o gv29
talon --f config.csv --db /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/gv29.db --build GRCh38 --cov 0.95 --identity 0.95 --threads 25 --o all_conditions
for db in rpmi saur calb lps polyic; do talon_filter_transcripts --db /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/gv29.db -a GENCODEv29 --datasets $db --minCount 5 --o /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/whitelists_min5/"$db"_whitelist.csv; done
talon_abundance --db /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/gv29.db -a GENCODEv29 -b GRCh38 --o /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/unfiltered
rm /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/whitelists_min5/*
cat /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/whitelists_min5/*_whitelist.csv > /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/whitelists_min5/whitelist.csv
talon_abundance --db /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/gv29.db -a GENCODEv29 -b GRCh38 --whitelist /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/whitelists_min5/whitelist.csv --o /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/min5
talon_create_GTF --db /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/gv29.db -a GENCODEv29 -b GRCh38 --whitelist /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/whitelists_min5/whitelist.csv --observed --o /data/host_pathogen_PID/pbmc_isoseq/analysis/talon/filtered
python /data/host_pathogen_PID/host_pathogen_interactions/reannotation/transcript_dominance.py --ab min5_talon_abundance_filtered.tsv --out min5_talon_abundance_filtered_genetotals.tsv