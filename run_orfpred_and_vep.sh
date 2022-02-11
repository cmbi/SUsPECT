#!/usr/bin/env bash

grep 'novel' $talongtf | grep 'protein_coding' > observed_novel_pc.gtf
perl $TransDecoderpath/util/gtf_genome_to_cdna_fasta.pl observed_novel_pc.gtf $genomefile > observed_novel_pc.fasta
perl $TransDecoderpath/util/gtf_to_alignment_gff3.pl observed_novel_pc.gtf > observed_novel_pc.gff3
rm -r observed_novel_pc.fasta.transdecoder_dir*
TransDecoder.LongOrfs -t observed_novel_pc.fasta
cd observed_novel_pc.fasta.transdecoder_dir
pyfasta split -n 40 longest_orfs.pep
rm longest_orfs.pep.flat
rm longest_orfs.pep.gdx
parallel --jobs 20 -I% hmmsearch --cpu 2 --domtblout %.domtblout $pfamhmm % ::: longest_orfs.pep.*
cat *.domtblout > longest_orfs.domtblout
cd ..
TransDecoder.Predict -t observed_novel_pc.fasta --retain_pfam_hits observed_novel_pc.fasta.transdecoder_dir/longest_orfs.domtblout
perl $TransDecoderpath/util/cdna_alignment_orf_to_genome_orf.pl observed_novel_pc.fasta.transdecoder.gff3 observed_novel_pc.gff3 observed_novel_pc.fasta > observed_novel_pc.fasta.transdecoder.genome.gff3
perl $TransDecoderpath/util/gff3_file_to_proteins.pl --gff3 observed_novel_pc.fasta.transdecoder.genome.gff3 --fasta $genomefile > observed_novel_pc.orfs.faa
sed -i 's/*//g' observed_novel_pc.orfs.faa
grep -v "#" observed_novel_pc.fasta.transdecoder.genome.gff3 | awk NF | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > LR_transdecoder.gff.gz
tabix -p gff LR_transdecoder.gff.gz
vep --cache --merged --assembly GRCh38 --cache_version 104 --max_af --sift b --polyphen b --protein --domains --plugin Conservation,$pathtobigwig --gff LR_transdecoder.gff.gz --pick_allele --pick_order rank,biotype,length --fasta $genomefile --no_stats --offline --fork 4 --force -i $variantfile -o LR_all_patients.tab
python /mnt/xomics/renees/data/host_pathogen_PID/host_pathogen_interactions/reannotation/parse_vep.py --vep LR_all_patients.tab