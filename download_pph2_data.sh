#!/usr/bin/env bash

wget http://genetics.bwh.harvard.edu/pph2/dokuwiki/_media/polyphen-2.2.2r405c.tar.gz
tar xf polyphen-2.2.2r405c.tar.gz
wget ftp://genetics.bwh.harvard.edu/pph2/bundled/polyphen-2.2.2-databases-2011_12.tar.bz2
wget ftp://genetics.bwh.harvard.edu/pph2/bundled/polyphen-2.2.2-alignments-mlc-2011_12.tar.bz2
wget ftp://genetics.bwh.harvard.edu/pph2/bundled/polyphen-2.2.2-alignments-multiz-2009_10.tar.bz2
tar xf polyphen-2.2.2-databases-2011_12.tar.bz2
tar xf polyphen-2.2.2-alignments-mlc-2011_12.tar.bz2
tar xf polyphen-2.2.2-alignments-multiz-2009_10.tar.bz2
cd polyphen-2.2.2
wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.2.31/ncbi-blast-2.2.31+-x64-linux.tar.gz
tar xf ncbi-blast-2.2.31+-x64-linux.tar.gz
mv ncbi-blast-2.2.31+/* blast
rm -rf ncbi-blast-2.2.31+*
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat -O bin/blat
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa -O bin/twoBitToFa
chmod +x bin/blat bin/twoBitToFa
cd nrdb
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref100/uniref100.fasta.gz
gunzip uniref100.fasta.gz
../update/format_defline.pl uniref100.fasta >uniref100-formatted.fasta
../blast/bin/makeblastdb -in uniref100-formatted.fasta -dbtype prot -out uniref100 -parse_seqids
rm -f uniref100.fasta uniref100-formatted.fasta
cd ..
rsync -rltv --delete-after --port=33444 rsync.wwpdb.org::ftp/data/structures/divided/pdb/ wwpdb/divided/pdb/
rsync -rltv --delete-after --port=33444 rsync.wwpdb.org::ftp/data/structures/all/pdb/ wwpdb/all/pdb/
sleep 8h && rsync -rltvz --delete-after rsync://rsync.cmbi.ru.nl/dssp/ dssp/
cd ..
sudo mkdir /export/apps/polyphen
sudo chown aorth:aorth /export/apps/polyphen
cp -r polyphen-2.2.2 /export/apps/polyphen/2.2.2
export PPH=/export/apps/polyphen/2.2.2
cd $PPH/src
make download
make clean
make
make install
cd $PPH
./configure