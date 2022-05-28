#!/usr/bin/env python3

# script adapted from make_pacbio_cds_gtf.py in Long Read Proteogenomics pipeline 
# make GTF file that includes the ORF regions (as CDS features)
# input - gtf of novel seqs ('Novel_ProteinCoding.gtf'), orf calls from CPAT ('Novel_proteincoding.ORF_prob.tsv')
# output - gtf with added "cds" features (orfs)

# %%

import pandas as pd
from collections import defaultdict
import copy
import argparse
import gtfparse
import logging
# %%

def get_first_block_index(orf_coord, cblens):
    # get the index corresponding to the first block containing the orf start
    # return index, and the dela (spacing upstream of end)
    for i, cblen in enumerate(cblens):
        if orf_coord <= cblen:
            delta = cblen - orf_coord
            return i, delta

def make_coords_trimmed_to_orf_range(i1, delta1, i2, delta2, coords):
    orf_coords = copy.deepcopy(coords)
    orf_coords = orf_coords[i1: i2+1]
    # trim ends to orf start/end
    orf_coords[0][0] = orf_coords[0][1] - delta1
    orf_coords[-1][1] = orf_coords[-1][1] - delta2
    return(orf_coords)

def make_coords_trimmed_to_orf_range_neg_strand(i1, delta1, i2, delta2, coords):
    orf_coords = copy.deepcopy(coords)
    orf_coords = orf_coords[i1: i2+1] 
    # trim ends to orf start/end
    orf_coords[0][1] = orf_coords[0][0] + delta1
    orf_coords[-1][0] = orf_coords[-1][0] + delta2
    return(orf_coords)

def make_pacbio_cds_gtf(sample_gtf, cpat_orfs):
    """
    sample_gtf : filename
    cpat_orfs : filename
    """
    # import gtf, only exon info.
    # only move forward with representative pb isoform (for same-protein groups)
    gtf = gtfparse.read_gtf(sample_gtf)
    gtf['transcript_id']=gtf['transcript_id'].str.upper()
    talon_gene = pd.Series(gtf.gene_id.values, index=gtf.transcript_id).to_dict()


    gtf = gtf[['seqname', 'feature', 'start', 'end', 'strand', 'transcript_id']]
    gtf = gtf[gtf['feature'] == 'exon']
    gtf.columns = ['chr', 'feat', 'start', 'end', 'strand', 'acc']


    # pb coords into dict
    pbs = defaultdict(lambda: ['chr', 'strand', [], [], []]) # pb -> [chr, strand, [start, end], [block lengths], [cum. block lengths]]
    # PB.1.1 -> ['chr1', '+', [[100,150], [200,270]], [50, 70], [50, 120], [150-200]]
    for i, row in gtf.iterrows():
        chr, feat, start, end, strand, acc = row
        pbs[acc][0] = chr
        pbs[acc][1] = strand
        pbs[acc][2].append([int(start), int(end)])
    # sort all coords, calc blocks
    def make_cumulative_blens(blocks):
        cblocks = []
        cbl = 0 # cumulative block length
        for b in blocks:
            cbl += b
            cblocks.append(cbl)
        return cblocks
    for acc, infos in pbs.items():
        strand = infos[1]
        if strand == '+':
            infos[2] = sorted(infos[2])
        elif strand == '-':
            infos[2] = sorted(infos[2], reverse=True)
        infos[3] = [end-start+1 for [start, end] in infos[2]]
        infos[4] = make_cumulative_blens(infos[3])


    # read in the ranges of orf on pb transcripts
    ranges = pd.read_table(cpat_orfs)
    ranges= ranges[ranges['Coding_prob']>0.364][['ID','ORF_start','ORF_end']] # human cutoff
    ranges['acc']= ranges['ID'].apply(lambda x: x.split('_')[0])
    ranges['orf_nr']= ranges['ID'].apply(lambda x: x.split('_')[2])


    with open('novel_cds.gtf', 'w') as ofile:
        for i, row in ranges.iterrows():
            oid, orf_start, orf_end, acc, orf_nr = row
            if acc in pbs:
                gene = talon_gene[acc]
                infos = pbs[acc]
                chr, strand, coords, blens, cblens = infos
                # NOTE - uncomment to do chr22-oly test 
                # if chr != 'chr22': continue
                i1, delta1 = get_first_block_index(orf_start, cblens)
                i2, delta2 = get_first_block_index(orf_end, cblens)
                if strand == '+':
                    orf_coords = make_coords_trimmed_to_orf_range(i1, delta1, i2, delta2, coords)
                elif strand == '-':
                    orf_coords = make_coords_trimmed_to_orf_range_neg_strand(i1, delta1,i2, delta2, coords)
                # write out the coordinates
                out_acc = 'gene_id "{}"; transcript_id "{}"; orf_number "{}";'.format(gene, acc, orf_nr)
                for [start, end] in coords:
                    ofile.write('\t'.join([chr, 'hg38_canon', 'exon', str(start), str(end), '.', strand, '.', out_acc]) + '\n')
                for [start, end] in orf_coords:
                    ofile.write('\t'.join([chr, 'hg38_canon', 'CDS', str(start), str(end), '.', strand, '.', out_acc]) + '\n')


def main():
    parser = argparse.ArgumentParser("IO file locations for make pacbio cds gtf")
    parser.add_argument("--sample_gtf", action="store", dest = "sample_gtf") # only observed, novel, pc!
    parser.add_argument("--cpat_orfs", action="store", dest="cpat_orfs") # cpat table output.ORF_prob.tsv
    results = parser.parse_args()
    make_pacbio_cds_gtf(results.sample_gtf, results.cpat_orfs)

if __name__ == "__main__":
    main()
