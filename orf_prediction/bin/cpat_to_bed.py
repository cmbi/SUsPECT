#!/usr/bin/env python

# import modules used here
import sys
import argparse
# import pandas as pd
from biopj import relbed, bedio


# Gather our code in a main() function
def main():
    parser = argparse.ArgumentParser(description='Converts CPAT output to fasta.', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_cpat_file', metavar='INPUT_CPAT', type=str,
                                         help='Input file in the CPAT output file format')
    parser.add_argument('input_bed_file', metavar='INPUT_BED', type=str,
                                         help='Transcripts in the BED format')
    parser.add_argument('output_file', metavar='OUTPUT', type=str,
                                         help='The output fasta file')
    # parser.add_argument('-c', dest='cutoff', type=float, default=0.364, 
    #                                      help='The CPAT score cutoff to use')

    args = parser.parse_args()

    out_f = open(args.output_file, 'w')
    # write track line
    out_f.write('track name=CPAT description="CPAT predicted ORFs" useScore=1\n')
    out = bedio.BedFile(out_f)

    print("Loading BED file")
    parser = relbed.RelBEDParser(args.input_bed_file)
    print("Converting CPAT file coordinates")
    with open(args.input_cpat_file) as in_cpat:
        for l in in_cpat:
            p = l.rstrip().split("\t")#.values.flatten().tolist()
            if p[9] != 'Coding_prob' and len(p[5]) > 0:
                tid=p[0].split('_')[0]
                if float(p[9])>0.364:
                    p_bl = bedio.BedLine("{}\t{}\t{}\t{}".format(tid, int(p[4])-1, p[5], p[0]))
                    p_bl_parsed = parser.parse_line(p_bl)
                    p_bl_parsed.score = float(p[8])*1000
                    p_bl_parsed.thick_start = p_bl_parsed.chrom_start
                    p_bl_parsed.thick_end = p_bl_parsed.chrom_end
                    p_bl_parsed.item_rgb = "0,0,0"
                    out.write(p_bl_parsed)
        out.close()

# Standard boilerplate to call the main() function to begin
# the program. Program can be used as module.
if __name__ == '__main__':
    main()
