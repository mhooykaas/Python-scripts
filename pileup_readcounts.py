#!/usr/bin/env python

import re
from Bio import SeqIO
import csv
import os
import argparse
from argparse import RawTextHelpFormatter 

parser = argparse.ArgumentParser(description="""
Created on Thu Dec 20 2018, @author: mhooykaas

Script counts read coverage at each location in genome.
Similar information can be obtained with samtools mpileup but mpileup by default filters out certain alignments based on flags (and -ff flag to undo this did not work for me), whereas this script keeps all (incl eg secondary alignments)

Input: reference sequence (multi-fasta format), read alignment file (sam format), path for output 
Output: csv file containing read counts
Example call: python pileup_readcounts.py -ref assemblyx.fasta -sam assemblyx_longreads_minimap2.sam -output longreads_pileup.csv

Example of scoring:
      +   +   +   +   +   +   +   +   +   +               +   + -> these positions get +1 in span_pileup csv file (deduced from CIGAR M match, ie including both matches and mismatches between ref and read)
ref   1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,  
      A   G   C   A   T   G   C   A   G   C   G   T   A   G   C
      |   |   |   |   |   |   |   |   |   |               |   |
read  A   G   C   T   T   G   C   A   G   C   -   -   -   G   C
      1,  2,  3,  4,  5,  6,  7,  8,  9,  10                  11  

""", formatter_class=RawTextHelpFormatter)
parser.add_argument('-ref',    help='provide reference sequence in fasta format',   required  = True)
parser.add_argument('-sam',    help='provide alignment file in sam format',         required  = True)
parser.add_argument('-output', help='provide path for output',                      required  = True)
args = vars(parser.parse_args())

ref    = args['ref']
sam    = args['sam']
output = args['output']

sam                 = open(sam)
ref                 = list(SeqIO.parse(ref, "fasta"))
pileup_dic          = {}     # dictionary with all genome positions + counts

# fill up all pileup dictionary positions with 0's
for contig in ref:     # for all contigs in fasta
    pileup_dic[contig.id]=[]
    for j in range(len(contig.seq)):     # at all positions
        pileup_dic[contig.id].append(0)

for l in sam:
    l = l.strip()   # remove line endings
    if not l.startswith("@"):   # ignore headings; process lines describing the alignments
        # split columns of this line into list
        lsplit = l.split("\t")
        # we could extract all columns but we only need contig name (rname), position and cigar string, and for endreads: qname
        # [0]  qname
        # [1]  flag
        # [2]  rname
        # [3]  start position
        # [4]  maping quality
        # [5]  cigar string
        # [6]  rnext
        # [7]  pnext
        # [8]  tlen
        # [9]  sequence
        # [10] quality

        rname           = lsplit[2]
        start_position  = int(lsplit[3])
        cigar           = lsplit[5]

        # if there is no cigar string, skip read (so quit current iteration and continue)
        if cigar == '*':
            continue
        
        # initialize start and end positions of read within reference contig
        # all nucleotide positions between start and end of stretch should get +1 count, excluding first nucleotide OF EACH SEGMENT but including end position(s)
        from_position = start_position 
        to_position = start_position    # within loop add length of segment to obtain end position

        # process cigar string 
        # cigarstring: length of segment + one of options: MX=DNISHP:
        # M	Match; can be either an alignment match or mismatch! The nucleotide is present in the reference.
        # X	Sequence Mismatch; the nucleotide is present in the reference
        # =	Sequence Match; the nucleotide is present in the reference
        # D	Deletion; the nucleotide is present in the reference but not in the read
        # N	Skipped region; a region of nucleotides is not present in the read
        # I	Insertion; the nucleotide is present in the read  but not in the rference.
        # S	Soft Clipping;  the clipped nucleotides are present in SEQ
        # H	Hard Clipping; the clipped nucleotides are not present in SEQ
        # Hard masked bases do not appear in the SEQ string, soft masked bases do
        # P	padding (silent deletion from padded reference)
        
        '''
        example:
        cigar = '10M2D10M4I20M2D10M'
        there are 7 segments; will loop through all of them, testing whether segment matches any of below regular expression patterns
        if match or mismatch: loop through all positions in this DNA stretch and add 1 to spanpileup dictionary except all first positions of segment
        if deletion in read: shift to new position in reference but do not count in span pileup dictionary
        '''
        # compile regular expressions:
        pattern_all = re.compile("[0-9]+[MX=DNISHP]")       # all, to count how many segments alignment consists of
        pattern_match = re.compile("[0-9]+[MX=]")           # (mis)matched: progress ref, also count in pileup
        pattern_del = re.compile("[0-9]+[DN]")              # deletion from read: progress ref but do not count in pileup
        pattern_clip = re.compile("[0-9]+[ISH]")            # insertion in read: do not progress ref, do not count in pileup
        # do not include P: Padding; not sure what it means and I did not see it in my own sam files so far
        # find all 'segments' in cigar string (combinations of number + character)
        all = pattern_all.findall(cigar)
        
        for i in all:
            # if match (or mismatch): progress position and count in pileup
            if pattern_match.match(i):
                to_position += int(i[:-1])
                # from pos (start) till pos + ref_length: add 1 to pileup
                for p in range(from_position,to_position):      # since there is no evidence for spanning to first position (whether first nt of read or after gap), do not add +1 to first position of each segment
                    pileup_dic[rname][p-1] += 1                 # p-1 because Python lists start at 0 whereas sam file positions start at 1
                from_position = to_position
            
            # if deletion in read: progress position in ref, but do not count in pileup
            if pattern_del.match(i):
                to_position += int(i[:-1])
                from_position = to_position
            else:
                continue
sam.close()   

#write lines to csv file
output = open(output,'w', newline='')

csv_writer = csv.writer(output)
for contig in pileup_dic:
    for position in range(len(pileup_dic[contig])):
        contig = contig
        position_real = position + 1    # correct all positions: positions should start at 1 (Python lists start at 0, so everything is shifted)
        count = pileup_dic[contig][position]
        row=[contig, position_real, count]
        csv_writer.writerow(row)

output.close()