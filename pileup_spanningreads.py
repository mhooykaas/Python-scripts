#!/usr/bin/env python
# adapted from pileup_readcounts.py

import re
from Bio import SeqIO
import csv
import os
import argparse
from argparse import RawTextHelpFormatter 

parser = argparse.ArgumentParser(description="""
Created on Thu Dec 29 2018, @author: mhooykaas

Pileup-like script to examine support for a genome assembly sequence by sequencing reads. In contrast to pileup_readcounts.py, count only reads spanning two consecutive nucleotides.
(where in the output e.g. position 5 represents nt 4->5, position 6 represents nt 5->6 etc).
Goal of the script is to reveal (lack of) evidence at certain positions (e.g. where many read alignments suddenly are clipped).
Also checks the start of a contig: how many reads span end-start of contig? (This is expected if contig is circular. If there are 0 reads spanning to position 1, there is either overlap between start and end, or a gap.

Input: reference sequence (multi-fasta format), read alignment file (sam format), path for output  
Output: csv file containing read counts. column 1: contig name; column 2: position within contig, column 3: count (number of reads spanning 'position minus 1' to 'position')
Example call: python pileup_spanningreads.py -ref assemblyx.fasta -sam assemblyx_longreads_minimap2.sam -output spanpileup.csv

Example of counting a very short read aligned to reference: 
          +   +   +   +   +   +   +   +   +                   + -> these positions get +1 in span_pileup csv file (deduced from CIGAR M match, ie including both matches and mismatches between ref and read)
ref   1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,  
      A   G   C   A   T   G   C   A   G   C   G   T   A   G   C
      |   |   |   |   |   |   |   |   |   |               |   |
read  A   G   C   T   T   G   C   A   G   C   -   -   -   G   C
      1,  2,  3,  4,  5,  6,  7,  8,  9,  10                  11  

In order to count position 1 of each contig: count reads aligning to both start and end of contig and check whether there is no additional sequence (insertion) between end and start in read compare to reference
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
spanpileup_dic      = {}     # dictionary with all genome positions + counts

# keep track of potentially end-spanning reads that are clipped at end or start of contig. 
# these reads are aligned at two positions, thus appear twice in sam file
# collect these start- and endreads while looping over sam file and count overlapping ones at the end
endreads            = {} # align at end of contig
startreads          = {} # align at start of contig

# fill up all span pileup dictionary positions with 0's (because final dictionary should cover all positions incl ones without read coverage)
for contig in ref:                          # for all contigs in fasta
    spanpileup_dic[contig.id]=[]
    for j in range(len(contig.seq)):        # at all positions
        spanpileup_dic[contig.id].append(0)
    endreads[contig.name] = []              # need separate list of reads for each contig
    startreads[contig.name] = []

# loop over lines of sam file
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

        qname           = lsplit[0]
        rname           = lsplit[2]
        start_position  = int(lsplit[3])
        cigar           = lsplit[5]
        seq             = lsplit[9]
        
        # if there is no cigar string, skip read (so quit current iteration and continue)
        if cigar == '*':
            continue
        
        # initialize start and end positions of read within reference contig
        # all nucleotide positions between start and end of stretch should get +1 count, excluding first nucleotide OF EACH SEGMENT
        from_position = start_position 
        to_position = start_position    # within loop add length of segment to obtain end position
        read_length_aligned = 0         # I need this to find out which nucleotide of an end-read is the last to match the reference. This is the sum of (mis)matched segments and internal clips/insertions in read compared to reference.
                                        # if read is truly end-spanning, this length is equal to length of part 5' clipped from start-read, because in that case the last end-read nt is adjacent to first start-read nt.

        # I need to know length of contig, to check if read reaches the end -> in that case collect in endreads dictionary:
        for contig in ref:
            if contig.name == rname:        # not sure of SeqIO attribute 'name' is always the same as attribute 'id'
                lencontig = len(contig.seq)

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
        
        # if read aligns to start add to start-read list; may be end-spanning
        if start_position == 1 and pattern_clip.match(all[0]):      # if alignment starts at nt 1 of contig and read is 5' clipped: =start-read. In case of a true end-spanning read, a 'start'-read's first aligning nt is adjacent to 'end'-read's last aligning nt to the contig end. 
            startreads[rname].append((qname,int(all[0][:-1])))      # tuple of read name and length of clipped part; both should match end-read

        # loop over all segments i of CIGAR
        for i in all:
            # if match (or mismatch): progress position and count in pileup
            if pattern_match.match(i):
                read_length_aligned += int(i[:-1])    
                to_position += int(i[:-1])
                # from pos (start) till pos + ref_length: add 1 to pileup
                for p in range(from_position+1,to_position):            # since there is no evidence for spanning to first position (whether first nt of read or first nt after a gap), do not add +1 to first position of each segment
                    spanpileup_dic[rname][p-1] += 1                     # p-1 because Python lists start at 0 whereas sam file positions start at 1
                if p == lencontig and pattern_clip.match(all[-1]):      # if we have reached the last position of the contig this is an end-read and if contig is circular it is likely clipped. Save name and clipped length to compare to startreads list. 
                    endreads[rname].append((qname,read_length_aligned)) # tuple of read name and length of aligned part; should be equal to part of same start-read sticking out at 5' start
                from_position = to_position                             # in the next iteration start from current end position
            
            # if sequence segment is absent from ref: no need to progress position in ref, but do need to add length to read_length_aligned (needed if this is an end-read)
            if pattern_clip.match(i):
                read_length_aligned += int(i[:-1])   
            
            # if deletion in read: progress position in ref, but do not count in pileup
            if pattern_del.match(i):
                to_position += int(i[:-1])
                from_position = to_position                             # in the next iteration start from current end position
            else:
                continue
sam.close()   

# now separately add counts for end-start-spanning reads to pileup dictionary
for contig in ref:
    # per contig: count matching reads in endreads and startreads lists
    endspanning_count = len(set(endreads[contig.name])&set(startreads[contig.name]))
    print(startreads[contig.name])
    print(endreads[contig.name])
    spanpileup_dic[contig.name][0] += endspanning_count

#write lines to csv file
output = open(output,'w', newline='')

csv_writer = csv.writer(output)
for contig in spanpileup_dic:
    for position in range(len(spanpileup_dic[contig])):
        contig = contig
        position_real = position + 1    # correct all positions: positions should start at 1 (Python lists start at 0, so everything is shifted)
        count = spanpileup_dic[contig][position]
        row=[contig, position_real, count]
        csv_writer.writerow(row)

output.close()