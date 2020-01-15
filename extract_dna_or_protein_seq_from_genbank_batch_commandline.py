#!/usr/bin/env python
"""
Created 2019, @author: Marjolein Hooykaas
Purpose: Extracts DNA or protein sequences from genbank files based annotation (gene name search terms). Then appends these to a multi-fasta file. Long (incl or excl genbank filename) or short names can be written.
Instructions: Save all genbank files to extract sequences from in one folder, choose output type (DNA/protein), optionally specify details about names to be used in fasta file

NB: There is variation in genbank file formatting. Eg, Rast annotations in genbank format have accession numbers (same as 'locus' (contigs names)). Prokka genbank annotations contains only 'locus', but accession field is empty and 'locus' does not seem to be a dictionary key
Rast genbank feature keys ['db_xref', 'translation', 'product', 'transl_table'] or ['db_xref', 'translation', 'product', 'EC_number', 'transl_table']
Prokka genbank feature keys: ['locus_tag', 'inference', 'codon_start', 'transl_table', 'product', 'translation'],['organism', 'mol_type', 'strain'] (list 1 for true gene features, list2 for contigs?)

example of genbank structure:
     gene            complement(943..1656)
                     /locus_tag="Atu3001"
                     /old_locus_tag="AGR_L_3592"
     CDS             complement(943..1656)
                     /locus_tag="Atu3001"
                     /old_locus_tag="AGR_L_3592"
                     /codon_start=1
                     /transl_table=11
                     /product="hypothetical protein"
                     /protein_id="AAK90379.1"
                     /translation="MKYIFLLSAFLIFLFSAFIYFFPIFGAVVCPPCFGFVEVEGGIY
                     IDASLNADLEKDKILSNLDAAHLLLREVYGEVEAPLPAIFLCVSKNCATYLGRRGEKA
                     SSFGHWAIVVYRDGNNSGILAHELSHIEIGFRLGFYNMSSVPIWFDEGVAVVASQDRR
                     YLNVDPSGRLSCKEGVSGAVIADLDEWWRRASIGDVGIYSAAACEVMKWMDRRGNEGS
                     LVRLLDLLRSGESFDVAFE"
example 2:
     CDS             13958..14305
                     /db_xref="SEED:fig|358.213.peg.13"
                     /translation="MIGSSGNVRVYLACGVTDMRRGIDGLSALVETVVKEAPGSGAIF
                     GFRGKRADRIKLLWWDGQGFCLFYKILERGYFPWPTAKEGVAHLTQAQLSMLVEGIDW
                     RRPAWTSAPDRTG"
                     /product="Mobile element protein"
                     /transl_table=11

example 3:
     gene            complement(790..2031)
                     /gene="acs"
                     /locus_tag="pTiBo001"
     CDS             complement(790..2031)
                     /gene="acs"
                     /locus_tag="pTiBo001"
                     /note="agrocinopine synthase"
                     /codon_start=1
                     /transl_table=11
                     /product="Acs"
                     /protein_id="AAZ50392.1"
                     /translation="MWELEWDLPAGTSVSEVLARYSTPNLLQKLDEKLDVQVVEHRGM
                     FNLGKGIQECTKTAILSAIGEGHRNLCEIDIALTADGVPIVAHEFNVFRVAALDEDKP
"""

from Bio import SeqIO
import os
import argparse
from argparse import RawTextHelpFormatter 

parser = argparse.ArgumentParser(description="""
Created in 2019, by Marjolein Hooykaas

Extracts DNA or protein sequences of specified [genes] from genbank files. Writes a multi-fasta file. If file already exists, appends the sequences to this file.

Example call: $ python extract_dna_or_protein_seq_from_genbank_batch_commandline.py -folder /d/genbankfiles/ -genes ["RecA","recombinase A"] -format DNA --filename
""", formatter_class=RawTextHelpFormatter)
parser.add_argument('-folder',    help='location of genbank files to be parsed',   required  = True)
#parser.add_argument('-genes',    help='provide bracketed list of gene names to search by (NB: list should contain >1 name), eg ["RecA","recombinase A"]', required  = True)
parser.add_argument('-genes',    help='provide text file containing a list of gene names to search for (NB: list should contain >1 name), eg RecA, recombinase A (on different lines)', required  = True)
parser.add_argument('-format', choices=['DNA', 'protein'], help='output format: DNA or protein', required  = True)
parser.add_argument('-f', '--filename', action='store_true', help='whether to include genbank filename in fasta names (not compatible with -s)')
parser.add_argument('-s' , '--short_names', action='store_true', help='whether to make fasta names extra short')
args = vars(parser.parse_args())

folder    = args['folder']
genes    = args['genes']
#genes = genes[1:-1].split(',') # this does not always work, therefore switch to reading in a file with a list of names
# removed help='provide bracketed list of gene names to search by (NB: list should contain >1 name), eg ["RecA","recombinase A"]', required  = True)
genes = [line.rstrip('\n') for line in open(genes)]
print(genes)

write_type = args['format']

os.chdir(folder)
output_file_name = f"{genes[0]}.fasta" # generate named based on first gene name in list: append later runs on other strains to this file
output_file_name = "".join([x if x.isalnum() or x in "._-()" else "_" for x in output_file_name]) # filter out unwanted characters from name

if not os.path.exists(f"{folder}/{write_type}/"):
    os.mkdir(f"{folder}/{write_type}/")
    
if os.path.exists(f"{folder}/{write_type}/{output_file_name}"):
    print(f"Filename {output_file_name} already exists; new sequences will be appended to file")

output = f"{folder}/{write_type}/{output_file_name}"

genes = [gene.lower() for gene in genes] # better chance for a hit when search is case insensitive -> change both query and subject to lowercase
print(f"search for: ", genes)

addition = '' # if file name should not be in fasta name lines, write empty string (for all output)

for input_file in os.listdir(folder):
    if input_file.endswith(('.gbff','gb', 'gbk')):
        print(f"\nStart searching in {input_file}")
        print(f"if successful append sequence to {output}")
        number_of_hits = 0 # per genbank file start at 0 again
        featuresFound = set()
        #if include_filename == 'TRUE':      # for concatenating alignments, it may be convenient to include file names in fasta names, to match sequences to each other
        if args['filename']:
            addition = f"{input_file} | "
        for rec in SeqIO.parse(input_file, "gb"): # calls the record for the genbank file and SeqIO (BioPython module) to parse it
            try:
                organism = rec.annotations['organism'] # defines your organism ID
            except:
                organism = 'unknown organism'
            try:
                acc = rec.annotations['accessions'][0] # accession numbers given in Rast genbank, but not prokka genbank -> try:except statement required
            except:
                acc = 'contig xx'
            for feature in rec.features: # parse all features in the genbank file
                if feature.type in ['CDS','rRNA','tRNA']:  # only check CDS, rRNA, tRNA features and skip 'gene' because it seems to be redundant (and never contains translation key)
                    for gene in genes: # look for all provided gene names    
                        descriptions = [feature.qualifiers.get('gene',[""])[0],feature.qualifiers.get('product',[""])[0]] # optionally change product to note #either in gene or product field; however, one or both fields may not exist -> with dict.get method avoid many if statements, provide default ("") in case of not found
                        for description in descriptions:
                            #if gene in description.lower() and (not ("like" or "family" or "superfamily") in description): # try to filter out bad hits by discarding those with 'like'or 'family' (assume these hits are similar but not orthologs) 
                            if gene in description.lower() and not(True in [x in description.lower() for x in ["like","family"]]): #try to filter out bad hits by discarding those with 'like'or 'family' (assume these hits are similar but not orthologs) 
                                print(f'gene found: {gene} in {description}, write to fasta:')
                                if write_type == 'DNA':                         
                                    fastaContent = feature.extract(rec.seq)
                                else:
                                    fastaContent = feature.qualifiers.get('translation',[""])[0]
                                with open(output, "a") as ofile: # opens the output file and "a" designates it for appending
                                    geneIDs = set([feature.qualifiers.get('locus_tag',['ID_not_found'])[0],feature.qualifiers.get('db_xref',['ID_not_found'])[0]]) # most genbank files seen contain either locus_tag or db_xref keys for gene identifiers, but eg C58 contains both; because data type is set, will never contain 2x 'ID not found'
                                    geneIDs.discard('ID_not_found') 
                                    if len(featuresFound & geneIDs)==0: # if no overlap between features already found and currente geneID: write sequence
                                        print(f">{geneIDs} | {organism} | accession_{acc} | {description}")
                                        if not args['short_names']:
                                            ofile.write(f">{addition}{list(geneIDs)[0]} | {organism} | accession_{acc} | {description}\n{fastaContent}\n\n")
                                        else:
                                            ofile.write("_".join(f">{list(geneIDs)[0]}_{organism.split()[-1]}_{description}".split()))
                                            ofile.write(f"\n{fastaContent}\n\n")
                                        number_of_hits += 1
                                        featuresFound.update(geneIDs)
                                        featuresFound.discard('ID_not_found') # purge 'ID_not_found', otherwise most further features will not be written out due to overlap over 'ID_not_found'
        print(f"\nNumber of sequences written for this genome: {number_of_hits}") # for each each genome print if gene was found