"""
Functions to convert genbank file to different formats:
- write DNA sequences of all CDS's to a multi-fasta file
- write protein sequences of all CDS's to a multi-fasta file (two functions, a simple one and one for customized titles)
- convert genbank file to gff format (not only CDS features!)
- convert genbank file to fasta format (not per feature, just entire DNA sequence)

adapted from:
https://bioinformatics.stackexchange.com/questions/4365/how-to-extract-the-protein-fasta-file-from-a-genbank-file
https://www.biostars.org/p/230441/

When these functions are useful: 
    1. NCBI PGAP pipeline provides no list of proteins, only genbank file and gff. 
    2. to customize the headers of multi-fasta files (>....) -> possible to adapt the out.write line
"""

from Bio import SeqIO
from BCBio import GFF
import os
import re

# a loop to convert genbank files in folder to all formats
def convert(folder):
    # loop to collect files to be converted in a list
    filenames = []
    for filename in os.listdir(folder):
        if filename.endswith(".gb") or filename.endswith(".gbk"): # NB: output files will have two .. before their extensions if they were named *.gbk
            filenames.append(filename)
    print(filenames)
    # start a new table to collect metadata (will overwrite any existing file!)
    gbtable = os.path.join(folder, "gbinfotable.txt")
    keylist = ['filename', 'name', 'id', 'description', 'dbxreds', 'molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment']
    with open(gbtable, 'w') as outfile:
        for k in keylist:
            outfile.write(f"{k}\t")
        outfile.write(f"\n")
    # convert gb files to all formats and write metadata to infotable
    # NB: alternative versions if these functions are available below
    for f in filenames:
        print(f"{f}: genbank_to_faa2")
        genbank_to_faa2(f)
        print(f"{f}: genbank_to_ffn")
        genbank_to_ffn(f)
        print(f"{f}: genbank_to_fasta")
        genbank_to_fasta(f)
        print(f"{f}: genbank_to_gff")
        genbank_to_gff(f)
        print(f"{f}: genbank_to_infotable")
        genbank_to_infotable(f, gbtable)

############################### CONVERSION FUNCTIONS ###################################

# for nucleic acid sequences   
def genbank_to_ffn(file_name):
    out = f'{file_name[:-3]}.ffn'
    with open(out, 'w') as out:
        for rec in SeqIO.parse(file_name, "genbank"):
            for feature in rec.features:
                if feature.type == "CDS": #and not feature.qualifiers.get('pseudo'): # if pseudogene than no protein_id is present       
                    geneL = ""
                    gene = feature.qualifiers.get('gene',["unknown"])[0]
                    geneL = f"[gene={gene}] "
                    locus_tagL = ""
                    locus_tag = ""
                    if feature.qualifiers.get('locus_tag'):
                        locus_tag = feature.qualifiers.get('locus_tag',[""])[0]
                        locus_tagL = f"[locus_tag={locus_tag}] "
                    if feature.qualifiers.get('db_xref'):   #in case there is no locus_tag attribute, use db_xref
                        locus_tag = feature.qualifiers.get('db_xref',[""])[0]
                        locus_tagL = f"[locus_tag={locus_tag}] "
                    
                    productL = ""
                    product = feature.qualifiers.get('product',["unknown"])[0]
                    productL = f"[product={product}] "
                    product = product.replace(" ","-") # if using 'product' in the titles before the first space: need to get rid of internal spaces

                    protein_id = feature.qualifiers.get('protein_id',[locus_tag])[0]
                    protein_idL = f"[protein_id={protein_id}]"
             
                    pseudoL = ""
                    if feature.qualifiers.get('pseudo'):
                        pseudoL = f"[pseudo=true] "
                    print(feature.location)
                    p=r'\d+:>*\d+'
                    location = re.findall(p, str(feature.location))[0] # finds the only or first location (in case of a gene spanning de contig ends; two regions are indicated)
                    location = location.replace(':','-')
                    location = location.replace('>','')

                    title = f'{rec.name}_{product}_{location}'
                    out.write(f">{title} {geneL}{locus_tagL}{protein_idL}{pseudoL}{feature.location} {productL} [gbkey=CDS]\n{feature.location.extract(rec).seq}\n") #mimicking 'real' NCBI files
                    #out.write(f">lcl|{rec.name}_cds_{protein_id} {geneL}{locus_tagL}{protein_idL}{pseudoL}{feature.location} {productL} [gbkey=CDS]\n{feature.location.extract(rec).seq}\n") #mimicking 'real' NCBI files

def genbank_to_ffnoriginal(file_name):
    out = f'{file_name[:-3]}.ffn'
    with open(out, 'w') as out:
        for rec in SeqIO.parse(file_name, "genbank"):
            for feature in rec.features:
                if feature.type == "CDS": #and not feature.qualifiers.get('pseudo'): # if pseudogene than no protein_id is present       
                    #print(feature.qualifiers)
                    geneL = ""
                    if feature.qualifiers.get('gene'):
                        gene = feature.qualifiers.get('gene',[""])[0]
                        geneL = f"[gene={gene}] "
                    
                    locus_tagL = ""
                    locus_tag = ""
                    if feature.qualifiers.get('locus_tag'):
                        locus_tag = feature.qualifiers.get('locus_tag',[""])[0]
                        locus_tagL = f"[locus_tag={locus_tag}] "
                        
                    if feature.qualifiers.get('db_xref'):   #in case there is no locus_tag attribute, use db_xref
                        locus_tag = feature.qualifiers.get('db_xref',[""])[0]
                        locus_tagL = f"[locus_tag={locus_tag}] "
                    
                    productL = ""
                    if feature.qualifiers.get('product'):
                        product = feature.qualifiers.get('product',[""])[0]
                        productL = f"[product={product}] "
                        product = product.replace(" ","-") # if using 'product' in the titles before the first space: need to get rid of internal spaces
                    
                    protein_id = locus_tag
                    protein_idL = ""
                    if feature.qualifiers.get('protein_id'):
                        protein_id = feature.qualifiers.get('protein_id',[""])[0]
                        protein_idL = f"[protein_id={protein_id}] "
                    
                    pseudoL = ""
                    if feature.qualifiers.get('pseudo'):
                        pseudoL = f"[pseudo=true] "
                    out.write(f">{rec.name}_{protein_id}_{product}_{feature.location} {geneL}{locus_tagL}{protein_idL}{pseudoL}{feature.location} {productL} [gbkey=CDS]\n{feature.location.extract(rec).seq}\n") #mimicking 'real' NCBI files
                    #out.write(f">lcl|{rec.name}_cds_{protein_id} {geneL}{locus_tagL}{protein_idL}{pseudoL}{feature.location} {productL} [gbkey=CDS]\n{feature.location.extract(rec).seq}\n") #mimicking 'real' NCBI files
       
def genbank_to_fasta(file_name):
    SeqIO.convert(file_name, "genbank", f"{file_name[:-3]}.fasta", "fasta")

# same result as genbank_to_fasta:
def genbank_to_fasta2(file_name):
    out = f'{file_name[:-3]}.fasta'
    for rec in SeqIO.parse(file_name, "genbank"):
        SeqIO.write(rec, out, 'fasta')

# this costum function seems to give the same result as genbank_to_faa()
def genbank_to_faa2(file_name):
    out = f'{file_name[:-3]}.faa'
    with open(out, 'w') as out:
        for rec in SeqIO.parse(file_name, "genbank"):
            for feature in rec.features:
                if feature.type == "CDS": #and not feature.qualifiers.get('pseudo'): # if pseudogene than no protein_id is present       
                    geneL = ""
                    gene = feature.qualifiers.get('gene',["unknown"])[0]
                    geneL = f"[gene={gene}] "

                    locus_tagL = ""
                    locus_tag = ""
                    if feature.qualifiers.get('locus_tag'):
                        locus_tag = feature.qualifiers.get('locus_tag',[""])[0]
                        locus_tagL = f"[locus_tag={locus_tag}] "
                    if feature.qualifiers.get('db_xref'):   #in case there is no locus_tag attribute, use db_xref
                        locus_tag = feature.qualifiers.get('db_xref',[""])[0]
                        locus_tagL = f"[locus_tag={locus_tag}] "

                    productL = ""
                    product = feature.qualifiers.get('product',["unknown"])[0]
                    productL = f"[product={product}] "
                    product = product.replace(" ","-") # if using 'product' in the titles before the first space: need to get rid of internal spaces

                    protein_id = feature.qualifiers.get('protein_id',[locus_tag])[0]
                    protein_idL = f"[protein_id={protein_id}]"
           
                    pseudoL = ""
                    if feature.qualifiers.get('pseudo'):
                        pseudoL = f"[pseudo=true] "
                    print(feature.location)
                    p=r'\d+:>*\d+'
                    location = re.findall(p, str(feature.location))[0] # finds the only or first location (in case of a gene spanning de contig ends; two regions are indicated)
                    location = location.replace(':','-')
                    location = location.replace('>','')
                    #print(location)
                    title = f'{rec.name}_{product}_{location}'
                    #out.write(f">{protein_id} {product}\n{feature.qualifiers.get('translation')[0]}\n") # for short fastsa tiles, like genbank_to_faa()
                    if len(feature.qualifiers.get('translation',[''])[0])>0: # do not write pseudo genes without translation
                        out.write(f">{title} {geneL}{locus_tagL}{protein_idL}{pseudoL}{feature.location} {productL} [gbkey=CDS]\n{feature.qualifiers.get('translation',[''])[0]}\n") # optionally change for different fasta titles (like genbank_to_faa())

def genbank_to_faaoriginal(file_name):
    out = f'{file_name[:-3]}.faa'
    with open(out, 'w') as out:
        for rec in SeqIO.parse(file_name, "genbank"):
            for feature in rec.features:
                if feature.type == "CDS": #and not feature.qualifiers.get('pseudo'): # if pseudogene than no protein_id is present       

                    geneL = ""
                    if feature.qualifiers.get('gene'):
                        gene = feature.qualifiers.get('gene',[""])[0]
                        geneL = f"[gene={gene}] "

                    locus_tagL = ""
                    locus_tag = ""
                    if feature.qualifiers.get('locus_tag'):
                        locus_tag = feature.qualifiers.get('locus_tag',[""])[0]
                        locus_tagL = f"[locus_tag={locus_tag}] "
                    if feature.qualifiers.get('db_xref'):   #in case there is no locus_tag attribute, use db_xref
                        locus_tag = feature.qualifiers.get('db_xref',[""])[0]
                        locus_tagL = f"[locus_tag={locus_tag}] "

                    productL = ""
                    if feature.qualifiers.get('product'):
                        product = feature.qualifiers.get('product',[""])[0]
                        productL = f"[product={product}] "
                        product = product.replace(" ","-") # if using 'product' in the titles before the first space: need to get rid of internal spaces
                    
                    if feature.qualifiers.get('protein_id'):
                        protein_id = feature.qualifiers.get('protein_id',[locus_tag])[0]
                        protein_idL = f"[protein_id={protein_id}]"
          
                    pseudoL = ""
                    if feature.qualifiers.get('pseudo'):
                        pseudoL = f"[pseudo=true] "
                    
                    location = str(feature.location)[1:-4].replace(':','_')
                    print(location)                    
                    #out.write(f">{protein_id} {product}\n{feature.qualifiers.get('translation')[0]}\n") # for short fastsa tiles, like genbank_to_faa()
                    if len(feature.qualifiers.get('translation',[''])[0])>0: # do not write pseudo genes without translation
                        out.write(f">{rec.name}_{product}_{location} {geneL}{locus_tagL}{protein_idL}{pseudoL}{feature.location} {productL} [gbkey=CDS]\n{feature.qualifiers.get('translation',[''])[0]}\n") # optionally change for different fasta titles (like genbank_to_faa())

# for amino acid sequences
def genbank_to_faa(file_name):
    # stores all the CDS entries
    all_entries = []
    with open(file_name, 'r') as gb:
        cdss = SeqIO.InsdcIO.GenBankCdsFeatureIterator(gb)
        for cds in cdss:
            #print(cds)
            if cds.seq is not None:
                all_entries.append(cds)
    SeqIO.write(all_entries, f'{file_name[:-3]}.faa', 'fasta') # short fasta titles
                   
def genbank_to_gff(file_name):
    with open(file_name, 'r') as genbank, open(f'{file_name[:-3]}.gff', 'w') as gffoutput:
        GFF.write(SeqIO.parse(genbank, "genbank"), gffoutput)
    
def genbank_to_infotable(file_name, gbtable):
    #out = f'{file_name[:-3]}.txt'
    dic = {}
    # this keylist should be equal to the one above, where the header of the infotable is written:
    keylist = ['filename', 'name', 'id', 'description', 'dbxreds', 'molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment']
    for k in keylist:
        dic[k] = "" 
    with open(gbtable, 'a+') as outfile:
        for rec in SeqIO.parse(file_name, "genbank"):
            # not sure yet if records always have the properties name, id, description and dbxrefs; perhaps I need try/if statements?
            dic['filename'] = file_name
            dic['name'] = rec.name            
            dic['id'] = rec.id
            dic['description'] = rec.description
            dic['dbxrefs'] = rec.dbxrefs
            dic['molecule_type'] = rec.annotations.get('molecule_type',"")
            dic['topology'] = rec.annotations.get('topology',"")
            dic['data_file_division'] = rec.annotations.get('data_file_division',"")
            dic['date'] = rec.annotations.get('date',"")
            dic['accessions'] = rec.annotations.get('accessions',[""])[0]
            dic['sequence_version'] = rec.annotations.get('sequence_version',"")
            dic['keywords'] = rec.annotations.get('keywords',[""])[0]
            dic['source'] = rec.annotations.get('source',"")
            dic['organism'] = rec.annotations.get('organism',"")
            dic['taxonomy'] = rec.annotations.get('taxonomy',"")
            dic['references'] = rec.annotations.get('references',"")
            dic['comment'] = rec.annotations.get('comment',"").replace("\n"," ")
        for k in keylist:
            outfile.write(f"{dic[k]}\t")
        outfile.write(f"\n")
        print(dic)
