# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 15:38:55 2016

@author: ellisrj2
"""
from Bio import Entrez
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import csv

Entrez.email = "randalljellis@gmail.com"
Tool = "SoularBLASTer"

genes = [gene.rstrip('\n') for gene in open('genes.txt')]
accession_numbers = []
homologies = []
for gene in genes:
    handle = Entrez.esearch(db="protein", term="Mus musculus[Orgn] AND " + gene + "[Gene]")
    record = Entrez.read(handle)
    get_fasta = Entrez.efetch(db="protein", id=record["IdList"][0], rettype="fasta", retmode="text")
    filename = record["IdList"][0] + '.fasta'
    out_handle = open(filename, "w")
    out_handle.write(get_fasta.read())
    out_handle.close()
    get_fasta.close()
    
    blastp_cline = NcbiblastpCommandline(query=filename, db="nr", evalue=0.001, outfmt=5, out=filename[:-5] + ".xml", num_alignments=1, entrez_query="human[Organism]", remote=True)
    
    stdout, stderr = blastp_cline()
    
    blastp_cline = open(filename[:-5] + ".xml")
    blast_record = NCBIXML.read(blastp_cline)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            ident = float(hsp.identities)/alignment.length
            accession = alignment.accession
            accession_numbers.append(accession)
            homologies.append(ident)
    
rows = zip(genes, accession_numbers, homologies)
with open('homologies.csv', 'wb') as thefile:
    writer = csv.writer(thefile)
    writer.writerow(['Gene', 'Homolog Accession', 'Homology percentage'])
    for row in rows:
        writer.writerow(row)    
    
    
    