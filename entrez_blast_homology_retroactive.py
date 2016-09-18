# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 11:00:10 2016
@author: ellisrj2
"""
def homology(dirpath):
    import os, time
    from Bio.Blast import NCBIXML
    import csv
    import xml.etree.ElementTree as ET
    
    genes = []
    query_accession_numbers = []
    accession_numbers = []
    homologies = []
    descriptions = []
    # get all entries in the directory w/ stats
    path = os.listdir(dirpath)
    file_list = [entry for entry in path if entry[-3:] == 'xml']
    
    for entry in file_list:
        blastp_cline = open(entry)
        if '_NO_FASTA_SEQ' not in entry:
            blast_record = NCBIXML.read(blastp_cline)
            if len(blast_record.alignments) > 0:
                for alignment in blast_record.alignments:
                    hsp = alignment.hsps[0]
                    ident = float(hsp.identities)/hsp.align_length
                    homolog_accession = alignment.accession
                    accession_numbers.append(homolog_accession)
                    homologies.append(ident)
                    
                    # accession number
                    tree = ET.parse(entry)
                    root = tree.getroot()
                    gene_info = root[5].text
                    idx_first_pipe = gene_info.find('|')
                    idx_second_pipe = min(gene_info[idx_first_pipe+1:].find(i) for i in "|")
                    idx_last_pipe = max(gene_info.rfind(i) for i in "|")
                    descriptions.append(gene_info[idx_last_pipe+2:])
                    accession = gene_info[idx_first_pipe+1:idx_second_pipe+idx_first_pipe+1]
                    query_accession_numbers.append(accession)
                    genes.append(entry[:-4])
            else:
                # accession number
                tree = ET.parse(entry)
                root = tree.getroot()
                gene_info = root[5].text
                idx_first_pipe = gene_info.find('|')
                idx_second_pipe = min(gene_info[idx_first_pipe+1:].find(i) for i in "|")
                idx_last_pipe = max(gene_info.rfind(i) for i in "|")
                descriptions.append(gene_info[idx_last_pipe+2:])
                accession = gene_info[idx_first_pipe+1:idx_second_pipe+idx_first_pipe+1]
                query_accession_numbers.append(accession)
                accession_numbers.append("No homologs")
                homologies.append("No homologs")
                genes.append(entry[:-4])
                continue
                time.sleep(1)
        else:
            tree = ET.parse(entry)
            root = tree.getroot()
            descriptions.append(root[0][0].text)
            query_accession_numbers.append(entry[:-4]) #"No FASTA seq"
            accession_numbers.append(root[0][0].text) #"No FASTA seq"
            homologies.append(root[0][0].text) #"No FASTA seq"
            underscore_index = entry.find('_')
            genes.append(entry[:underscore_index])
            
    rows = zip(genes, query_accession_numbers, descriptions, accession_numbers, homologies)
    with open('GES_Master_homologies1.csv', 'wb') as thefile:
        writer = csv.writer(thefile)
        writer.writerow(['Gene', 'Gene Accession', 'Description',  'Homolog Accession', 'Homology percentage'])
        for row in rows:
            writer.writerow(row)
            