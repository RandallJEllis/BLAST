# BLAST

A repository for BLAST-related tools.

entrez_blast_homology.py reads a list of gene symbols in a text file, retrieves a FASTA sequence for each gene according to inputted species and database parameters (e.g., gene, protein) in the Entrez.efetch function, queries the FASTA sequence in BLAST according to the inputted parameters in the NcbiblastpComomandline function, and outputs an XML file with the top homologous hit (changing the num_alignments parameter can give you more hits). 

The gene names, accessions for the queried gene and top homolog, and ident values are then all outputted to a CSV file. 

entrez_blast_homology_retroactive.py can be used after retrieving XML outputs from BLAST using entrez_blast_homology.py but failing to store the accessions and ident values to write to a CSV. This way you can separate the tasks of retrieving the BLAST outputs as XML files and parsing the XML files for the top homologs. 

RandallJEllis@gmail.com
