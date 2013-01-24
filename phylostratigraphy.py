#!/usr/bin/python
#opt/local/bin/python2.7

from PhyloStratiphytUtils import *
import pandas as pd
import numpy as np
import datetime
from Bio import SeqIO
from cogent.parse.ncbi_taxonomy import NcbiTaxonomyFromFiles

# Read sequences used for the analyses
print return_time(), "Reading sequences used for blast search"
sequence_length = {}
handle = open("example/dysbindin.fa", "rU")
for record in SeqIO.parse(handle, "fasta") :
    sequence_length[record.id] = len(record.seq)
handle.close()


## Read blast output results
print return_time(), "Reading BLAST output file"
# Define the columns names
colnames = ['query_id', 'subject_id', 'identity', 'alignment_length', 'query_length', 'subject_length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']
# read blast output skipping commented lines
file1 = 'example/blastp_Pbar_ant_vs_nr_e10.txt_nohash'
file2 = 'example/dysbindin.blast_out.txt'
table = pd.read_table(file1, "r", delimiter="\t", header=None, names = colnames)
print return_time(), "   '-> [", len(table), "] results read"
# reduce the blast output only to the sequences present on the fasta file
#table = table[table['query_id'].isin(sequence_length.keys()) ]
table['subject_gi'] = table['subject_id'].str.split('|').apply(lambda x: x[1])
table['subject_gi']=table['subject_gi'].astype(np.int64)

print return_time(), "Reading sequence GI to Taxonomy Id mapping file"
tp = pd.read_table('gi_taxid_prot.dmp.gz', chunksize=10000,compression="gzip",header=None,sep="\t",names=['gi','taxid'])

print return_time(), "   '-> Matching blast results GIs with sequences Taxonomy Ids"
gimap = pd.concat([pd.merge(table,chunk,left_on='subject_gi', right_on='gi',how="inner") for chunk in tp], ignore_index=True)
print return_time(), "   '-> [", len(gimap), "] BLAST hits with identified taxid"

print return_time(), "Reading NBCI Taxonomy files"
tree = NcbiTaxonomyFromFiles(open('/Users/inti/Installed_Apps/biosql/scripts/taxdata/nodes.dmp'), open('/Users/inti/Installed_Apps/biosql/scripts/taxdata/names.dmp'))

print return_time(), "Identifying last common ancestor between query specie and blast hits"
lcas = get_lca_for_list(list(set(gimap['taxid'])),tree)
gimap2 = pd.merge(gimap,lcas,left_on='taxid',right_on='subject_taxid', how='inner')
print gimap2[gimap2.columns[[0,1,15,16,17,19,20,21]]][:5]

exit(0)