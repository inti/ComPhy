#!/usr/bin/python
#opt/local/bin/python2.7

from PhyloStratiphytUtils import *
import pandas as pd



## Read blast output results
file1 = 'example/blastp_Pbar_ant_vs_nr_e10.txt_nohash'
file2 = 'example/dysbindin.blast_out.txt'
table = load_blast_results(file2)

sequence_length = read_fasta_seqs('example/dysbindin.fa')
table = filter_table_by_ids(table,sequence_length.keys())

table = load_gi_taxid_mapping('gi_taxid_prot.dmp.gz',table)
#table = load_gi_taxid_mapping('head.txt.gz',table)

tree = load_ncbi_taxonomy_tree('/Users/inti/Installed_Apps/biosql/scripts/taxdata/nodes.dmp','/Users/inti/Installed_Apps/biosql/scripts/taxdata/names.dmp')

lcas = get_lca_for_list(list(set(table['taxid'])),tree)
table = pd.merge(table,lcas,left_on='taxid',right_on='subject_taxid', how='inner')


table.ix[(table['evalue'] < 1e-3) & (table['lca_rank'] == min(table['lca_rank'])),'lca_name' ]

exit(0)