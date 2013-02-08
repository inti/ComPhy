#!/usr/bin/python
#opt/local/bin/python2.7

import PhyloStratiphytUtils as psutils
import TreeTaxonomyUtils as ttutils
import pandas as pd

p_barbatus_taxid = 144034
d_mell_taxid = 7227
h_sapiens_taxid = 9606
evalue_cutoff = 1e-3
ref_species = p_barbatus_taxid

# Load taxonomy information
tax_nodes, tax_names = ttutils.load_ncbi_tax_files('ncbi_tax_data/nodes.dmp','ncbi_tax_data/names.dmp')


#h5_store = store_gi_to_taxid_mapping_to_h5(gi2taxid='gi_taxid_prot.dmp.gz')
#gi_2_taxid_table = h5_store.root.gi_2_taxid.table


## Read blast output results
file1 = 'example/blastp_Pbar_ant_vs_nr_e10.txt_nohash'
file2 = 'dysbindin.blast_out.txt_nohash'
file3 = 'dmel-all-translation-r5.44.fasta.nr.e10.sw.txt_nohash'
table = psutils.load_blast_results(file3,evalue_limit=evalue_cutoff)

psutils.store_gi_to_taxid_mapping_to_h5(gi2taxid='gi_taxid_prot.dmp.gz',store_file='gi_taxid.h5', close=True,override=False)
df_unique_gi = psutils.add_taxid_from_h5file(table = table,file = 'gi_taxid.h5')

# identify LCA of ref specie with every specie with a sequence search hit
lcas = ttutils.get_lca_for_list(list(set(df_unique_gi['taxid'])),tax_nodes,tax_names,ref_species_id=ref_species)
lcas = pd.merge(df_unique_gi,lcas,left_on='taxid',right_on='subject_taxid')
table = pd.merge(table,lcas,on='subject_gi', how='inner')

# map each query sequence to the oldest LCA with the target species
table_mrca = psutils.map_to_oldest_stratum(table)
# count number of genes by phylostrata
phylostratum_counts = psutils.count_genes_per_phylostrata(table_mrca,tax_nodes, tax_names,ref_species_id=ref_species)
print phylostratum_counts
b = psutils.bootstrap(table,B=100)
b

