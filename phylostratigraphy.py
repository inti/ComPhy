#!/usr/bin/python
#opt/local/bin/python2.7

import PhyloStratiphytUtils as psutils
import TreeTaxonomyUtils as ttutils
import pandas as pd

d_mell_taxid = 7227
h_sapiens_taxid = 9606
evalue_cutoff = 1e-3
ref_species = d_mell_taxid

# Load taxonomy information
tax_nodes, tax_names = ttutils.load_ncbi_tax_files('ncbi_tax_data/nodes.dmp','ncbi_tax_data/names.dmp')


#h5_store = store_gi_to_taxid_mapping_to_h5(gi2taxid='gi_taxid_prot.dmp.gz')
#gi_2_taxid_table = h5_store.root.gi_2_taxid.table


## Read blast output results
file1 = 'example/blastp_Pbar_ant_vs_nr_e10.txt_nohash'
file2 = 'dysbindin.blast_out.txt_nohash'
file3 = 'dmel-all-translation-r5.44.fasta.nr.e10.sw.txt_nohash'
table = psutils.load_blast_results(file2,evalue_limit=evalue_cutoff)
table['taxid'] = -1

psutils.store_gi_to_taxid_mapping_to_h5(gi2taxid='gi_taxid_prot.dmp.gz',store_file='gi_taxid.h5', close=True,override=False)
df_unique_gi = psutils.add_taxid_from_h5file(table = table,file = 'gi_taxid.h5')
#table=table.ix[table['taxid']>0]
#add_tax_id(table,gi_2_taxid_table)
#load_gi_taxid_mapping('gi_taxid_prot.dmp.gz',table)

# identify LCA of ref specie with every specie with a sequence search hit
lcas = ttutils.get_lca_for_list(list(df_unique_gi['gi']),tax_nodes,tax_names,ref_species_id=ref_species)
table = pd.merge(table,lcas,left_on='taxid',right_on='subject_taxid', how='inner')

# map each query sequence to the oldest LCA with the target species
table_mrca = psutils.map_to_oldest_stratum(table)
# count number of genes by phylostrata
phylostratum_counts = psutils.count_genes_per_phylostrata(table_mrca,tax_nodes, tax_names,ref_species_id=ref_species)
c = psutils.count_genes_per_phylostrata(table_mrca,tax_nodes,tax_names,ref_species_id=ref_species)


exit(0)
