#!/usr/bin/python
#opt/local/bin/python2.7

import PhyloStratiphytUtils as psutils
import TreeTaxonomyUtils as ttutils
import pandas as pd
import numpy as np
from ete2 import NCBITaxa
ncbi = NCBITaxa()


A_mellifera_taxid = 7460
d_mell_taxid = 7227
h_sapiens_taxid = 9606
evalue_cutoff = 1e-3
ref_species = A_mellifera_taxid

# Load taxonomy information
tax_nodes, tax_names = ttutils.load_ncbi_tax_files('ncbi_tax_data/nodes.dmp','ncbi_tax_data/names.dmp')


#h5_store = store_gi_to_taxid_mapping_to_h5(gi2taxid='gi_taxid_prot.dmp.gz')
#gi_2_taxid_table = h5_store.root.gi_2_taxid.table


## Read blast output results
file1 = 'example/blastp_Pbar_ant_vs_nr_e10.txt_nohash'
file2 = 'dysbindin.blast_out.txt_nohash'
file3 = 'dmel-all-translation-r5.44.fasta.nr.e10.sw.txt_50k_nohash' #'dmel-all-translation-r5.44.fasta.nr.e10.sw.txt_nohash'
file4 = 'example/Apis_mellifera.Amel_2.0.13.pep.all.fa.nr.e10.sw.txt'
table = psutils.load_blast_results(file4,evalue_limit=evalue_cutoff)

psutils.store_gi_to_taxid_mapping_to_h5(gi2taxid='gi_taxid_prot.dmp.gz',store_file='gi_taxid.h5', close=True,override=False)
df_unique_gi = psutils.add_taxid_from_h5file(table = table,file = 'gi_taxid.h5')
df_unique_gi = df_unique_gi[df_unique_gi["taxid"] > 0]
# identify LCA of ref specie with every specie with a sequence search hit
lcas = ttutils.get_lca_for_list_with_ETE2( pd.unique(df_unique_gi.taxid),ref_species_id=ref_species)

df_unique_gi = pd.merge(df_unique_gi,lcas,left_on = "taxid",right_on="subject_taxid")

table = pd.merge(table,df_unique_gi,on="subject_gi")

# map each query sequence to the oldest LCA with the target species
# old deprecated method psutils.map_to_oldest_stratum(table)
table_mrca = table.groupby("query_id").apply(lambda x: x.ix[x.lca_rank == np.min(x.lca_rank),["lca_name","lca_tax_id","lca_rank"]].drop_duplicates() ).reset_index().drop("level_1",1) 

# count number of genes by phylostrata
# old deprecated method psutils.count_genes_per_phylostrata(table_mrca,tax_nodes, tax_names,ref_species_id=ref_species)
phylostratum_counts = table_mrca.groupby(["lca_name","lca_tax_id","lca_rank"]).size().reset_index(name="sum").sort("lca_rank")
phylostratum_counts["Phylogenetic_Rank"] = phylostratum_counts.lca_tax_id.apply(lambda x: ncbi.get_rank([x])[x] )

timetree_db = pd.read_table("Hedges_MBE_Supporting_Tables_Excel_rev.txt",sep="\t")
tt_taxids = pd.Series(ncbi.get_name_translator(timetree_db["Node Name"].astype(str))).reset_index(name = "ncbi_taxid")
timetree_db = pd.merge(timetree_db,tt_taxids,left_on="Node Name",right_on="index",how="outer")


phylostratum_counts = pd.merge(phylostratum_counts,timetree_db[["Node Name","Adjusted Time","Study Count","95% CI Lower","95% CI Upper","ncbi_taxid"]],left_on="lca_tax_id",right_on="ncbi_taxid",how="left")

exit(0)
