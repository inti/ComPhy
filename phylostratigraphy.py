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



class phylostatigraphy():
	def __init__(self,ref_specie,alns,evalue_limit=1e-6,h5_file="gi_taxid.h5"):
		from ete_dev import  ncbiquery
		ncbiquery.connect_database()
		import pandas as pd
		#set some defaults
		self.pandas = pd
		self.ncbiquery = ncbiquery
		self.get_sp_lineage = ncbiquery.get_sp_lineage
		self.ref_specie = ref_specie
		self.alignments = alns
		self.alignments_evalue_limit = evalue_limit
		self.h5_file = h5_file
		## compute some initial values
		self.ref_specie_info = { "lineage": self.get_sp_lineage(self.ref_specie),
					 "name": self.ncbiquery.get_taxid_translator([self.ref_specie])[self.ref_specie]
					}
		self.ref_specie_info["ranks"] = dict([(level,i) for i,level in enumerate( self.ncbiquery.translate_to_names(self.ref_specie_info["lineage"]))])

	def pipeline(self):
		# load table with blast results
		self._get_aln_results()
		# find the taxid for each subject_gi
		self._get_unique_subject_gi()
		# get the lineages of the hits
		self._get_subjects_lineages()
		# get LCAs
		self._get_lcas()
			
 	def _get_aln_results(self):
		self.alignments_table = psutils.load_blast_results(self.alignments,evalue_limit=self.alignments_evalue_limit)
	
	def _get_unique_subject_gi(self):
		self.df_subject_taxid = psutils.add_taxid_from_h5file(table = self.alignments_table,file = self.h5_file)

	def _get_subjects_lineages(self):
		self._subject_lineages = [self.get_sp_lineage(taxon) for taxon in self.df_subject_taxid.taxid.unique()]
		
	def _get_lca_ete2(self,lineage1,lineage2):
		lca = -1
		if lineage1[-1] == lineage2[-1]: return lineage1[-1]
		for c,(a,b) in enumerate(zip(lineage1,lineage2)):
    			if a != b:
        			lca = lineage1[c-1]
        			break
		return lca

	def _get_lcas(self):
		subject_lineages_lcas = dict( (subject_taxid,lca) for subject_taxid,lca in  
							zip([l[-1]  for l in self._subject_lineages],
							    [self._get_lca_ete2(self.ref_specie_info['lineage'],l) for l in self._subject_lineages])
						   )
		self.df_subject_taxid['lca'] = self.df_subject_taxid.apply(lambda x: subject_lineages_lcas[x['taxid']],axis=1)
		self.df_subject_taxid['subject_binomial_name'] = self.ncbiquery.translate_to_names(self.df_subject_taxid.taxid.values)
		self.df_subject_taxid['lca_binomial_name'] = self.ncbiquery.translate_to_names(self.df_subject_taxid.lca.values)
		self.df_subject_taxid['lca_ranks'] = self.df_subject_taxid.lca_binomial_name.apply(lambda x: self.ref_specie_info["ranks"][x])
		self.alignments_table = self.pandas.merge(self.alignments_table,self.df_subject_taxid,on="subject_gi")
	
	def score_table(self):
		self.ps_prot_wise_scores = self.pandas.pivot_table(	self.alignments_table,
						rows=["query_id"],
						cols=["lca_binomial_name","lca_ranks","subject_binomial_name"],
						values="bit_score",aggfunc=np.sum)\
						.T.groupby(level=[0,1]).apply(lambda x: np.sum(x)).T


ps = phylostatigraphy(ref_specie=9606,alns="example/dysbindin.blast_out.txt")
ps.pipeline()

# identify LCA of ref specie with every specie with a sequence search hit
lcas = ttutils.get_lca_for_list(list(df_unique_gi['gi']),tax_nodes,tax_names,ref_species_id=ref_species)
table = pd.merge(table,lcas,left_on='taxid',right_on='subject_taxid', how='inner')

# map each query sequence to the oldest LCA with the target species
table_mrca = psutils.map_to_oldest_stratum(table)
# count number of genes by phylostrata
phylostratum_counts = psutils.count_genes_per_phylostrata(table_mrca,tax_nodes, tax_names,ref_species_id=ref_species)
c = psutils.count_genes_per_phylostrata(table_mrca,tax_nodes,tax_names,ref_species_id=ref_species)


exit(0)
