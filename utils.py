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
ps.score_table()


