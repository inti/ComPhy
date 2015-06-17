import PhyloStratiphytUtils as psutils
import TreeTaxonomyUtils as ttutils
import pandas as pd
import numpy as np

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
			
		self.score_table()
		self.score_weight_table()
		self.phyloStrat_profile(self.ps_prot_wise_W_scores,bootstrap=10)


 	def _get_aln_results(self):
		self.alignments_table = psutils.load_blast_results(self.alignments,evalue_limit=self.alignments_evalue_limit)
	
	def _get_unique_subject_gi(self):
		self.df_subject_taxid = psutils.add_taxid_from_h5file(table = self.alignments_table,file = self.h5_file)
		self.df_subject_taxid = self.df_subject_taxid[ self.df_subject_taxid.taxid >0 ]

	def _get_subjects_lineages(self):
		self._subject_lineages = [ self.get_sp_lineage(taxon) for taxon in self.df_subject_taxid.taxid.unique() ]
		where_is_none = np.where([lineage == None for lineage in self._subject_lineages])[0]
		if np.sum(where_is_none) > 0:
			self._subject_lineages = np.delete(self._subject_lineages, where_is_none )
			self.df_subject_taxid = self.df_subject_taxid.drop(where_is_none)
		
		where_is_not_cellular =  np.where([lineage[1] != 131567 for lineage in self._subject_lineages])[0]
		if np.sum(where_is_not_cellular) > 0:
			self._subject_lineages = np.delete(self._subject_lineages, where_is_not_cellular )
			self.df_subject_taxid = self.df_subject_taxid.drop(where_is_not_cellular)
		self.df_subject_taxid = self.pandas.merge(self.df_subject_taxid, self.pandas.DataFrame( [lineage[-1] for lineage in self._subject_lineages ],columns=["taxid"]
										     ),on="taxid",how="inner")	


	def _get_lca_ete2(self,lineage1,lineage2):
		lca = -1
		if lineage1[-1] == lineage2[-1]: return lineage1[-1]
		diff_len = np.abs(len(lineage1) - len(lineage2))
		if len(lineage2) > len(lineage1):
			if lineage1[-1] == lineage2[-1 + -1*diff_len]: return lineage1[-1]
		if len(lineage1) > len(lineage2):
			if lineage2[-1] == lineage1[-1 + -1*diff_len]: return lineage2[-1]

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
						.T.groupby(level=[0,1])\
						.apply(lambda x: np.sum(x))\
						.reset_index()\
						.sort("lca_ranks",ascending=1)\
						.set_index(["lca_binomial_name","lca_ranks"]).T.fillna(0)
	def score_weight_table(self):
		self.ps_prot_wise_W_scores = self.ps_prot_wise_scores.apply(lambda x: x/np.sum(x),axis=1) 

	
	def phyloStrat_profile(self,table, bootstrap=0, b_fraction=0.99, agg_func="sum",b_n=0):
		self.phyloStrat_profile = table.sum(axis=0)/table.shape[0]
		self.phyloStrat_profile = self.phyloStrat_profile.reset_index(name="phyloStrat_sum")
		if bootstrap > 0:
                        nrow,ncol = (table.shape)
			if b_n > 0:
				if b_n > nrow: 
					print "sorry you want to leave out all the data set!! b_n > nrows"
					return -1
				else:
					b_fraction= 1.0 - np.ceil(1000*float(b_n)/nrow)/1000
					if b_fraction == 1:
						 b_fraction=1.0 - np.floor(1000*float(b_n)/nrow)/1000
			print "Resampling ", bootstrap, "time the", b_fraction,"fraction of the data."
			B = pd.DataFrame(np.zeros((bootstrap,ncol)))
			for i in xrange(bootstrap):
				B.ix[i,:] = self._bootstrap_table(table,b_fraction,agg_func)
			B /= nrow
			self.phyloStrat_profile.ix[:,"phyloStrat_boot_mean"]   = B.mean(axis=0).values
			self.phyloStrat_profile.ix[:,"phyloStrat_boot_q2_5"]   = np.percentile(B,2.5,axis=0)
                        self.phyloStrat_profile.ix[:,"phyloStrat_boot_q15_89"] = np.percentile(B,15.89,axis=0)
                        self.phyloStrat_profile.ix[:,"phyloStrat_boot_q50"]    = np.percentile(B,50.0,axis=0)
                        self.phyloStrat_profile.ix[:,"phyloStrat_boot_q84_1"]  = np.percentile(B,84.1,axis=0)
                        self.phyloStrat_profile.ix[:,"phyloStrat_boot_q97_5"]  = np.percentile(B,97.5,axis=0)

	def _bootstrap_table(self,table,fraction,agg_func="sum"):
		nrow,ncol = (table.shape)
		back = np.zeros((1,ncol))
		n_sampled = 0
		for i in xrange(nrow):
			if np.random.uniform() < fraction:
				back += table.ix[i,:].values 
				n_sampled += 1

		if agg_func == "mean":
			back /= n_sampled
		return back





#ps = phylostatigraphy(ref_specie=9606,alns="example/dysbindin.blast_out.txt")
#ps = phylostatigraphy(ref_specie=144034,alns="example/blastp_Pbar_ant_vs_nr_e10.txt")
ps = phylostatigraphy(ref_specie=44477,alns="example/Apis_mellifera.Amel_2.0.13.pep.all.fa.nr.e10.sw.txt")
ps.pipeline()


