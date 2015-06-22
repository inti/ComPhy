
import pandas as pd
import numpy as np
from ete2 import NCBITaxa
ncbi = NCBITaxa()
from ..utils import ttutils
from ..utils import psutils

class phylostratigraphy():
	def __init__(self,ref_specie,alns,evalue_limit=1e-6,h5_gi_taxid_file="../dbs/gi_taxid.h5",n_bootstrap = 0, fraction_rows_to_bootstrap= 0.1, n_rows_to_bootstrap = 0, ancestor_times = False):
		#set some defaults
		self.pandas = pd
		self.ncbi_taxonomy = ncbi
		self.alignments = alns
		self.alignments_evalue_limit = evalue_limit
		self.h5_gi_taxid_file = h5_gi_taxid_file
		## compute some initial values
		self.ref_specie_info = { "lineage": 0, "name": 0 , "taxid" : ref_specie}
		self.ref_specie_info["name"],self.ref_specie_info["lineage"] = ttutils._get_commonName_and_lineage(self.ref_specie_info["taxid"] , self.ncbi_taxonomy)
		self.n_bootstrap = n_bootstrap
		self.fraction_rows_to_bootstrap = fraction_rows_to_bootstrap
		self.n_rows_to_bootstrap = n_rows_to_bootstrap
		self.ancestor_times = ancestor_times


	def pipeline(self):
		# load table with blast results
		self._get_aln_results()
		# find the taxid for each subject_gi
		self._get_unique_subject_gi()
		# get LCAs
		self._get_lca_ete2()
		self._phyloStrat_table()

		#self.score_table()
		#self.score_weight_table()
		if self.n_bootstrap > 0:
			self.phyloStrat_profile(self.phyloStrat_table_wide,bootstrap=self.n_bootstrap)
		if self.ancestor_times:
			self._get_ancestor_times()		

 	def _get_aln_results(self):
		self.alignments_table = psutils.load_blast_results(self.alignments,evalue_limit=self.alignments_evalue_limit)
	
	def _get_unique_subject_gi(self):
		self.df_subject_taxid = psutils.add_taxid_from_h5file(table = self.alignments_table,file = self.h5_gi_taxid_file)
		self.df_subject_taxid = self.df_subject_taxid[ self.df_subject_taxid.ncbi_taxid >0 ]

	def _get_lca_ete2(self):
		self.hits_lcas = ttutils.get_lca_for_list_with_ETE2( self.pandas.unique( self.df_subject_taxid.ncbi_taxid),self.ref_specie_info["taxid"],self.ncbi_taxonomy)
		self.df_subject_taxid  = self.pandas.merge( self.df_subject_taxid,self.hits_lcas,left_on = "ncbi_taxid",right_on="subject_taxid")
		self.alignments_table = self.pandas.merge(self.alignments_table, self.df_subject_taxid ,on="subject_gi")

	def _phyloStrat_table(self):
		self.phyloStrat_table = self.alignments_table.groupby("query_id").apply(lambda x: x.ix[x.lca_rank == np.min(x.lca_rank),["lca_name","lca_tax_id","lca_rank"]].drop_duplicates() ).reset_index().drop("level_1",1)
		self.phyloStrat_table_wide = self.pandas.pivot_table(self.phyloStrat_table,index=["query_id"],columns=["lca_tax_id"])
		self.phyloStrat_table_wide[np.isnan(self.phyloStrat_table_wide) == False ] = 1
		self.phyloStrat_table_wide = self.phyloStrat_table_wide.fillna(0)
		self.phyloStrat_table_summary = self.phyloStrat_table.groupby(["lca_name","lca_tax_id","lca_rank"]).size().reset_index(name="sum").sort("lca_rank")
		self.phyloStrat_table_summary["Phylogenetic_Rank"] = self.phyloStrat_table_summary.lca_tax_id.apply(lambda x: self.ncbi_taxonomy.get_rank([x])[x] )
	
#	def _phyloStrat_score_table(self):
#		self.ps_prot_wise_scores = self.pandas.pivot_table(	self.alignments_table,
#						rows=["query_id"],
#						cols=["lca_binomial_name","lca_ranks","subject_binomial_name"],
#						values="bit_score",aggfunc=np.sum)\
#						.T.groupby(level=[0,1])\
#						.apply(lambda x: np.sum(x))\
#						.reset_index()\
#						.sort("lca_ranks",ascending=1)\
#						.set_index(["lca_binomial_name","lca_ranks"]).T.fillna(0)
#
#	def score_weight_table(self):
#		self.ps_prot_wise_W_scores = self.ps_prot_wise_scores.apply(lambda x: x/np.sum(x),axis=1) 

	
	def phyloStrat_profile(self,table, bootstrap=0, b_fraction=0.99, agg_func="sum",b_n=0):
		self.phyloStrat_table_wide_bootstrap = table.sum(axis=0)/table.shape[0]
		self.phyloStrat_table_wide_bootstrap = self.phyloStrat_table_wide_bootstrap.reset_index(name="phyloStrat_sum")
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
			self.phyloStrat_table_wide_bootstrap.ix[:,"phyloStrat_boot_mean"]   = B.mean(axis=0).values
			self.phyloStrat_table_wide_bootstrap.ix[:,"phyloStrat_boot_q2_5"]   = np.percentile(B,2.5,axis=0)
                        self.phyloStrat_table_wide_bootstrap.ix[:,"phyloStrat_boot_q15_89"] = np.percentile(B,15.89,axis=0)
                        self.phyloStrat_table_wide_bootstrap.ix[:,"phyloStrat_boot_q50"]    = np.percentile(B,50.0,axis=0)
                        self.phyloStrat_table_wide_bootstrap.ix[:,"phyloStrat_boot_q84_1"]  = np.percentile(B,84.1,axis=0)
                        self.phyloStrat_table_wide_bootstrap.ix[:,"phyloStrat_boot_q97_5"]  = np.percentile(B,97.5,axis=0)

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

	def get_ancestor_times(self,timetreedb):
		self.timeTree_db = self.pandas.read_table(timetreedb,sep="\t")
		tt_taxids = self.pandas.Series(self.ncbi_taxonomy.get_name_translator(self.timeTree_db["Node Name"].astype(str))).reset_index(name = "ncbi_taxid")
		self.timeTree_db = pd.merge(self.timeTree_db,tt_taxids,left_on="Node Name",right_on="index",how="outer")


class ComPhy():
        def __init__(self, muti_species_file, evalue_limit=1e-6, h5_gi_taxid_file="dbs/gi_taxid.h5", n_bootstrap = 0, fraction_rows_to_bootstrap= 0.1, n_rows_to_bootstrap = 0, ancestor_times = False):
                #set some defaults
                self.pandas = pd
                self.ncbi_taxonomy = ncbi
                self.alignments_evalue_limit = evalue_limit
                self.h5_gi_taxid_file = h5_gi_taxid_file
                ## compute some initial values
                self.n_bootstrap = n_bootstrap
                self.fraction_rows_to_bootstrap = fraction_rows_to_bootstrap
                self.n_rows_to_bootstrap = n_rows_to_bootstrap
                self.ancestor_times = ancestor_times
		self.muti_species_file = muti_species_file
		self.muti_species_table = self.pandas.read_table(muti_species_file,names=["ncbi_taxid","alignment_file"])
		self.phylostratigraphy = phylostratigraphy
		self.multi_phy = self.muti_species_table.apply(lambda x: self.phylostratigraphy(ref_specie=x["ncbi_taxid"],alns=x["alignment_file"], n_bootstrap=self.n_bootstrap, h5_gi_taxid_file=self.h5_gi_taxid_file, fraction_rows_to_bootstrap = self.fraction_rows_to_bootstrap, n_rows_to_bootstrap = self.n_rows_to_bootstrap, ancestor_times = self.ancestor_times),axis=1).values

	def run(self):
		for taxa_ps in self.multi_phy:
			taxa_ps.pipeline()


