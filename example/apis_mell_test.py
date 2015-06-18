
import ComPhy as comphy

comphy.utils.psutils.store_gi_to_taxid_mapping_to_h5(gi2taxid='ComPhy/ncbi_tax_data/gi_taxid_prot.dmp.gz',store_file='ComPhy/dbs/gi_taxid.h5', close=True,override=False)

ps = comphy.phylostratigraphy(ref_specie=7460,alns="ComPhy/example/Apis_mellifera.Amel_2.0.13.pep.all.fa.nr.e10.sw.txt",h5_gi_taxid_file="ComPhy/dbs/gi_taxid.h5",n_bootstrap=10)
ps.pipeline()

print ps.phyloStrat_table_summary
print ps.phyloStrat_table_wide.head()

print ps.phyloStrat_table_wide_bootstrap.head()


