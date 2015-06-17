
import ComPhy as comphy
from ComPhy.PhyloStratiphy import phylostatigraphy

ps = comphy.phylostratigraphy(ref_specie=7460,alns="example/Apis_mellifera.Amel_2.0.13.pep.all.fa.nr.e10.sw.txt",n_bootstrap=10)
ps.pipeline()

print ps.phyloStrat_table_summary
print ps.phyloStrat_table_wide.head()

print ps.phyloStrat_table_wide_bootstrap.head()


