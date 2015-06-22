import pandas as pd
from sqlalchemy import create_engine 
import datetime as dt 
def set_ETE2_NCBI_db(update=False):
	from ete2 import NCBITaxa
	ncbi = NCBITaxa()
	if update:
		ncbi.update_taxonomy_database()
	return ncbi

def set_NCBI_GI_to_TAX_db(gi_taxid_dump,db_name,chunk_size):
	print dt.datetime.now(), "Openning file [",gi_taxid_dump,"]"
	file = pd.read_table(gi_taxid_dump,compression="gzip",chunksize=chunk_size,iterator=True,names=["gi","ncbi_taxid"])

	print dt.datetime.now(), "Starting to load data on database [",gi_taxid_dump,"] on chunks of [",chunk_size,"] lines"

	chunk = file.get_chunk()
	chunk.to_hdf(db_name,"gi_taxid",more="w",format="table", data_columns=["gi"],append=False,complib="blosc",index=False)

	for i,chunk in enumerate(file):
    		chunk.to_hdf(db_name,"gi_taxid",more="a",format="table",data_columns=["gi"],append=True,complib="blosc",index=False)

	print dt.datetime.now(), "   '-> done loading into [",gi_taxid_dump,"]"
	print dt.datetime.now(), "Indexing data"
	h5_store = pd.HDFStore(db_name)
	h5_store.create_table_index('gi_taxid',columns=["gi"], kind='full')
	h5_store.close()
	print dt.datetime.now(), "   '-> done "
	# Commented: input tables is already sorted by gi so no need to build a CSI
	#print dt.datetime.now(), "Building CSI index on gi column"
	#db_name_sorted = ''.join((db_name.strip("h5"),"sorted.h5"))
	#cmd = ''.join(("ptrepack --complib=blosc --chunkshape=auto --sortby=gi ",db_name," ", db_name_sorted)) #ptrepack --chunkshape=auto --sortby=gi gi_taxid_prot.h5 gi_taxid_prot.sorted.h5 --complib=blosc --complevel=6  --dont-regenerate-old-indexes --overwrite-nodes
	#os.sys(cmd)
        #print dt.datetime.now(), "   '-> done "


