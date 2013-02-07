
### define some functions
__name__ = 'PhyloStratiphytUtils'

__version__ = '0.01'

import datetime
import pandas as pd
import numpy as np
import os.path
import TreeTaxonomyUtils as ttutils
import progressbar 
from tables import *
from Bio import SeqIO
from cogent.parse.ncbi_taxonomy import NcbiTaxonomyFromFiles

def add_taxid_from_h5file(table,file):
    print return_time(), "Reading HDF5 file with GI to TaxIds file [",file,"]"
    h5file = pd.HDFStore(file, mode = "r")
    table_grouped = table.groupby('subject_gi')
    df_unique_gi = pd.DataFrame(table_grouped.groups.keys(),columns = ['subject_gi'])
    n_unique_gi = table_grouped.ngroups
    print return_time(), "   ... [",n_unique_gi,"] GIs to search for"
    b_chunk_size = 1000000
    nrows_b = h5file.root.gi_2_taxid.table.nrows
    df_lot = []
    pbar = progressbar.ProgressBar(maxval=nrows_b+1, term_width=50,widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()]).start()
    for b in xrange(int(nrows_b / b_chunk_size) + 1):
        b_start_i = b * b_chunk_size
        b_stop_i = min((b + 1) * b_chunk_size, nrows_b)
        
        B = h5file.select('gi_2_taxid', start = b_start_i, stop = b_stop_i)
        intersection = df_unique_gi.ix[ df_unique_gi['subject_gi'].isin( B['gi'] ),'subject_gi']
        if len(intersection) > 0:
            df_lot.append(pd.merge(df_unique_gi,B,left_on='subject_gi',right_on='gi',how='inner'))
        pbar.update(b_stop_i-1)
    pbar.finish()
    df_unique_gi = pd.concat(df_lot)
    print return_time(), "   ... found TaxIds for [",len(df_unique_gi),"] GIs"
    print return_time(), "   ... done"
    return df_unique_gi

def set_taxid(gi,taxid):
    table.ix[table_grouped.groups[gi],'taxid']=taxid

def add_taxid_from_h5file_old(table,file):
    filename = file
    print return_time(), "Reading HDF5 file with GI to TaxIds file [",file,"]"
    h5file = pd.HDFStore(file, mode = "r")
    unique_gi = pd.unique(table['subject_gi'])
    n_unique_gi = len(unique_gi)
    print return_time(), "   ... [",n_unique_gi,"] GIs to search for"
    pbar = progressbar.ProgressBar(maxval=n_unique_gi, term_width=50,widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()]).start()
    i = 0
    for sbj_gi in unique_gi:
        i = i+1
        pbar.update(i)
        taxid = h5file.root.gi_2_taxid.table.readWhere('(gi == sbj_gi)')
        if len(taxid) > 0:
            taxid = taxid[0][2]
            table.ix[table['subject_gi'] == sbj_gi,'taxid' ] = taxid
    pbar.finish()

    print return_time(), "   ... found TaxIds for [",len(pd.unique(table.ix[table['taxid'] > 0,'taxid'])),"] GIs"
    print return_time(), "   ... assigned TaxIds to [",len(table.ix[table['taxid'] > 0])," out of [",len(table),"] sequence search results"
    print return_time(), "   ... done"


def make_gi_tp_taxonomy_hd5_file(gi2taxid,store_file):
    print return_time(), "Reading GIs to TaxIds file and storing data onto [",store_file,"]"
    hd_store =  pd.HDFStore(store_file)
    gi_2_taxid = pd.read_table(gi2taxid, chunksize=1000000,compression='gzip',header=None,sep="\t",names=['gi','taxid'])
    for chunk in gi_2_taxid:
        hd_store.append('gi_2_taxid', chunk,data_columns=True)
    print "   ... done"
    return hd_store

def store_gi_to_taxid_mapping_to_h5(gi2taxid, store_file='gi_taxid.h5', close=False,override=False):
    if os.path.exists(store_file):
        print return_time(), "HDF5 file exists"
        if override == False:
            print return_time(), "   ... reading HDF5 file [",store_file,"]"
            hd_store = pd.HDFStore(store_file)
        else:
            print return_time(), "   ... will override file [",store_file,"]"
            hd_store = make_gi_tp_taxonomy_hd5_file(gi2taxid,store_file)

    else:
        hd_store = make_gi_tp_taxonomy_hd5_file(gi2taxid,store_file)

    if close == True:
        print return_time(), "   ... done"
        hd_store.close()
    else:
        print return_time(), "   ... done"
        return hd_store

def add_tax_id(data,h5_tbl):
    print return_time(), "Fecthing TaxIds from HDF5 file"
    data['taxid'] = data['subject_gi'].apply(lambda x: h5_tbl.readWhere('(gi == x)')['taxid'])
    empty_cells = data['taxid'].apply(len) == 0
    data.ix[empty_cells,'taxid'] = np.nan
    data['taxid'] = np.hstack(data['taxid'])
    print return_time(), "   ... done"


# map sequence search query sequences to the oldest LCA with eny of the hits passing the e-value threshold
def map_to_oldest_stratum(data,evalue_threshold=1e-3):
    data = data.ix[data['evalue'] <= evalue_threshold]
    grouped = data.groupby('query_id')
    back = grouped.apply(get_oldest_hit)
    return back

# rounting to get the oldest of the hits from the table. It is based on the rank of the phylostratum and assumes that lowest rank means oldest
def get_oldest_hit(data):
    where_lowest = data['lca_rank'] == min(data['lca_rank'])
    back = list( set(data.ix[ where_lowest  ,'lca_name']) )[0]
    return back

# returns a summary wit the number of genes mapped to each phylostratum
def count_genes_per_phylostrata(data,nodes_table, names_table, ref_species_id=9606):
    print return_time(), "Counting genes assign to each phylostrata"
    print return_time(), "   ... finding reference specie taxonomy"
    ref_lineage = ttutils.find_lineage_ById( ref_species_id, nodes_table)
    ref_lineage_names = ttutils.get_names_for_list_ById(ref_lineage,names_table)
    print return_time(), "   ... counting genes mapped to each ancestor of ref specie"
    counts = pd.DataFrame(data.values,columns=['phylostratum'],index = data.index).groupby('phylostratum').size()
    counts = counts.reindex(ref_lineage_names,fill_value=0)
    ranks = range(len(ref_lineage_names) - 1,-1,-1)
    counts = pd.DataFrame({'phylostratum':counts,'rank':ranks},index = ref_lineage_names)
    print return_time(), "   ... done"
    return counts


# Load species tree from NCBI Taxonomy
def load_ncbi_taxonomy_tree(nodes,names):
    print return_time(), "Reading NBCI Taxonomy files"
    tree = NcbiTaxonomyFromFiles( open(nodes) , open(names))
    print return_time(), "   ... done"
    return tree


#load gi to tax id mapping
def load_gi_taxid_mapping(file,table,chunk_size=1000000,compression = 'gzip'):
    print return_time(), "Reading sequence GI to Taxonomy Id mapping file [", file,"]"
    gi_2_taxid = pd.read_table(file, chunksize=chunk_size,compression=compression,header=None,sep="\t",names=['gi','taxid'], index_col=0)
    print return_time(), "   '-> Matching blast results GIs with sequences Taxonomy Ids and filtering blast result table"
    table['taxid'] = np.nan
    for chunk in gi_2_taxid:
        # get array with subject_id that overlap between two tables
        tmp_table = table.ix[chunk.index,'subject_gi'].values
        if len(tmp_table) > 0:
            table.ix[tmp_table,'taxid'] = chunk.ix[tmp_table,'taxid']
    table_size = len(table)
    without_matches = sum(np.isnan(table['taxid']))
    print return_time(), "     ... [", table_size - without_matches, "] of [",table_size,"] sequences had identified TaxIds"
    print return_time(), "   ... done"


# Read sequences used for the analyses
def read_fasta_seqs(file):
    print return_time(), "Reading sequences used for blast search"
    sequence_length = {}
    handle = open("example/dysbindin.fa", "rU")
    for record in SeqIO.parse(handle, "fasta") :
        sequence_length[record.id] = len(record.seq)
    handle.close()
    return sequence_length


def filter_table_by_ids(table, ids):
    print return_time(), "Filtering table"
    start_size = len(table)
    table = table[table['query_id'].isin(ids) ]
    end_size = len(table)
    print return_time(), "   '-> Filtered out [", start_size - end_size, "] rows of the table"
    print return_time(), "   ... done"
    return table

# function to load blast table
def load_blast_results(file,format='blastp_table',evalue_limit=1e-3):
    print return_time(), "Reading sequence search results from [", file, "]"
    if format == 'blastp_table':
        colnames = ['query_id', 'subject_id', 'identity', 'alignment_length', 'query_length', 'subject_length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']
        table = pd.read_table(file, "r", delimiter="\t", header=None, names = colnames)
        print return_time(), "   '-> [", len(table), "] results read"
        table = table.ix[table['evalue']<= evalue_limit]
        print return_time(), "        ... [", len(table), "] with e-values <= [",evalue_limit,"]"
        # reduce the blast output only to the sequences present on the fasta file
        #table = table[table['query_id'].isin(sequence_length.keys()) ]
        table['subject_gi'] = table['subject_id'].str.split('|').apply(lambda x: x[1])
        table['subject_gi'] = table['subject_gi'].astype(np.int64)
        #table.index = table['subject_gi']
        print return_time(), "   ... done"
        return table


# function to print out a string
def print_OUT(x):
    print return_time(), x


## simple function to return time for printing
def return_time():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

# class and function to read a text file with commented lines
class FileWrapper(file):
    def __init__(self, comment_literal, *args):
        super(FileWrapper, self).__init__(*args)
        self._comment_literal = comment_literal

    def next(self):
        while True:
            line = super(FileWrapper, self).next()
            if not line.startswith(self._comment_literal):
                return line





