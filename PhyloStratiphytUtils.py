
### define some functions
__name__ = 'PhyloStratiphytUtils'

__version__ = '0.01'

import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO
from cogent.parse.ncbi_taxonomy import NcbiTaxonomyFromFiles


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
def count_genes_per_phylostrata(data):
    counts = pd.DataFrame(data.values,columns=['phylostratum'],index = data.index).groupby('phylostratum').size()
    return counts


# Load species tree from NCBI Taxonomy
def load_ncbi_taxonomy_tree(nodes,names):
    print return_time(), "Reading NBCI Taxonomy files"
    tree = NcbiTaxonomyFromFiles( open(nodes) , open(names))
    print return_time(), "   ... done"
    return tree


#load gi to tax id mapping
def load_gi_taxid_mapping(file,table,chunk_size=10000,compression = 'gzip',filter='true'):
    start_size = len(table)
    print return_time(), "Reading sequence GI to Taxonomy Id mapping file [", file,"]"
    gi_2_taxid = pd.read_table(file, chunksize=chunk_size,compression=compression,header=None,sep="\t",names=['gi','taxid'])
    if filter == 'true':
        how_filter = 'inner'
    else:
        how_filter = 'outer'

    print return_time(), "   '-> Matching blast results GIs with sequences Taxonomy Ids and filtering blast result table"
    tmp_table = pd.concat([pd.merge(table,chunk,left_on='subject_gi', right_on='gi',how=how_filter) for chunk in gi_2_taxid], ignore_index=True)
    if filter == 'true':
            end_size = len(tmp_table)
            if len(tmp_table) == 0:
                print return_time(), "Something went wrong. I could map any sequence search result's hit to TaxIds. Check data type of columns with GIs. Returning original table"
                return table
            print return_time(), "   '-> Filtered out [", start_size - end_size, "] rows of the table"
            print return_time(), "   '-> [", end_size, "] blast results had mapable TaxIds"

    print return_time(), "   ... done"
    return tmp_table

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
def load_blast_results(file,format='blastp_table'):
    print return_time(), "Reading sequence search results from [", file, "]"
    if format == 'blastp_table':
        colnames = ['query_id', 'subject_id', 'identity', 'alignment_length', 'query_length', 'subject_length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']
        table = pd.read_table(file, "r", delimiter="\t", header=None, names = colnames)
        print return_time(), "   '-> [", len(table), "] results read"
        # reduce the blast output only to the sequences present on the fasta file
        #table = table[table['query_id'].isin(sequence_length.keys()) ]
        table['subject_gi'] = table['subject_id'].str.split('|').apply(lambda x: x[1])
        table['subject_gi']=table['subject_gi'].astype(np.int64)
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

# from a tree node object it return the lineage of the node.
# it assumes the species tree is stored on 'tree' and available on the working enviroment
def get_lineage(node,tree,value = 'name'):
    """returns lineage information in my_ranks order"""
    curr = node
    if value == 'name':
        lineage = [curr.Name]
    else:
        lineage = [curr.TaxonId]
    while curr.Parent is not None:
        curr = curr.Parent
        if value == 'name':
            lineage.append(curr.Name)
        else:
            lineage.append(curr.TaxonId)
        parent = tree.ById[curr.ParentId]
    return lineage

# provide two lineages and return their last common ancestor
def get_lca(lineage1,lineage2):
    for i in range(1,len(lineage1)):
        if lineage1[-i] != lineage2[-i]:
            lca = lineage1[-(i - 1)]
            return [lca,i-1]
            break

# for a list of ids it returns the last common ancestor (and additional information) between the taxid on the list and a refid.
def get_lca_for_list(coming_in,tree,ref_species_id=9606):
    print return_time(), "Identifying last common ancestor between query specie and blast hits"
    back = {}
    ref_node = tree.ById[ref_species_id]
    ref_lineage = get_lineage(ref_node,tree,'id')
    for number in coming_in:
        if np.isnan(number):
            continue
        else:
            sbj_node = tree.ById[number]
            if sbj_node.TaxonId == ref_node.TaxonId:
                back[number] = [sbj_node.Name, number, sbj_node.Name, sbj_node.TaxonId,len(ref_lineage)]
                continue

            sbj_node_lineage = get_lineage(sbj_node,tree,'id')
            qry_sbj_lca_id, qry_sbj_lca_rank = get_lca(ref_lineage,sbj_node_lineage)
            qry_sbj_lca_node = tree.ById[qry_sbj_lca_id]
            #print ref_node.Name," :: ", sbj_node.Name, " => ", qry_sbj_lca_node.Name," <=\n",sbj_node_lineage,"\n",ref_lineage
            back[number] = [sbj_node.Name, number, qry_sbj_lca_node.Name, qry_sbj_lca_id,qry_sbj_lca_rank]
    back = pd.DataFrame(back.values(), columns=['subject_name','subject_taxid','lca_name','lca_tax_id','lca_rank'])
    print return_time(), "   ... done"
    return back



