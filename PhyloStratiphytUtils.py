
### define some functions
__name__ = 'PhyloStratiphytUtils'

__version__ = '0.01'

import datetime
import pandas as pd
from numpy import isnan

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
def get_lca_for_list(coming_in,tree):
    back = {}
    ref_node = tree.ById[9606]
    ref_lineage = get_lineage(ref_node,tree,'id')
    for number in coming_in:
        if isnan(number):
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
    return back
