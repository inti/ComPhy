### define some functions
__name__ = 'TreeTaxonomyUtils'

__version__ = '0.01'

import pandas as pd
import numpy as np
import PhyloStratiphytUtils as psutils
import progressbar

def my_strip(text):
    try:
         return text.strip(' \t')
    except AttributeError:
        return text

# wrapper over funtions to read files
def load_ncbi_tax_files(nodes_file,names_file):
    print psutils.return_time(), "Reading NCBI Taxonomy file"
    print psutils.return_time(), "   ... [", nodes_file,"]"
    nodes_table = read_nodes_file(nodes_file)
    print psutils.return_time(), "   ... [", names_file,"]"
    names_table = read_names_file(names_file)
    print psutils.return_time(), "   ... done"
    return [ nodes_table, names_table ]


#read names file
def read_names_file(file):
    tax_names = pd.read_table(file,
                          header=None,sep="\|",
                          names=['taxid','name_txt','unique_name','name_class','extra'],
                          converters= {'taxid': my_strip,
                                        'name_txt':my_strip,'unique_name':my_strip,'name_class':my_strip})
    tax_names = tax_names.drop(['extra'],axis=1)
    tax_names['taxid'] = tax_names['taxid'].astype(np.int64)
    return tax_names

# read nodes file
def read_nodes_file(file):
    colnames = ['taxid','parent_tax_id','rank','embl_code','division_id','inherited_div_flag','genetic_code_id','inherited_GC','mitochondrial_genetic_code_id','inherited_MGC','GenBank_hidden','hidden_subtree_root','comments','extra']
    tax_nodes = pd.read_table(file, header=None,sep="\|", names=['taxid','parent_tax_id','rank','embl_code','division_id','inherited_div_flag','genetic_code_id','inherited_GC','mitochondrial_genetic_code_id','inherited_MGC','GenBank_hidden','hidden_subtree_root','comments','extra'], converters= {'taxid': my_strip,'parent_tax_id': my_strip,'rank': my_strip,'embl_code': my_strip,'division_id': my_strip,'inherited_div_flag': my_strip,'genetic_code_id': my_strip,'inherited_GC': my_strip,'mitochondrial_genetic_code_id': my_strip,'inherited_MGC': my_strip,'GenBank_hidden': my_strip,'hidden_subtree_root': my_strip,'comments': my_strip })
    tax_nodes = tax_nodes.drop(['extra'],axis=1)
    tax_nodes['taxid'] = tax_nodes['taxid'].astype(np.int64)
    tax_nodes['parent_tax_id'] = tax_nodes['parent_tax_id'].astype(np.int64)
    return tax_nodes


def find_lineage_ById(taxid,nodes_table):
    lineage = []
    taxid = np.int64(taxid)
    lineage.append(taxid)
    parent_id = taxid
    while parent_id > 1:
        parent_id = nodes_table.ix[ nodes_table['taxid'] == parent_id ,"parent_tax_id"]
#        parent_id = np.int64(parent_id.values[0])
        parent_id = np.int64(parent_id)
        lineage.append(parent_id)
    return lineage

##### OLD ROUTINE
# from a tree node object it return the lineage of the node.
# it assumes the species tree is stored on 'tree' and available on the working enviroment
#def get_lineage(node,tree,value = 'name'):
#    """returns lineage information in my_ranks order"""
#    curr = node
#    if value == 'name':
#        lineage = [curr.Name]
#    else:
#        lineage = [curr.TaxonId]
#    while curr.Parent is not None:
#        curr = curr.Parent
#        if value == 'name':
#            lineage.append(curr.Name)
#        else:
#            lineage.append(curr.TaxonId)
#        parent = tree.ById[curr.ParentId]
#    return lineage

def get_name_ById(taxid,names_table):
    name = names_table.ix[(names_table['taxid'] == taxid) & (names_table['name_class'] == "scientific name"),"name_txt"]
    name = name.values[0]
    return name

def get_names_for_list_ById(id_list, names_table):
    all_names = [ get_name_ById(id,names_table) for id in id_list]
    return all_names

# provide two lineages and return their last common ancestor

def get_lca_for_list(coming_in,nodes_table,names_table,ref_species_id=9606):
    print psutils.return_time(), "Identifying last common ancestor between query specie and blast hits"
    ref_lineage = find_lineage_ById( ref_species_id, nodes_table)
    ref_lineage = {ref_lineage[i]: i   for i in xrange(len(ref_lineage))}
    lca_stack = {}
    inner_nodes_lca_stack = {}
    for ancestor in ref_lineage:
        inner_nodes_lca_stack[ancestor] = ancestor
    
    pbar = progressbar.ProgressBar(maxval=len(coming_in)+1, term_width=50,widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()]).start()
    i=0
    for sbj_taxid in coming_in:
        pbar.update(i)
        i = i + 1
        if np.isnan(sbj_taxid):
            continue
    
        if ref_species_id == sbj_taxid:
            subject_name = get_name_ById(sbj_taxid,names_table)
            lca_stack[sbj_taxid] = { 'subject_name' : subject_name , 'subject_taxid': sbj_taxid, 'lca_name': subject_name, 'lca_tax_id': sbj_taxid, 'lca_rank' : len(ref_lineage) }
            continue

        parent_id = np.int64(sbj_taxid)
        lineage = []
        lineage.append(sbj_taxid)
        while parent_id > 1:
            parent_id = nodes_table.ix[ nodes_table['taxid'] == parent_id ,"parent_tax_id"]
            parent_id = np.int64(parent_id)
            if not parent_id:
                break

            lineage.append(parent_id)
            if parent_id == 10239: # if this sequence is a virus
                for viral_ancestor in lineage:
                    inner_nodes_lca_stack[viral_ancestor] = 'virus'
                break
    
            if parent_id in ref_lineage.keys():
                subject_name = get_name_ById(sbj_taxid,names_table)
                lca_name = get_name_ById(parent_id,names_table)
                lca_rank = ref_lineage[parent_id]
                lca_stack[sbj_taxid] = { 'subject_name' : subject_name , 'subject_taxid': sbj_taxid, 'lca_name': lca_name, 'lca_tax_id': parent_id, 'lca_rank' : lca_rank }
                for ancestor in lineage:
                    inner_nodes_lca_stack[ancestor] = parent_id
                break

            if parent_id in inner_nodes_lca_stack.keys():
                lca_tax_id = inner_nodes_lca_stack[ parent_id ]
                if lca_tax_id == 'virus':
                    for viral_ancestor in lineage:
                        inner_nodes_lca_stack[ viral_ancestor ] = 'virus'
                    break
                    
                subject_name = get_name_ById(sbj_taxid,names_table)
                lca_name = get_name_ById(lca_tax_id,names_table)
                lca_stack[sbj_taxid] = { 'subject_name' : subject_name , 'subject_taxid': sbj_taxid, 'lca_name': lca_name, 'lca_tax_id': lca_tax_id, 'lca_rank' : len(ref_lineage) - len(lineage) }
                for ancestor in lineage:
                    inner_nodes_lca_stack[ancestor] = lca_tax_id
                break

    pbar.finish()
    lca_stack = pd.DataFrame(lca_stack).T
    print psutils.return_time(), "   ... done"
    return lca_stack

def get_lca(lineage1,lineage2):
    for i in range(1,1+len(lineage1)):
        if lineage1[-i] != lineage2[-i]:
            lca = lineage1[-(i - 1)]
            return [lca,i-1]
            break

# for a list of ids it returns the last common ancestor (and additional information) between the taxid on the list and a refid.
def get_lca_for_list_old(coming_in,nodes_table,names_table,ref_species_id=9606):
    print psutils.return_time(), "Identifying last common ancestor between query specie and blast hits"
    back = {}
    ref_lineage = find_lineage_ById( ref_species_id, nodes_table)
    ref_species_name = get_name_ById( ref_species_id, names_table)
    for sbj_taxid in coming_in:
        if np.isnan(sbj_taxid):
            continue
        else:
            if sbj_taxid == ref_species_id:
                back[sbj_taxid] = [ref_species_name, ref_species_id, ref_species_name, ref_species_id,len(ref_lineage)]
                continue
            
            sbj_name = get_name_ById(sbj_taxid,names_table)
            sbj_node_lineage = find_lineage_ById(sbj_taxid,nodes_table)
            #print ref_node.Name," :: ", sbj_node.Name, "\n",sbj_node_lineage,"\n",ref_lineage
            qry_sbj_lca_id, qry_sbj_lca_rank = get_lca(ref_lineage,sbj_node_lineage)
            qry_sbj_lca_name = get_name_ById(qry_sbj_lca_id,names_table)
            back[sbj_taxid] = [sbj_name, sbj_taxid, qry_sbj_lca_name, qry_sbj_lca_id,qry_sbj_lca_rank]

    back = pd.DataFrame(back.values(), columns=['subject_name','subject_taxid','lca_name','lca_tax_id','lca_rank'])
    print psutils.return_time(), "   ... done"
    return back



