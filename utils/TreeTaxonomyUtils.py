### define some functions
__name__ = 'TreeTaxonomyUtils'

__version__ = '0.01'

import pandas as pd
import numpy as np
#import PhyloStratiphytUtils as psutils
import datetime 

def return_time():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

def my_strip(text):
    try:
         return text.strip(' \t')
    except AttributeError:
        return text

# wrapper over funtions to read files
def load_ncbi_tax_files(nodes_file,names_file):
    print return_time(), "Reading NCBI Taxonomy file"
    print return_time(), "   ... [", nodes_file,"]"
    nodes_table = read_nodes_file(nodes_file)
    print return_time(), "   ... [", names_file,"]"
    names_table = read_names_file(names_file)
    print return_time(), "   ... done"
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
        parent_id = np.int64(parent_id.values[0])
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
def get_lca(lineage1,lineage2):
    for i in range(1,1+len(lineage1)):
        if lineage1[-i] != lineage2[-i]:
            lca = lineage1[-(i - 1)]
            return [lca,i-1]
            break


# for a list of ids it returns the last common ancestor (and additional information) between the taxid on the list and a refid.
def get_lca_for_list(coming_in,nodes_table,names_table,ref_species_id=9606):
    print return_time(), "Identifying last common ancestor between query specie and blast hits"
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
    print return_time(), "   ... done"
    return back

def get_lca_with_root2tip_lineages(lineage1,lineage2,ncbi):
    # cap_size store the min length between the two lineages since the LCA must be between the root and this rank
    size_cap = np.minimum(len(lineage1),len(lineage2))
    # if they do not share the root of the tree there cannot be a LCA
    if lineage1[0] != lineage2[0]: return None, None, None
    # if tip is the same then it is the same specie
    if lineage1[-1] == lineage2[-1]: return lineage1[-1], len(lineage1), ncbi.get_taxid_translator([lineage1[-1]])[lineage1[-1]] 
    # it could happend that a sequence is annotated to a subspecie of the specie being analysed. On this case we will map the LCA to the specie, since by definition 
    # the subspecie is not a differente specie and it is on the same *lineage*.
    if lineage1[size_cap-1] == lineage2[size_cap-1]: return lineage1[size_cap-1], len(lineage1), ncbi.get_taxid_translator([lineage1[size_cap-1]])[lineage1[size_cap-1]]    

    lca_rank = np.where(np.subtract(lineage1[:size_cap],lineage2[:size_cap]) != 0)[0].min() - 1
    lca_taxid = lineage1[lca_rank]
    lca_common_name = ncbi.get_taxid_translator([lca_taxid])[lca_taxid]
    return lca_taxid, lca_rank, lca_common_name


def _get_commonName_and_lineage(taxid,ncbi):
    return ncbi.get_taxid_translator([taxid])[taxid], ncbi.get_lineage( taxid )

# for a list of ids it returns the last common ancestor (and additional information) between the taxid on the list and a refid.
def get_lca_for_list_with_ETE2(coming_in,ref_species_id,ncbi):
    from ete2 import NCBITaxa
    ncbi = NCBITaxa()
    print return_time(), "Identifying last common ancestor between query specie and blast hits"
    back = {}
    
    ref_species_name, ref_lineage = _get_commonName_and_lineage(ref_species_id,ncbi)
    for sbj_taxid in coming_in:
        if np.isnan(sbj_taxid):
            continue
        else:
            sbj_name, sbj_node_lineage = _get_commonName_and_lineage( sbj_taxid , ncbi)
	    # some times there are erros on the taxid and one needs to find the common name a the real stable taxid from the name
	    if len(sbj_node_lineage) < 2:
		sbj_new_taxid = ncbi.get_name_translator([sbj_name])[sbj_name]
		sbj_name, sbj_node_lineage = _get_commonName_and_lineage( sbj_new_taxid , ncbi)
		if len(sbj_node_lineage) < 2:
			back[sbj_taxid] = [sbj_name, sbj_taxid, None, None, None]
			break
            lca_taxid, lca_rank, lca_common_name = get_lca_with_root2tip_lineages(ref_lineage,sbj_node_lineage,ncbi)
	    back[sbj_taxid] = [sbj_name, sbj_taxid, lca_common_name, lca_taxid,lca_rank]

    back = pd.DataFrame(back.values(), columns=['subject_name','subject_taxid','lca_name','lca_tax_id','lca_rank'])
    print return_time(), "   ... done"
    return back


