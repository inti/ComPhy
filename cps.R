library(reshape2)
library(plyr);
library(stringr);

# read NCBI taxonomy data
#gi_taxid_nucl = read.table("data/gi_taxid_nucl.dmp",header=F,sep="\t");
gi_taxid_prot = read.table("data/gi_taxid_prot.dmp",header=F,sep="\t");
names(gi_taxid_prot) = c("gi","tax_id");

names.dmp = read.table("data/names.dmp_v2",sep="|",fill=T);
colnames(names.dmp) = c("tax_id","name_txt","unique_name","name class");

nodes.dmp = read.table("data/nodes.dmp_v2",sep="|",fill=T);
colnames(nodes.dmp) = c("tax_id","parent_tax_id","rank","embl_code","division_id","inherited_div_flag","genetic_code_id","inherited_GC_flag","mitochondrial_genetic_code_id","inherited_MGC_flag","GenBank_hidden_flag","hidden_subtree_root_flag","comments");

bp = read.table("../../ComparaGenomics/PhyloStratiphy/Ants/blastp__pbar.genome.OGS.1.2.maker.proteins.fasta__nr.e10.sw.txt",sep="\t")
colnames(bp) = c("query_id", "subject_id", "perct_identity", "alignment_length", "query_length", "subject_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score")


bp$gi = str_extract(bp$subject_id, "(\\d+)");
matches_gis = data.frame(gi = unique(bp$gi));
matches_gis = merge(matches_gis,gi_taxid_prot,by="gi")

bp = merge(bp, matches_gis, by="gi")
gene_at_taxon_score = ddply(bp, .(query_id,tax_id), summarise, gene_at_taxon_score =  sum(  ( alignment_length/query_length)* -log(evalue)) ,.progress="text",.parallel=T)
total_score_gene_across_taxons = ddply(gene_at_taxon_score,.(query_id), summarise, total_score_per_gene = sum(gene_at_taxon_score))
gene_at_taxon_score$norm_gene_at_taxon_score = gene_at_taxon_score$gene_at_taxon_score/gene_at_taxon_score$total_score_per_gene

species_taxon_score = ddply(gene_at_taxon_score, .(tax_id), summarise = sum(norm_gene_at_taxon_score,na.rm=T))
species_taxon_score = species_taxon_score[order(-species_taxon_score$norm_taxon_score),]


get_path_to_root = function(nodes_info, tax_id,names_info) {
	parent_id = 2; 
	id = tax_id;
	while ( parent_id > 1  ){
		parent_id = nodes_info[which(nodes_info$tax_id == id ),"parent_tax_id"];
		id = parent_id;
		cat(str_c(id,names_info[which(names_info$tax_id == id || names_info$`name class` == "scientific name"),"name_txt"],sep=" "),"\n");
	}
}

get_path_to_root(nodes.dmp,9606,names.dmp)

