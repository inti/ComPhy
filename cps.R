require(reshape2)
require(plyr);
require(stringr);

# load required functions
source("r_code/functions.R");
print_out("Loading local functions");
# read NCBI taxonomy data

tax_data = read_tax_data();

# read blast results
bp = read.table("example/blastp_Pbar_ant_vs_nr_e10.txt",sep="\t");
colnames(bp) = c("query_id", "subject_id", "perct_identity", "alignment_length", "query_length", "subject_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score")

# extrat gi id from blast hit identifier
bp$gi = str_extract(bp$subject_id, "(\\d+)");
# make data frame with set of hits
matches_gis = data.frame(gi = unique(bp$gi));
# extract tax information for those hits
matches_gis = merge(matches_gis,tax_data$gi_taxid_prot,by="gi")

# combine blast results with tax information
bp = merge(bp, matches_gis, by="gi");
# calculate hit score for each gene on each taxon
# i.e., ∑_1_J (w_ij * log(p_ij)), where the p_ij is the p-value from the blast result of the i genes with its j hit and w_ij is the coverage of the hit on the query.
gene_at_taxon_score = ddply(bp, .(query_id,tax_id), summarise, gene_at_taxon_score =  sum(  ( alignment_length/query_length)* -log(evalue),na.rm=T) ,.progress="text",.parallel=T);
# calcute total score for gene across all taxons
# ∑_1_T ∑_1_J (w_ij * log(p_ij)) with T taxons. 
total_score_gene_across_taxons = ddply(gene_at_taxon_score,.(query_id), summarise, total_score_per_gene = sum(gene_at_taxon_score,na.rm=T),.progress="text",.parallel=T);

# merge the two calculation results
gene_at_taxon_score = merge(gene_at_taxon_score,total_score_gene_across_taxons,by="query_id")

# normalize the per txon score to sum 1.
# equivalent to calculate ∑_1_J (w_ij * log(p_ij)) / [ ∑_1_T ∑_1_J (w_ij * log(p_ij)) ]
gene_at_taxon_score$norm_gene_at_taxon_score = gene_at_taxon_score$gene_at_taxon_score/gene_at_taxon_score$total_score_per_gene

# for each taxon calculate the genome-wide score
species_taxon_score = ddply(gene_at_taxon_score, .(tax_id), summarise ,species_taxon_score= sum(norm_gene_at_taxon_score,na.rm=T))
species_taxon_score = arrange(species_taxon_score,desc(species_taxon_score))

get_path_to_root(tax_data$nodes.dmp,9606,tax_data$names.dmp)

