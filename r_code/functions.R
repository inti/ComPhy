print_out<-function(string,LOG){
  print(paste(date(),"      ",string),quote = FALSE)
}

get_path_to_root = function(nodes_info, tax_id,names_info) {
  parent_id = 9999999; 
  id = tax_id;
  path = data.frame(tax_id = id,name_text = nodes_info[which(nodes_info$tax_id == id ),"name_txt"]);
  c <- 2;
  while ( parent_id > 1  ){
    path[counter,c("tax_id","name_text")]parent_id = nodes_info[which(nodes_info$tax_id == id ),c("parent_tax_id","name_txt")];
    c <- c +1;
    id = parent_id;
    #cat(str_c(id,names_info[which(names_info$tax_id == id || names_info$`name class` == "scientific name"),"name_txt"],sep=" "),"\n");
  }
  return(path);
}
read_tax_data_and_compress = function () {
        print_out("Reading gi to tax id mapping");
        print_out("This may take a few minutes");
        #gi_taxid_nucl = read.table("gi_taxid_nucl.dmp",header=F,sep="\t");
        gi_taxid_prot = read.table("gi_taxid_prot.dmp",header=F,sep="\t");
        names(gi_taxid_prot) = c("gi","tax_id");
        names.dmp = read.table("names.dmp_v2",sep="|",fill=T);
        colnames(names.dmp) = c("tax_id","name_txt","unique_name","name class");
        nodes.dmp = read.table("nodes.dmp_v2",sep="|",fill=T);
        colnames(nodes.dmp) = c("tax_id","parent_tax_id","rank","embl_code","division_id","inherited_div_flag","genetic_code_id","inherited_GC_flag","mitochondrial_genetic_code_id","inherited_MGC_flag","GenBank_hidden_flag","hidden_subtree_root_flag","comments");
	print_out("compressing files");
	save(gi_taxid_prot, file = "gi_taxid_prot.RData", compress="xz");
        save(names.dmp, file = "names.RData", compress="xz");
        save(nodes.dmp, file = "nodes.RData", compress="xz");
        print_out("done reading and compressing tax information");
        return(0);
}

read_tax_data = function () {
        print_out("Reading gi to tax id mapping");
	print_out("This may take a few minutes"); 
	load("data/gi_taxid_prot.RData");
	load("data/names.RData");
	load("data/nodes.RData");
	print_out("done reading tax information");
	return(list(nodes = nodes.dmp, names = names.dmp, gi_taxid_prot = gi_taxid_prot));
}
