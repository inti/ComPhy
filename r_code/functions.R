print_OUT<-function(string,LOG){
  print(paste(date(),"      ",string),quote = FALSE)
}

get_path_to_root = function(nodes_info, tax_id,names_info) {
  parent_id = 2; 
  id = tax_id;
  while ( parent_id > 1  ){
    parent_id = nodes_info[which(nodes_info$tax_id == id ),"parent_tax_id"];
    id = parent_id;
    cat(str_c(id,names_info[which(names_info$tax_id == id || names_info$`name class` == "scientific name"),"name_txt"],sep=" "),"\n");
  }
}