DATE=28-12-2018
for group in Archaea Bacteria Eukaryota Viruses
do 
	echo `date` Working on $group
 	mkdir $group
	awk -F"\t" '{print $1,"_",$3}' uniprot_ref_prot_${group}_${DATE}.tab|  sed "s/ //g" | xargs -P5 -I{} echo "wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/${group}/{}.fasta.gz &>{}.log" > ${group}.cmd 
	echo `date` Finished making download commands
	sleep 5
	echo `date` Pulling files down 
	wc -l ${group}.cmd
	cd $group
	~/anaconda2/envs/ngs/bin/parallel -j 5 --no-notice -a ../${group}.cmd --progress {}   
	cd ..
done
