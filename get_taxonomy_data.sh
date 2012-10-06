echo Starting to download NCBI taxonomy db info
out_folder=ncbi_tax_data

mkdir $out_folder

echo moving into folder [ $out_folder ]
cd $out_folder

tax_ftp="ftp://ftp.ncbi.nih.gov/pub/taxonomy"
echo Data from FTP folder [ $tax_ftp ]

for file in taxdump.tar.gz taxdump_readme.txt
do 
	wget $tax_ftp/$file
done
			
echo Uncompressing
ls *.gz | xargs gunzip -v
ls *.tar | perl -ne 'chomp($_); system "tar xvf $_";'

echo Making binary file for [ gi_taxid_prot.dmp ]
perl -e 'use Bio::LITE::Taxonomy::NCBI::Gi2taxid qw/new_dict/; new_dict (in => "gi_taxid_prot.dmp", out => "gi_taxid_prot.bin");'

echo Making binary file for [ gi_taxid_nucl.dmp ]
perl -e 'use Bio::LITE::Taxonomy::NCBI::Gi2taxid qw/new_dict/; new_dict (in => "gi_taxid_nucl.dmp", out => "gi_taxid_nucl.bin");'

rm citations.dmp delnodes.dmp division.dmp gc.prt gencode.dmp merged.dmp taxdump.tar

echo Files updated on [ `date` ] > files_update_info.txt

cd ..

echo Done  


