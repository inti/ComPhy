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
			
rm citations.dmp delnodes.dmp division.dmp gc.prt gencode.dmp merged.dmp taxdump.tar

echo Files updated on [ `date` ] > files_update_info.txt

cd ..

echo Done  


