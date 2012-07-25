echo Starting to download NCBI taxonomy db info
out_folder=data

mkdir $out_folder

echo moving into folder [ $out_folder ]
cd $out_folder

tax_ftp="ftp://ftp.ncbi.nih.gov/pub/taxonomy"
echo Data from FTP folder [ $tax_ftp ]

for file in gi_taxid_prot.dmp.gz gi_taxid_nucl.dmp.gz taxcat.tar.gz taxdump.tar.gz taxdump_readme.txt gi_taxid.readme taxcat_readme.txt 
do 
	wget $tax_ftp/$file
done

echo Uncompressing
ls *.gz | xargs gunzip -v
ls *.tar | perl -ne 'chomp($_); system "tar xvf $_";'

echo Done  

