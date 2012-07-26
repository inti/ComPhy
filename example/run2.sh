#/usr/bin/bash

# run test of script.

cd ..

perl phylostratiphy.pl \
    -blast example/dysbindin.blast_out.txt \
    -tax_folder ncbi_tax_data/ \
    -query_taxon 9606 \
	 -gi_tax_id example/dysbindin.tax_info.csv  \
    -hard 1e-3 \
    -virus_list example/virus_list.txt \
    -out example/test_phylostratiphy

cd example

