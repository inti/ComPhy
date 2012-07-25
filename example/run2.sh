#/usr/bin/bash

# run test of script.

perl phylostratiphy.pl \
    -blast dysbindin.blast_out.txt \
    -tax_folder data/ \
    -query_taxon 9606 \
	 -gi_tax_id dysbindin.tax_info.csv  \
    -hard 1e-3 \
    -virus_list virus_list.txt \
    -out test_phylostratiphy

