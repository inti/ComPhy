#/usr/bin/bash

# run test of script.

my_perl phylostratiphy.pl \
    -blast dysbindin.blast_out.txt 
    -tax_folder data/ \
    -prot_only \
    -query_taxon 9606 \
    -virus_list tmp.virus.txt \
    -out test_phylostratiphy


