#/usr/bin/bash

# run test of script.

my_perl phylostratiphy.pl \
    -blast dysbindin.blast_out.txt \ #blastp__Aech_v3.8.pep__nr.e10.sw.txt_old \
    -tax_folder data/ \
    -prot_only \
    -query_taxon 9606 \ #103372 \
    -virus_list tmp.virus.txt \
    -out test_phylostratiphy
