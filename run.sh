#/usr/bin/bash

# run test of script.

my_perl phylostratiphy.pl \
    -blast blastp__Aech_v3.8.pep__nr.e10.sw.txt_old \
    -tax_folder data2/ \
    -prot_only \
    -query_taxon 103372 \
    -virus_list tmp.virus.txt \
    -out test_phylostratiphy

