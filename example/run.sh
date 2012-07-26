#/usr/bin/bash

# run test of script.

my_perl ../phylostratiphy.pl \
    -blast example/blastp__Aech_v3.8.pep__nr.e10.sw.txt_old \
    -tax_folder ncbi_tax_data/ \
    -prot_only \
    -query_taxon 103372 \
    -out test_phylostratiphy


