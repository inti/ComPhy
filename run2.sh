#/usr/bin/bash

# run test of script.

perl phylostratiphy.pl \
    -blast dysbindin.blast_out.txt \
    -tax_folder data/ \
    -query_taxon 9606 \
    -hard 1e-3 \
    -virus_list tmp.virus.txt \
    -out test_phylostratiphy


