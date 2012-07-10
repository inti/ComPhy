#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::DB::Taxonomy;
# local modules
use PhyloStratiphytUtils;

my $tax_folder = $ARGV[0];

print_OUT("   '-> Reading phylogenetic tree and species information");
my $nodesfile = $tax_folder . "nodes.dmp";
my $namefile = $tax_folder . "names.dmp";
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -nodesfile => $nodesfile, -namesfile => $namefile);

my $taxon= $db->get_taxon(10239);
print_OUT("Starting to identify viral taxa");
my %all_viral_taxa = ();
my @taxons = $db->each_Descendent($taxon);
do {
	my @new_taxons = ();
	foreach my $taxa (@taxons) {
		push @new_taxons, $db->each_Descendent($taxa);
	}
	@taxons = ();
	@taxons = @new_taxons;
} until (scalar @taxons == 0);

print_OUT("Printing out");

open (OUT,">$ARGV[1]") or die $!;
foreach my $tx (keys %all_viral_taxa){
	print OUT $tx->id,"\t",$tx->scientific_name,"\n";
}
close(OUT);

print_OUT("Done");

exit;
