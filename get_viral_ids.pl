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

my %all_viral_taxa = ();
my @taxons = $db->each_Descendent($taxon);
do {
	print scalar localtime,"\t","Descendent [ ", scalar @taxons," ] and total [ ",scalar (keys %all_viral_taxa ) ," ]\n";
	my @new_taxons = ();
	foreach my $taxa (@taxons) {
		print $taxa->scientific_name,"\n";
		print "\thas as descendents :";
		foreach my $new_taxa ( $db->each_Descendent($taxa) ) {
			print "\t",$new_taxa->scientific_name,"\n"
		}
		push @new_taxons, $db->each_Descendent($taxa);
	}
	@taxons = ();
	@taxons = @new_taxons;
} until (scalar @taxons == 0);

print_OUT("Done");


exit;
