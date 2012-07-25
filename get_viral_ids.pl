#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Bio::DB::Taxonomy;
# local modules
use PhyloStratiphytUtils;

our (   $help, $man, $tax_folder, $out, $tax_id);

GetOptions(
    'help' => \$help,
    'man' => \$man,
    'tax_folder=s' => \$tax_folder,
    'out|o=s' => \$out,
    'tax_id=i' => \$tax_id,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $man);

defined $tax_id or $tax_id = 10239;

print_OUT("Reading phylogenetic tree and species information");
my $nodesfile = $tax_folder . "nodes.dmp";
my $namefile = $tax_folder . "names.dmp";
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -nodesfile => $nodesfile, -namesfile => $namefile);


my $node = $db->get_Taxonomy_Node(-taxonid => $tax_id);

print_OUT("Fetching extant taxa of tax id [ " . $node->id . " ] which is a [ " . $node->scientific_name . " ] of rank [ " . $node->rank . " ] ");
# to only get children that are of a particular rank in the taxonomy test if their rank is 'species' for example

my @extant_children = grep { $_->rank eq 'species' } $db->get_all_Descendents($node);
print_OUT("Printing out");
open (OUT,">$out") or die $!;
for my $child ( @extant_children ) {
    print OUT $child->ncbi_taxid,",",$child->id, ",", $child->rank, ",",$child->scientific_name, "\n";
}
close(OUT);
print_OUT("Done");
exit;
