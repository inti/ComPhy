#!/usr/bin/perl -w
use strict;
use Bio::TreeIO; 
use Bio::DB::Taxonomy;
use Bio::Phylo::IO 'parse';
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use PDL;
use MLDBM qw(DB_File Storable);
use Fcntl;


# local modules

use PhyloStratiphytUtils;

our (   $help, $man, $tax_folder, $blast_out, $blast_format, $query_taxon,
        $nucl_only, $prot_only, $seq_to_gi, $tax_info, $use_coverage);

GetOptions(
    'help' => \$help,
    'man' => \$man,
    'blast=s@' => \$blast_out,
    'tax_folder=s' => \$tax_folder,
    'blast_format=s' => \$blast_format,
    'nucl_only' => \$nucl_only,
    'prot_only' => \$prot_only,
    'seq_to_gi=s' => \$seq_to_gi,
    'tax_info=s' => \$tax_info,
    'use_coverage' => \$use_coverage,
    'query_taxon=s' => \$query_taxon,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $man);

$seq_to_gi ||= "seq_to_gi";
$tax_info ||= "tax_info";

print_OUT("Parsing taxonomy information");

my $seq_id_files = 'both';
if ($nucl_only) {
    $seq_id_files = "nucl";
    print_OUT("Only working with protein sequences");
} elsif ($prot_only) {
    $seq_id_files = "prot";        
    print_OUT("Only working with nucleic acid sequences");
}
print_OUT("Mapping sequence ids to taxonomy ids");

my (%seq_to_tax_db,%tax_info_db,%tax_tree_db);
if (defined $tax_folder){
    print_OUT("Openning databases");
    print_OUT("   '-> [ $tax_folder/$seq_to_gi.db ]");

    tie %seq_to_tax_db, "DB_File", "$tax_folder/$seq_to_gi.db" or die "Cannot open db file [ $tax_folder/$seq_to_gi.db ]: $!\n";
    
} else {
    print_OUT("Taxonomy information will downloaded and stored on folder [ data ]");
    $tax_folder = "data";
    %seq_to_tax_db  = build_database($seq_to_gi,$tax_folder,$seq_id_files);   
}

print_OUT("   '-> Reading phylogenetic tree and species information");
my $nodesfile = $tax_folder . "nodes.dmp";
my $namefile = $tax_folder . "names.dmp";
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -nodesfile => $nodesfile, -namesfile => $namefile);

#map { 
#    print $_,"\t", $seq_to_tax_db{$_},"\n";
#    my $taxon =   $db->get_taxon(-taxonid => $seq_to_tax_db{$_});
#    print print "id is ", $taxon->id, "\n"; # 9606
#    print "rank is ", $taxon->rank, "\n"; # species
#    print "scientific name is ", $taxon->scientific_name, "\n"; # Homo sapiens
#    print "division is ", $taxon->division, "\n"; # Primates
#    getc; 
#} keys %seq_to_tax_db;


### define some variables to start storing the results

my %S = (); # hash will store to score for each species.


print_OUT("Starting to parse blast output");
$seq_to_tax_db{"322792145"} = 13686;
foreach my $file (@$blast_out){
    #print_OUT("   '-> [ $file ]");
    open (FILE,$file)or die $!;
    my %fields = ();
    while (my $line = <FILE>){
        if ($line =~ m/^#/){
            if ($line =~ m/^# Fields:/){
                # clean up a bit the line
                $line =~ s/^# Fields: //;
                $line =~ s/\s$//;
                $line =~ s/\.//g;
                $line =~ s/\%/percent/g;
                # split the fields names on the comas.
                my @fs = split(/\,\s/,$line);
                # store the field names and positions on a hash
                for (my $i = 0; $i < scalar @fs; $i++) {
                    $fs[$i] =~ s/\s+/_/g;
                    $fields{$fs[$i]} = $i;
                }
            }
            next;
        }
        chomp($line);
        my @data = split(/\t+/,$line);
        # remove the frame of the db hit, we only need to know about the frame of the query seq
        if (exists $fields{'query/sbjct_frames'}){
            $data[ $fields{'query/sbjct_frames'} ] =~ s/\/\w+$//;
        }
#print join "\n", keys %fields,"\n";
#        print $data[ $fields{'subject_id'}],"\n";
        my @subject_id = split(/\|/,$data[ $fields{'subject_id'}]);
        my $coverage = 1;
        if (defined $use_coverage) {
            $coverage = ($data[ $fields{'s_end'}] - $data[ $fields{'s_start'}])/$data[ $fields{'query_length'}]; 
        }
        next if (not defined $seq_to_tax_db{$subject_id[1]});
        $S{ $seq_to_tax_db{$subject_id[1]} } += (-log ($data[ $fields{'evalue'}])) * $coverage;


#        print $subject_id[1],"\t",$seq_to_tax_db{$subject_id[1]},"\n";
#        my $taxon =   $db->get_taxon(-taxonid => $seq_to_tax_db{$subject_id[1]});
#        print "id is ", $taxon->id, "\n"; # 9606
#        print "rank is ", $taxon->rank, "\n"; # species
#        print "scientific name is ", $taxon->scientific_name, "\n"; # Homo sapiens
#        print "division is ", $taxon->division, "\n"; # Primates
        # add hit score to the target specie
    }
}

print Dumper(%S);

print_OUT("Finished processing blast output");

my $main_taxon = $db->get_taxon(-taxonid => $query_taxon);

print_OUT("Starting to calculate PhyloStratus scores");
print_OUT("Identifiying last common ancestors between [ " . $main_taxon->scientific_name . " ] and [ " . scalar (keys %S) . " ] target species");

# get target species and add the query specie
my @species_names = map { $db->get_taxon(-taxonid => $_)->scientific_name;  } keys %S;
# obtain a tree containing only the species of interest
my $tree = $db->get_tree((@species_names,$main_taxon->scientific_name));
# remove redundant nodes, i.e., those with only one ancestor AND one descentdant.
$tree->contract_linear_paths;
#my $out = new Bio::TreeIO(-fh => \*STDOUT, -format => 'newick');
#$out->write_tree($tree),"\n";

my $qry_node = $tree->find_node($main_taxon->id);

my %LCA = (); # this hash stores the the last-common ancestor between the query species and the target specie. This node represent oldest node to which the score of the species needs to be addedd.
my %S_f = ();
foreach my $target_node ($tree->get_nodes){
    next if ($target_node->id eq $qry_node->id);
    next unless (scalar $target_node->each_Descendent() == 0);
    my $lca = $tree->get_lca(($target_node,$qry_node));
    print $qry_node->scientific_name," ",$target_node->scientific_name," ",$lca->scientific_name,"\n";
    $LCA{$target_node->id} = $lca->id;
    do {
        my $anc = $target_node->ancestor();
        $S_f{ $anc->id } += $S{ $target_node->id() };
        print "\t",$anc->scientific_name,"\t", $S_f{ $anc->id },"\n";
       next if ( $lca->id() eq $anc->id  );
    }
}

print_OUT("Done");


exit;

