#!/usr/bin/perl -w
use strict;
use Bio::TreeIO; 
use Bio::DB::Taxonomy;
use Bio::Phylo::IO 'parse';
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use PDL;
use PDL::Matrix;
use PDL::NiceSlice;
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

my %target_taxons = ();
my $seq_counter = 0;
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
            $seq_counter++;
            next;
        }
        chomp($line);
        my @data = split(/\t+/,$line);
        # remove the frame of the db hit, we only need to know about the frame of the query seq
        if (exists $fields{'query/sbjct_frames'}){
            $data[ $fields{'query/sbjct_frames'} ] =~ s/\/\w+$//;
        }
        my @subject_id = split(/\|/,$data[ $fields{'subject_id'}]);
        my $coverage = 1;
        if (defined $use_coverage) {
            $coverage = ($data[ $fields{'s_end'}] - $data[ $fields{'s_start'}])/$data[ $fields{'query_length'}]; 
        }
        if ($data[ $fields{'evalue'}] == 0){
            $S{ $data[ $fields{'query_id'}] }{ $seq_to_tax_db{$subject_id[1]} } += 500;
        } else {
            $S{ $data[ $fields{'query_id'}] }{ $seq_to_tax_db{$subject_id[1]} } += (-log ($data[ $fields{'evalue'}])) * $coverage;
        }
        $target_taxons{$seq_to_tax_db{$subject_id[1]}} = "";
    }
}
print_OUT("Finished processing blast output: [ $seq_counter ] sequences of which [ " . scalar (keys %S) . " ] have hits");

# for easy operation store the score of each gene on each specie on a matrix. Later internal nodes of the tree will be added as additional columns, that will make calculation of scores for new columns (internal nodes) faster.
# print score matrix;

my $S_g_f = [];
my %S_g_f_gene_idx = ();
my %S_g_f_taxon_idx = ();
my $gene_counter = 0;
my $taxon_counter = 0;

foreach my $spc (sort {$a cmp $b} keys %target_taxons){
    $S_g_f_taxon_idx{$spc} = $taxon_counter++;
}
my $max_leaf_idx = $taxon_counter;

while (my ($gene, $target_species) = each %S){
    $S_g_f_gene_idx{$gene} = $gene_counter++;
    foreach my $spc (sort {$a cmp $b} keys %target_taxons)  {
        if (not exists $target_species->{$spc}) { 
            $S_g_f->[ $S_g_f_gene_idx{$gene}  ] [ $S_g_f_taxon_idx{$spc} ] = 0;
        } else {
            $S_g_f->[ $S_g_f_gene_idx{$gene}  ] [ $S_g_f_taxon_idx{$spc} ] = $target_species->{$spc};
        }
    }
}

# create matrix in PDL format
$S_g_f = mpdl $S_g_f;
#$S_g_f /=  $S_g_f->xchg(0,1)->sumover;
print $S_g_f;

# get taxon information for the query taxon.
my $main_taxon = $db->get_taxon(-taxonid => $query_taxon);


print_OUT("Starting to calculate PhyloStratus scores");
print_OUT("Identifiying last common ancestors between [ " . $main_taxon->scientific_name . " ] and [ " . scalar (keys %target_taxons) . " ] target species");

# get target species and add the query specie
my @species_names = map { $db->get_taxon(-taxonid => $_)->scientific_name;  } keys %target_taxons;
# obtain a tree containing only the species of interest
my $tree = $db->get_tree((@species_names,$main_taxon->scientific_name));
# remove redundant nodes, i.e., those with only one ancestor AND one descentdant.
$tree->contract_linear_paths;
#my $out = new Bio::TreeIO(-fh => \*STDOUT, -format => 'newick');
#print_OUT("this is the contracted tree of interest");
#$out->write_tree($tree),"\n";

# get the node for the taxon of interest
my $qry_node = $tree->find_node($main_taxon->id);

my %LCA = (); # this hash stores the the last-common ancestor between the query species and the target specie. This node represent oldest node to which the score of the species needs to be addedd.
my %S_f = (); # here store the cummulative score for each stratum over all genes.
my %internal_descedents = ();
foreach my $target_node ($tree->get_nodes){
    # skip if it is the query specie
    next if ($target_node->id eq $qry_node->id); 
    # skip if the node does not have descendents, i.e., it the leaves. 
    next if ($target_node->is_Leaf() == 1);
    # get the last common ancestor of the two, the query and the target specie.
    # store the name of the LCA
    my $lca = $tree->get_lca(($target_node,$qry_node));
    $LCA{$target_node->id} = $lca;
    $internal_descedents{ $lca->id } = return_all_Leaf_Descendents($lca);
}
# using the matrix indexes of the leaf descendents for each internal node calculate the score for each internal node
my $I_g_f = mzeroes $gene_counter, scalar $tree->get_nodes - 1;
print $I_g_f;
foreach my $node_id ( sort { $tree->find_node($a)->height() <=> $tree->find_node($b)->height() }   keys %internal_descedents){
    my $tree_node = $tree->find_node($node_id);
    my $desc = $internal_descedents{$node_id};
    my $mat_idx = [];
    foreach my $taxon (@$desc){
        next if ($taxon->id == $main_taxon->id);
        push @{$mat_idx}, $S_g_f_taxon_idx{$taxon->id};
    }
    $S_g_f_taxon_idx{$node_id} = $taxon_counter++;
    $mat_idx = pdl $mat_idx;
    $S_g_f = $S_g_f->glue(1, $S_g_f(,$mat_idx)->xchg(0,1)->sumover );
    my $node_descendents_idx = [];
    foreach my $taxon ($tree_node->each_Descendent()){
        next if ($taxon->id == $main_taxon->id);
        push @{$node_descendents_idx},  $S_g_f_taxon_idx{$taxon->id};
    }
    $node_descendents_idx = pdl $node_descendents_idx;
    print $taxon_counter," ",scalar $tree_node->each_Descendent()," ",$node_descendents_idx,"\n";
#    print $S_g_f(,$S_g_f_taxon_idx{$node_id});
#    print $S_g_f( , $node_descendents_idx );
    getc;
    $I_g_f(,$taxon_counter) = $S_g_f(,$S_g_f_taxon_idx{$node_id}) - $S_g_f( , $node_descendents_idx );
#    print $I_g_f;
    getc;
}

print "I_g_f\n";
print $I_g_f;
print "S_g_f\n";
print $S_g_f,"\n";


print_OUT("Done");


exit;

sub return_all_Leaf_Descendents {
    my $taxon = shift;
    my @back = ();
    foreach my $d ($taxon->get_all_Descendents()){
        push @back, $d if ($d->is_Leaf == 1);
    }
    return(\@back);
}
