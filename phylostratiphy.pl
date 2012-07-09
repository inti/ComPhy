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

use constant e_constant => log(10);
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

# define the ids presents of the taxo-to-id mapping files
my $seq_id_files = 'both';
if ($nucl_only) {
    $seq_id_files = "nucl";
    print_OUT("Only working with protein sequences");
} elsif ($prot_only) {
    $seq_id_files = "prot";        
    print_OUT("Only working with nucleic acid sequences");
}
print_OUT("Mapping sequence ids to taxonomy ids");

# load sequence to gene-id (gi) mappings
my (%seq_to_tax_db,%tax_info_db,%tax_tree_db);
if (defined $tax_folder){
    if (not -e "$tax_folder/$seq_to_gi.db") {
        print_OUT("Moving into [ $tax_folder ] to create db");
        chdir($tax_folder);
        unlink("$seq_to_gi.db");
        %seq_to_tax_db  = build_database($seq_to_gi,$tax_folder,$seq_id_files);   
        chdir("../");
    }
    print_OUT("Openning databases");
    print_OUT("   '-> [ $tax_folder/$seq_to_gi.db ]");
    tie %seq_to_tax_db, "DB_File", "$tax_folder/$seq_to_gi.db" or die "Cannot open db file [ $tax_folder/$seq_to_gi.db ]: $!\n";
    
} else {
    print_OUT("Taxonomy information will downloaded and stored on folder [ data ]");
    $tax_folder = "data";
    unless (-d $tax_folder){ mkdir($tax_folder);}
    chdir($tax_folder);
    %seq_to_tax_db  = build_database($seq_to_gi,$tax_folder,$seq_id_files);   
    chdir("../");
}

# load tree of life information, both node' connections and names of nodes.
print_OUT("   '-> Reading phylogenetic tree and species information");
my $nodesfile = $tax_folder . "nodes.dmp";
my $namefile = $tax_folder . "names.dmp";
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -nodesfile => $nodesfile, -namesfile => $namefile);


### define some variables to start storing the results

my %S = (); # hash will store to score for each species.


print_OUT("Starting to parse blast output");
$seq_to_tax_db{"322792145"} = 13686; # this is to be removed 

my %target_taxons = ();
my $seq_counter = 0;

# loop over blast results.
# for each hit we will store the log(p-value) of the blast hit for each taxon of the tartget sequence.
foreach my $file (@$blast_out){
    #print_OUT("   '-> [ $file ]");
    open (FILE,$file)or die $!;
    my %fields = ();
    while (my $line = <FILE>){
        # if line has the header of the table. split the header and keep the columns names.
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
        # extract the indetifiers of the target sequence
        my @subject_id = split(/\|/,$data[ $fields{'subject_id'}]);
        next if (not exists $seq_to_tax_db{$subject_id[1]} ); # to be removed later
        next if (not defined $db->get_taxon(-taxonid => $seq_to_tax_db{$subject_id[1]})); # to avoid crash later when recovering the taxon objects.
        # define coverage as the fraction of the target sequence (subject) covered by the query sequence
        my $coverage = 1;
        if (defined $use_coverage) {
            $coverage = ($data[ $fields{'s_end'}] - $data[ $fields{'s_start'}])/$data[ $fields{'query_length'}]; 
        }
        # get the p-value for the hit from the e-value.
        my $e_value = $data[ $fields{'evalue'}];
        my $p_value = 1 - e_constant**(-$e_value);
        if ($p_value == 0){
            # if p-value is equal 0 the log p-value is NA. then give the equivalane of a p-value of 1x10^-500.
            # this number should be change later
            $S{ $data[ $fields{'query_id'}] }{ $seq_to_tax_db{$subject_id[1]} } += 500;
        } else {
            # otherwise store the log of the p-value
            $S{ $data[ $fields{'query_id'}] }{ $seq_to_tax_db{$subject_id[1]} } += (-log $p_value) * $coverage;
        }
        # add this taxon to the list of taxons targeted by the querys.
        $target_taxons{$seq_to_tax_db{$subject_id[1]}} = "";
    }
}
print_OUT("Finished processing blast output: [ $seq_counter ] sequences of which [ " . scalar (keys %S) . " ] have hits");

# for easy operation store the score of each gene on each specie on a matrix. Later internal nodes of the tree will be added as additional columns, that will make calculation of scores for new columns (internal nodes) faster.


# define some data structures to hold the data.
my $S_g_f = [];
my %S_g_f_gene_idx = ();
my %S_g_f_taxon_idx = ();
my $gene_counter = 0;
my $taxon_counter = 0;

# generate an internal taxon id which we will use for matrix indexes.
foreach my $spc (sort {$a cmp $b} keys %target_taxons){
    $S_g_f_taxon_idx{$spc} = $taxon_counter++;
}

while (my ($gene, $target_species) = each %S){
    # generate a gene counter that will be used for matrix indexes.
    $S_g_f_gene_idx{$gene} = $gene_counter++;
    # loop over the taxons targeted by all queries gene.
    foreach my $spc (sort {$a cmp $b} keys %target_taxons)  {
        if (not exists $target_species->{$spc}) { 
            # if gene had not hit on this taxon then the score is 0
            $S_g_f->[ $S_g_f_gene_idx{$gene}  ] [ $S_g_f_taxon_idx{$spc} ] = 0;
        } else {
            # otherwise the score is the sume of the log(p-values) of the hits on this taxon
            $S_g_f->[ $S_g_f_gene_idx{$gene}  ] [ $S_g_f_taxon_idx{$spc} ] = $target_species->{$spc};
        }
    }
}

# create matrix in PDL format to store the scores of each gene on each taxon.
$S_g_f = mpdl $S_g_f;
$S_g_f /=  $S_g_f->xchg(0,1)->sumover;

# get taxon information for the taxon of the query sequences.
my $main_taxon = $db->get_taxon(-taxonid => $query_taxon);

print_OUT("Starting to calculate PhyloStratum Scores");
print_OUT("Identifiying last common ancestors between [ " . $main_taxon->scientific_name . " ] and [ " . scalar (keys %target_taxons) . " ] target taxons");

# get target species and add the query specie
my @species_names = map { $db->get_taxon(-taxonid => $_)->scientific_name;  } keys %target_taxons;
# obtain a tree containing only the species of interest
my $tree = $db->get_tree((@species_names,$main_taxon->scientific_name));
# remove redundant nodes, i.e., those with only one ancestor AND one descentdant.
$tree->contract_linear_paths;

# get the node for the taxon of interest
my $qry_node = $tree->find_node($main_taxon->id);

my %LCA = (); # this hash stores the the last-common ancestor between the query species and the target specie. This node represent oldest node to which the score of the species needs to be addedd.

# get the root of the tree
my $tree_root = $tree->get_root_node;
# loop over each of the tree leaves, i.e., taxons with blast hits
foreach my $node_id (@{ return_all_Leaf_Descendents($tree_root) }){
    # skip if the leave is the taxon of interest
    next if ($node_id->id == $main_taxon->id);
    # get the last common ancestor of the target and taxon of interest
    my $lca = $tree->get_lca(($node_id,$qry_node));
    $LCA{$lca->id} = $lca;
    # get the path to the root of the target taxon
    my @path = reverse $tree->get_lineage_nodes($node_id);
    # loop from the target taxon to the LCA and add the target taxon score to each ancestor
    foreach my $ancestor ( @path){
        if (exists $S_g_f_taxon_idx{$node_id->ancestor->id}){
            my $matrix_idx = $S_g_f_taxon_idx{$node_id->ancestor->id};
            $S_g_f(,$matrix_idx) += $S_g_f(,$matrix_idx);
        } else {
            # record its position on the matrix
            $S_g_f_taxon_idx{$node_id->ancestor->id} = $taxon_counter++;
            # create new column
            $S_g_f = $S_g_f->glue(1, $S_g_f(,$S_g_f_taxon_idx{$node_id->id}) );
        }
        last if ($ancestor->id == $lca->id);
    }
}



my $last_recorded_lca = '';

my @path_to_taxon = reverse $tree->get_lineage_nodes($main_taxon);
for (my $i = 0; $i < scalar @path_to_taxon; $i++){
    $last_recorded_lca = $i if (exists $S_g_f_taxon_idx{ $path_to_taxon[$i]->id() });
    next if ($last_recorded_lca eq '');
    if (not exists $S_g_f_taxon_idx{ $path_to_taxon[$i]->id() }){
        $S_g_f_taxon_idx{ $path_to_taxon[$i]->id() } = $S_g_f_taxon_idx{ $path_to_taxon[ $last_recorded_lca ]->id };
    }
}

my $I_g_f;
for (my $i = 0; $i < scalar @path_to_taxon; $i++){
    next if (not exists $S_g_f_taxon_idx{ $path_to_taxon[$i]->id() });
    last if (not defined $path_to_taxon[$i]->ancestor);
    my $idx_current = $S_g_f_taxon_idx{ $path_to_taxon[$i]->id() };
    my $idx_previous = $S_g_f_taxon_idx{ $path_to_taxon[$i]->ancestor->id() };
    my $tmp_col = $S_g_f(,$idx_current) - $S_g_f(,$idx_previous);
    push @{ $I_g_f }, [$tmp_col->list];
}
$I_g_f = pdl $I_g_f;
my $S_f = $I_g_f->sumover; # here store the cummulative score for each stratum over all genes.

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

sub nr_array {
    my %tmp = ();
    map { $tmp{$_} = ''; } @_;
    return(keys %tmp); 
}
