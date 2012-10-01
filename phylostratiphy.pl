#!/usr/bin/perl -w
use strict;
use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;
use Bio::DB::EUtilities;
use Getopt::Long;
use Pod::Usage;
use PDL;
use PDL::Matrix;
use PDL::NiceSlice;
use Data::Dumper;

use constant E_CONSTANT => log(10);
my $EMAIL = 'intipedroso@gmail.com';
# local modules

use PhyloStratiphytUtils;

our (   $help, $man, $tax_folder, $blast_out, $blast_format, $user_provided_query_taxon_id, $out,
        $use_coverage, $hard_threshold, $soft_threshold, $gi_tax_id_info,$blastdbcmd,
        $seq_db, $not_use_ncbi_entrez );

GetOptions(
    'help' => \$help,
    'man' => \$man,
    'blast=s@' => \$blast_out,
    'tax_folder=s' => \$tax_folder,
    'blast_format=s' => \$blast_format,
    'use_coverage' => \$use_coverage,
    'query_taxon=s' => \$user_provided_query_taxon_id,
    'out|o=s' => \$out,
    'hard_threshold|hard=f' => \$hard_threshold,
    'soft_threshold|soft' => \$soft_threshold,
    'gi_tax_id=s' => \$gi_tax_id_info,  # dysbindin.tax_info.csv
    'blastdbcmd=s' => \$blastdbcmd,
    'seq_db|db=s' => \$seq_db,
    'no_ncbi_entrez' => \$not_use_ncbi_entrez,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $man);


defined $blastdbcmd or $blastdbcmd = `which blastdbcmd`;
chomp($blastdbcmd);
defined $seq_db or $seq_db = "nr"; # assuming proteins and that path to dbs is on a enviromental variable
defined $hard_threshold or $hard_threshold = 1e-3;
defined $soft_threshold and $hard_threshold = undef;
defined $blast_format or $blast_format = 'table';
print_OUT("Reading taxonomy information");

# load tree of life information, both node' connections and names of nodes.
print_OUT("   '-> Reading phylogenetic tree and species information");
my $nodesfile = $tax_folder . "nodes.dmp";
my $namefile = $tax_folder . "names.dmp";
my $taxNCBI = Bio::LITE::Taxonomy::NCBI->new( db=>"NCBI", names=> $namefile, nodes=>$nodesfile, dict=>"$tax_folder/gi_taxid_prot.bin");
my @ql = $taxNCBI->get_taxonomy( $user_provided_query_taxon_id);

### define some variables to start storing the results

print_OUT("Starting to parse blast output");

my %S = (); # hash will store to score for each species.
# loop over blast results.
# for each hit we will store the log(p-value) of the blast hit for each taxon of the tartget sequence.
foreach my $file (@$blast_out){
    #print_OUT("   '-> [ $file ]");
    if ($blast_format eq 'table'){
        my $parsed_blast_out = parse_blast_table($file);
        @S{keys %{$parsed_blast_out}} = values %{$parsed_blast_out};
    }
#    open (FILE,$file)or die $!;
#    my %fields = ();
#    while (my $line = <FILE>){
#        # if line has the header of the table. split the header and keep the columns names.
#        if ($line =~ m/^#/){
#            if ($line =~ m/^# Fields:/){
#                # clean up a bit the line
#                $line =~ s/^# Fields: //;
#                $line =~ s/\s$//;
#        $line =~ s/\.//g;
#        $line =~ s/\%/percent/g;
#        # split the fields names on the comas.
#        my @fs = split(/\,\s/,$line);
#        # store the field names and positions on a hash
#        for (my $i = 0; $i < scalar @fs; $i++) {
#            $fs[$i] =~ s/\s+/_/g;
#            $fields{$fs[$i]} = $i;
#        }
#        $seq_counter++;
#    }
#    next;
#}
#chomp($line);
#my @data = split(/\t+/,$line);
## remove the frame of the db hit, we only need to know about the frame of the query seq
#if (exists $fields{'query/sbjct_frames'}){
#    $data[ $fields{'query/sbjct_frames'} ] =~ s/\/\w+$//;
#}
## extract the indetifiers of the target sequence
#my @subject_id = split(/\|/,$data[ $fields{'subject_id'}]);
## define coverage as the fraction of the target sequence (subject) covered by the query sequence
#my $coverage = 1;
#if (defined $use_coverage) {
#    $coverage = ($data[ $fields{'s_end'}] - $data[ $fields{'s_start'}])/$data[ $fields{'query_length'}];
#}
## if using hard threshold then remove hits by e-value
#if (defined $hard_threshold){
#    next if ($data[ $fields{'evalue'}] > $hard_threshold);
#}
## get the taxid and lineage for the sequence.
#if (not defined $gi_taxData{ $subject_id[1] }){
#    my $sbjct_taxid = $taxNCBI->get_taxid( $subject_id[1] );
#    if ($sbjct_taxid == 0){
#        push @{ $ids_not_found{ $data[ $fields{'subject_id'}] }}, $data[ $fields{'query_id'}];
#    } else {
#        $gi_taxData{ $subject_id[1] } = { 'lineage' => [], 'tax_id' => $sbjct_taxid, 'lca_with_qry' => ""};
#    }
#}
#if (exists $gi_taxData{ $subject_id[1] }){
#    if (not exists $lineages{ $gi_taxData{ $subject_id[1] }->{'tax_id'} }){
#        my @sbjct_lineage = $taxNCBI->get_taxonomy( $gi_taxData{ $subject_id[1] }->{'tax_id'} );
#        if (scalar @sbjct_lineage == 0){
#            push @{ $ids_not_found{ $data[ $fields{'subject_id'}] }}, $data[ $fields{'query_id'}];
#        } else {
#            my $lca = get_lca_from_lineages(\@sbjct_lineage,\@ql); # need double checking on the ones that do not give match
#            next if ($lca eq "diff_root"); # exclude those that have a different root to cell organisms.
#            $gi_taxData{ $subject_id[1] }->{'lineage'} = \@sbjct_lineage;
#            $gi_taxData{ $subject_id[1] }->{'lca_with_qry'} = $lca;
#            $lineages{ $gi_taxData{ $subject_id[1] }->{'tax_id'} }->{'lca_with_qry'} = $lca;
#            $lineages{ $gi_taxData{ $subject_id[1] }->{'tax_id'} }->{'lineage'} = \@sbjct_lineage;
#        }
#    } else {
#        $gi_taxData{ $subject_id[1] }->{'lineage'} = $lineages{ $gi_taxData{ $subject_id[1] }->{'tax_id'} }->{'lineage'} ;
#        $gi_taxData{ $subject_id[1] }->{'lca_with_qry'} = $lineages{ $gi_taxData{ $subject_id[1] }->{'tax_id'} }->{'lca_with_qry'} ;
#    }
#}
## get the p-value for the hit from the e-value.
#my $e_value = pdl $data[ $fields{'evalue'}];
#my $p_value = pdl E_CONSTANT**(-$e_value);
#$p_value = pdl 1 - $p_value;
#$p_value = pdl $e_value if ($p_value == 0);
#my $score = -1*($p_value->log) * $coverage;
#push @{ $S{ $data[ $fields{'query_id'}] }  }, { 'subject_id' => $subject_id[1],
#    'score' => $score,
#    'p_value' => $p_value,
#    'e_value' => $e_value };
#}
}

print_OUT("Finished processing blast output with results for [ " . scalar (keys %S) . " ] sequences.");

######## FILTER BLAST HITS AND GET TAXONMY INFORMATION FOR TARGET SEQUENCES AND SPECIES.
my %gi_taxData = (); # store taxonomy information for each target gi
my %lineages = (); # each entry has a 'lca_with_qry' and a 'lineage'
my $seq_counter = 0;
my %ids_not_found = (); # store ids of target sequences without taxonomy information on local DBs

while (my ($qry_seq,$blast_subjects) = each %S){
    # loop over blast results for this query sequence
    my $target_counter = -1; # set to -1 so that first item sets it to 0.
    foreach my $target_seqs (@{$blast_subjects}){
        $target_counter++;
        ## Filter blast results.
        # if using hard threshold then remove hits by e-value
        if (defined $hard_threshold){
            if ($target_seqs->{'evalue'} > $hard_threshold){
                # flag this target sequence as not passing the filter
                $S{ $qry_seq }->[ $target_counter ]->{'pass_hard_thresold'} = 0;
                # if skip the taxonomy of this seq unless we are interested in doing soft-threshold method.
                unless (defined $soft_threshold ){
                    next;
                }
            } else {
                $S{ $qry_seq }->[ $target_counter ]->{'pass_hard_thresold'} = 1;
            }
        }
        if (not defined $gi_taxData{ $target_seqs->{'subject_id'} }){
            my $sbjct_taxid = $taxNCBI->get_taxid( $target_seqs->{'subject_id'} );
            if ($sbjct_taxid == 0){
                push @{ $ids_not_found{ $target_seqs->{'subject'} }}, $qry_seq;
            } else {
                $gi_taxData{ $target_seqs->{'subject_id'} } = { 'lineage' => [],
                                                                'tax_id' => $sbjct_taxid,
                                                                'lca_with_qry' => ''};
                if (not exists $lineages{ $sbjct_taxid }){
                    my @sbjct_lineage = $taxNCBI->get_taxonomy( $sbjct_taxid );
                    if (scalar @sbjct_lineage == 0){
                        push @{ $ids_not_found{ $target_seqs->{'subject'} }}, $qry_seq;
                    } else {
                        my $lca = get_lca_from_lineages(\@sbjct_lineage,\@ql);
                        if ($lca eq "diff_root") {# exclude those that have a different root to cell organisms.
                            $gi_taxData{ $target_seqs->{'subject_id'} }->{'lineage'} = \@sbjct_lineage;
                            $gi_taxData{ $target_seqs->{'subject_id'} }->{'lca_with_qry'} = $lca;
                            next;
                        }
                        $gi_taxData{ $target_seqs->{'subject_id'} }->{'lineage'} = \@sbjct_lineage;
                        $gi_taxData{ $target_seqs->{'subject_id'} }->{'lca_with_qry'} = $lca;
                        $lineages{ $sbjct_taxid }->{'lca_with_qry'} = $lca;
                        $lineages{ $sbjct_taxid }->{'lineage'} = \@sbjct_lineage;
                    }
                } else {
                    $gi_taxData{ $target_seqs->{'subject_id'} }->{'lineage'} = $lineages{ $sbjct_taxid }->{'lineage'} ;
                    $gi_taxData{ $target_seqs->{'subject_id'} }->{'lca_with_qry'} = $lineages{ $sbjct_taxid }->{'lca_with_qry'} ;
                }
            }
        }
    }

}


######## FIND TAXONOMY INFORMATION FOR IDS WITH INCOSISTENT INFORMATION ON LOCAL FILES ################
# compile all ids for which taxonomy information was not found with local files.
my @accs = ();
my %accs_to_gi = ();
foreach my $hidden (keys %ids_not_found){
    my @subject_id = split( /\|/,$hidden );
    push @accs, ($subject_id[3] =~ m/(.*)\.\d+$/);
    $accs_to_gi{ $accs[-1] } = $subject_id[1];
}
print_OUT("There were [ " . scalar @accs . " ] subject ids from blast search unmatched with local taxonomy DBs");
if (defined $not_use_ncbi_entrez){
    print_OUT("   '-> Printing this ids to [ $out.ids_without_taxonomy_information.txt ].");
    open(OUT,">$out.ids_without_taxonomy_information.txt") or die $!;
    print OUT join "\n", @accs;
    close(OUT);
} else {
    print_OUT("   '-> Using NCBI webservices to get taxonomy information on them.");
    if (scalar @accs > 0){
        my $factory = Bio::DB::EUtilities->new( -eutil => 'esearch',
                                                -email => $EMAIL,
                                                -db    => 'protein',
                                                -retmax     => 10*(scalar @accs),
                                                -term  => join(',',@accs),
                                                -usehistory => 'y');
        
        my $hist = $factory->next_History || die print_OUT("No history data returned from esearch");
        my @uids = $factory->get_ids;
        $factory->set_parameters(   -eutil => 'esummary',
                                    -db => 'protein',
                                    -id => \@uids,
                                    -email => $EMAIL,
                                    -history => $hist);
    #    $factory->reset_parameters(-eutil => 'esummary', -db => 'protein', -id => \@uids,-email => $EMAIL);
        while (my $ds = $factory->next_DocSum) {
            my ($taxid) = $ds->get_contents_by_name("TaxId");
            my ($caption) = $ds->get_contents_by_name("Caption");
            if (not defined $lineages{ $taxid }){
                my @sbjct_lineage = $taxNCBI->get_taxonomy( $taxid );
                my $lca = get_lca_from_lineages(\@sbjct_lineage,\@ql); # need double checking on the ones that do not give match
                next if ($lca eq "diff_root"); # exclude those that have a different root to cell organisms.
                $gi_taxData{ $accs_to_gi{ $caption } }->{'lineage'} = \@sbjct_lineage;
                $gi_taxData{ $accs_to_gi{ $caption } }->{'lca_with_qry'} = $lca;
                $lineages{ $taxid }->{'lineage'} =  \@sbjct_lineage;
                $lineages{ $taxid }->{'lca_with_qry'} = $lca;
            } else {
                $gi_taxData{ $accs_to_gi{ $caption } }->{'lineage'} = $lineages{ $taxid }->{'lineage'};
                $gi_taxData{ $accs_to_gi{ $caption } }->{'lca_with_qry'} = $lineages{ $taxid }->{'lca_with_qry'};
            }
        }
    }
    print_OUT("   '-> done with NCBI webservices.")
}

############## PROCEED TO FINISH CALCULATIONS ##################
# Get the oldest stratum for every gene
my %qry_ancestors_ranks = ();
my $c = 0;
map { $qry_ancestors_ranks{$_} = $c++;   } @ql;
my %PhyloStratum_scores = ();
foreach my $qry_seq (keys %S) {
    my $oldest_stratum = scalar @ql;
    foreach my $target_seq ( @{ $S{ $qry_seq } } ){
        #print Dumper($target_seq);
        #        getc;
        next if ($target_seq->{'pass_hard_thresold'} == 0);
        if (not exists $gi_taxData{ $target_seq->{'subject_id'}}){
            print_OUT("Something went wrong. I cannot find tax infor for this sequence after having done all the work. [ $target_seq->{'subject'} ]");
            next;
        }
        my $lca = $gi_taxData{ $target_seq->{'subject_id'} }->{'lca_with_qry'};
        next if ($lca eq "diff_root");
        my $target_seq_stratum = $qry_ancestors_ranks{ $lca };
        if ($target_seq_stratum < $oldest_stratum) {
            $oldest_stratum = $target_seq_stratum ;
        }
    }
    my @phyloScores = list zeroes scalar @ql;
    $phyloScores[ $oldest_stratum ] = 1;
    $PhyloStratum_scores{$qry_seq} = \@phyloScores;
}
print_OUT("Finished calculating hard coded scores");
print_OUT("Printing summary scores to [ $out.qry_node_phylostratumscores.txt ].");

# Print out the PhyloStratumScores as hardcoded
open(OUT,">$out.qry_node_phylostratumscores.txt") or die $!;
print OUT join "\t", @ql;
print OUT "\n";
print OUT join "\t", list sumover mpdl values %PhyloStratum_scores;
print OUT "\n";
close(OUT);

print_OUT("Printing gene phylostratum scores to [ $out.txt ].");

# Make table with scores for each gene
my $phyloScoresTable = join "\t", ("ID",@ql);
$phyloScoresTable .= "\n";

while (my ($qry_id,$qry_scores) = each %PhyloStratum_scores) {
    $phyloScoresTable .= join "\t", ($qry_id,@$qry_scores);
    $phyloScoresTable .= "\n";
}
open(OUT,">$out.txt") or die $!;
print OUT $phyloScoresTable;
close(OUT);

print_OUT("Done");

exit;

# TODO:
# 3- implement soft-coded
# 4- guess specie


__END__

=head1 NAME
 
 Perl implementation of the PhyloStratigraphy comparative genomics methods.
 
=head1 DESCRIPTION

B<This program> will perform a PhyloStratigraphy analysis. It provides a implementation and extensions of the original methodology of Domazet-Loso et. al. (2003).
 Please see README for information on Perl module requirments and some additional information.
 
 Comments, suggestions or complains you be addressed to IntiPedroso@gmail.com
 
 References:
  Domazet-Loso, T., and Tautz, D. (2003). An evolutionary analysis of orphan genes in Drosophila. Genome Res. 13, 2213-2219.

=head1 SYNOPSIS

    ./phylostratiphy.pl [options]

    General options
    -h, --help		print help message
    -m, --man		print complete documentation
 
    Input Files
    -blast         Output from blast (Required)
    -tax_folder    Folder with taxonomy information (Required)
    -blast_format  Format of blast output (Not yet functional)
    -query_taxon   Taxnomy identifier of query specie. (Required)
    -seq_db        Path or name to sequence db in blast format. (Optional)
    -blastdbcmd    Path to blastdbcmd executable (Optional)
    -gi_tax_id     Tabe with taxonomy information for the target sequences of the blast output (Optional)
 
    Output Files
    -o, --o        Name of output files (Required)

    Analysis options
    -hard, --hard_threshold    Switches to original methodology of Domazet-Loso et. al. (2003)
    -soft, --soft_threshold    Uses a softhreshols strategy that aims to account for coverage and uncertaintiy of blast results.
    -use_coverage  Means that scores are weightes by coverage of sequence alignment. Not compatible with -hard.
    -no_ncbi_entrez     Do not use NCBI Entrez API to get information on sequence ids without data on local DBs.

 
 
=head1 OPTIONS

=over 8

=item B<-help>

Print help message

=item B<-man>

Print complete documentation

=item B<-blast>
 Output from blast
 
=item B<-tax_folder>
 
 Folder with taxonomy information
 
=item B<-blast_format>  
 
 Format of blast output

=item B<-virus_list>    
 
 List with ids of viral taxons
 
=item B<-query_taxon>   
 
 Taxnomy identifier of query specie. User must use the id provided by the NCBI Taxonomy database, e.g, human = 9606 and Acromyrmex echinatior = 103372.
 
=item B<-seq_db> 
 
 Path or name to sequence db in blast format. If nor provided it will be set to 'nr' and it will be assumed that you have set the  ~/.ncbirc as indicated on the BLAST documentation with the path to the folder where the sequence databased are stored, i.e., the $BLASTDB enviromental variable has been set.

=item B<-blastdbcmd>

 Path to blastdbcmd executable. If not provided it will be assumed that it can found with the command `which blastdbcmd`

=item B<-gi_tax_id>
 
 Tabe with taxonomy information for the target sequences of the blast output. This table has the format obtained by extrating the taxonomy information the sequence db with the following command (as it is done internally by the software)
 >blastdbcmd -outfmt "%a,%g,%T,%L,%S" -entry_batch file_with_sequence_ids.txt -db seq_db -out output_file.csv
 where seq_db is whatever the -seq_db as been set to and output_file.csv is the file that you need to provide to the -gi_tax_id option. If this option is not provided this command will be run internally and you will get the file as "out".gi_tax_id.csv, where "out" is the whatever you provide to the option -o, --out

=item B<-o, --out>
 
 Name of output files

=item B<-hard, --hard_threshold>
 
 Switches to original methodology of Domazet-Loso et. al. (2003)

=item B<-soft, --soft_threshold>
 
 Uses a softhreshols strategy that aims to account for coverage and uncertaintiy of blast results.

=item B<-use_coverage>
 
 Means that scores are weightes by coverage of sequence alignment. Not compatible with -hard.
 
=item B<-no_ncbi_entrez>
 
 Do not use NCBI Entrez API to get information on sequence ids without data on local DBs.

=back
 
=head2 DESCRIPTION

TODO

=head3 Examples

=over 8

=item B<1. Basic Analysis>

To run the analyses with a small example do

>my_perl phylostratiphy.pl -blast dysbindin.blast_out.txt -tax_folder data/ -query_taxon 9606 -virus_list tmp.virus.txt -out test_phylostratiphy

The example consistes of human sequences and so -query_taxon 9606 corresponds to the human tax id.
MORE TO COME
 
 
=back



=cut

