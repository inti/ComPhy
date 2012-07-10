package PhyloStratiphytUtils;
use strict;
use warnings;
use Carp;
use Exporter qw (import);

eval { 
    use DB_File;
    use Fcntl;
    use Data::Dumper;
};
if ($@) { 
	print "Some libraries does not seem to be in you system. quitting\n";
	exit(1);
}

our (@EXPORT, @EXPORT_OK, %EXPORT_TAGS);

@EXPORT = qw( parse_gi_taxid_files print_OUT build_database progress_bar read_taxonomy_files);				# symbols to export by default
@EXPORT_OK = qw( parse_gi_taxid_files print_OUT build_database progress_bar read_taxonomy_files);			# symbols to export on request

sub read_taxonomy_files {
    my $tax_folder = shift;
    my (%tax_info_db, %tax_tree_db);
    
    print_OUT("   '-> Working on information linked to species id");
    # first work on the species names.
    print_OUT("   '-> [ names.dmp ]");
    open (DMP,"$tax_folder/names.dmp") or die $!;
    my $n_lines = get_file_number_of_lines("$tax_folder/names.dmp");
    my $counter = 0;
    while (my $line = <DMP>) {
        print scalar localtime," \t", progress_bar(++$counter,$n_lines);
        chomp($line);
        $line =~ s/\t+//g; 
        $line =~ s/|$//; 
        my ($tax_id, $name_txt,$unique_name,$name_class) = split(/\|/,$line);
        if (not defined $tax_info_db{$tax_id}) {
            $tax_info_db{$tax_id} = {   'synonym' => [], 
                'name_txt' => "",
                'rank' => "",
                'embl_code' => "",
                'division_id' => "",
                'comments' => ""
            };
            
        }
        if ( $name_class ne 'scientific name') { 
            @{ $tax_info_db{$tax_id}->{'synonym'} } = $name_txt; 
            next;
        }
        $tax_info_db{$tax_id}->{'name_txt'} = $name_txt;
        if ($unique_name eq ''){ 
            $tax_info_db{$tax_id}->{$name_class} = $name_txt;
        } else {
            $tax_info_db{$tax_id}->{$name_class} = $unique_name;
            push @{ $tax_info_db{$tax_id}->{'synonym'} }, $name_txt;
        }
    }
    print "\n"; # new line for to finish the progress bar.
    
    
    #####
    # now work on the taxonomy itself.
    print_OUT("   '-> Working on taxnomy tree");
    # create database for the taxonomy tree
    open (DMP,"$tax_folder/nodes.dmp") or die $!;
    print_OUT("   '-> [ nodes.dmp ]");
    $n_lines = get_file_number_of_lines("$tax_folder/nodes.dmp");
    $counter = 0;
    while (my $line = <DMP>) {
        print scalar localtime," \t", progress_bar(++$counter,$n_lines);
        chomp($line);
        $line =~ s/\t+//g; 
        $line =~ s/|$//; 
        my ($tax_id,                        #-- node id in GenBank taxonomy database
        $parent_tax_id,                  #-- parent node id in GenBank taxonomy database
        $rank,                           #-- rank of this node (superkingdom, kingdom, ...) 
        $embl_code,                      #-- locus-name prefix; not unique
        $division_id,                    #-- see division.dmp file
        $inherited_div_flag,             #(1 or 0)		-- 1 if node inherits division from parent
        $genetic_code_id,				#-- see gencode.dmp file
        $inherited_GC_flag,              #(1 or 0)		-- 1 if node inherits genetic code from parent
        $mitochondrial_genetic_code_id,	#-- see gencode.dmp file
        $inherited_MGC_flag,             #(1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
        $GenBank_hidden_flag,            #(1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
        $hidden_subtree_root_flag,       #(1 or 0)       -- 1 if this subtree has no sequence data yet
        $comments                       #-- free-text comments and citations 
        )= split(/\|/,$line);
        
        # store child parent relationship
        $tax_tree_db{$tax_id} = { 'parent' => $parent_tax_id };
        # store adition information on the tree.
        $tax_info_db{$tax_id}->{'rank'} = $rank;
        $tax_info_db{$tax_id}->{'embl_code'} = $embl_code;
        $tax_info_db{$tax_id}->{'division_id'} = $division_id;
        $tax_info_db{$tax_id}->{'comments'} = $comments;
    }
    print "\n"; # finish progress bar.
    return(\%tax_info_db,\%tax_tree_db);
}

sub build_database {
    print_OUT("Building database with taxonomy information");
    my $seq_to_tax_db_name = shift;
    my $tax_folder = shift;
    my $seq_id_files_to_parse = shift;

    $seq_id_files_to_parse ||= "both";
    print_OUT("   '-> Working on sequence id to species id mapping");
    # CREATE HAShes to store data
    # create the db on disk.
    print_OUT("Creating DBs on disk");
    
    my %seq_to_tax_db;
    # remove previous db file.
    unlink ("$seq_to_tax_db_name.db") if (-e "$seq_to_tax_db_name.db"); 
    my $SEQ = tie %seq_to_tax_db, "DB_File", "$seq_to_tax_db_name.db", O_RDWR|O_CREAT, 0640 or die "Cannot open file to create db [ $seq_to_tax_db_name ]: $!\n";
    
    # set the sequences id to taxonomy id database.
    my $seq_id_files = [];
    if ($seq_id_files_to_parse eq 'both') {
        $seq_id_files = ["gi_taxid_nucl.dmp","gi_taxid_nucl.dmp"];
    } elsif ($seq_id_files_to_parse eq 'nucl') {
        $seq_id_files = ["gi_taxid_nucl.dmp"];        
    } elsif ($seq_id_files_to_parse eq 'prot') {
        $seq_id_files = ["gi_taxid_prot.dmp"];        
    }
    # read file and store the mapping between ids.
    foreach my $f (@$seq_id_files){
        print_OUT("   '-> [ $f ]");
        open (DMP,$f) or die $!;
        my $n_lines = get_file_number_of_lines($f);
        my $counter = 0;
        while (my $line = <DMP>) {
            #print scalar localtime," \t", progress_bar($.,$n_lines);
            chomp($line);
            my($gi,$tax_id) = split(/[\s+\t+]/,$line);
            $seq_to_tax_db{$gi} = $tax_id;         
        }
        close(DMP);
        print "\n";
    }
    
    return(%seq_to_tax_db);
}

# wget-style. routine by tachyon
# at http://tachyon.perlmonk.org/
sub progress_bar {
    my ( $got, $total, $width, $char ) = @_;
    $width ||= 25; $char ||= '=';
    my $num_width = length $total;
    sprintf "   |%-${width}s| Done with [ %${num_width}s ] of [ %s (%.2f%%) ]\r", 
    $char x (($width-1)*$got/$total). '>', 
    $got, $total, 100*$got/+$total;
}


sub get_file_number_of_lines {
    my $name = shift;
    my $wc_output = `wc -l $name`;
    my @n = split(/\s+/,$wc_output);
    return($n[1]);
}


sub print_OUT {
	my $string = shift;
	my @file_handles = @_; 	
	print scalar localtime(), "\t$string\n";
	unless (scalar @file_handles == 0){
		foreach my $fh (@file_handles){
			print $fh scalar localtime(), "\t$string\n";
		}
	}
}

sub fetch_tax_ids_from_blastdb {
    my $seq_ids = shift;
    print_OUT("Getting taxons ids for hit sequences");
    print_OUT("   '-> Printing ids to temporary file");
    my $tmp_file = "tmp.seq_tax_ids.$$";
    open(IDS,">$tmp_file.txt") or die $!;
    foreach my $id (@$seq_ids){ print IDS $id,"\n"; }
    close(IDS);
    print_OUT("   '-> Running blastdbcmd to get sequences information");
    #`blastdbcmd -outfmt "%a,%g,%T,%L,%S" -entry_batch $tmp_file.txt -db nr -out $tmp_file.csv`;
    #open (TAX_IDS,"$tmp_file.csv") or die $!;
    print_OUT("   '-> Parsing sequence information");
    open (TAX_IDS,"gi_to_tax_id_from_blast.csv") or die $!;
    my %back_gi_to_taxinfo = ();
    my %target_taxons = ();
    my $taxon_counter = 0;
    my %seen_taxon = ();
    while (my $line = <TAX_IDS>){
        chomp($line);
        my @data = split(/,/,$line);
        $back_gi_to_taxinfo{$data[1]} = {'accession' => $data[0], 'gi' => $data[1],'taxid' => $data[2], 'common_tax_name' => $data[3],'scientific_name' => $data[4]};
        push @{ $target_taxons{ $back_gi_to_taxinfo{$data[1]}->{'taxid'} }->{'seqs'} }, $back_gi_to_taxinfo{$data[1]}->{'accession' };
        if (not exists $seen_taxon{ $back_gi_to_taxinfo{$data[1]}->{'taxid'} } ){
            $seen_taxon{$back_gi_to_taxinfo{$data[1]}->{'taxid'}} =  $taxon_counter++; 
        }
        $target_taxons{ $back_gi_to_taxinfo{$data[1]}->{'taxid'} }->{'tax_number'} = $seen_taxon{$back_gi_to_taxinfo{$data[1]}->{'taxid'}};
    }
    unlink("$tmp_file.txt");
    unlink("$tmp_file.csv");
    print_OUT("   '-> Temorary files removed");
    return(\%back_gi_to_taxinfo,\%target_taxons);
}

sub get_tree_from_taxids {
    my ($self, @species_ids) = @_;
    
    # the full lineages of the species are merged into a single tree
    my $tree;
    foreach my $ncbi_id (@species_ids) {
        if ($ncbi_id) {
            my $node = $self->get_taxon(-taxonid => $ncbi_id);
            #$node->name('supplied', $name);
            
            if ($tree) {
                $tree->merge_lineage($node);
            } else {
                $tree = Bio::Tree::Tree->new(-verbose => $self->verbose, -node => $node);
            }
        }
        else {
            $self->throw("No taxonomy database node for species ".$ncbi_id);
        }
    }
    
    return $tree;
}


1;
