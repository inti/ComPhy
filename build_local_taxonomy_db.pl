#!/usr/bin/perl 
use strict;
use warnings;
use MLDBM qw(DB_File Storable);
use Fcntl;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use LWP::Simple;
use Archive::Tar;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;

# local modules

use PhyloStratiphytUtils;

our ( $help, $man, $get_files, $tax_folder, $keep_files, $tax_ftp, $prot_only, $nucl_only, $p_bar);

GetOptions(
    'help' => \$help,
    'man' => \$man,
    'get_files' => \$get_files,
    'keep_files' => \$keep_files,
    'tax_folder=s' => \$tax_folder,
    'url=s' => \$tax_ftp,
    'nucl_only' => \$nucl_only,
    'prot_only' => \$prot_only,
    'progress_bar|bar' => \$p_bar
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $man);

$tax_ftp ||= "ftp://ftp.ncbi.nih.gov/pub/taxonomy";


my @files = qw (gi_taxid_prot.dmp.gz gi_taxid_nucl.dmp.gz taxdump.tar.gz);
if ($tax_folder){
    print_OUT("Setting taxonomy folder to [ $tax_folder ]");
    unless(-d $tax_folder) {
        mkdir($tax_folder);
    } 
    chdir("$tax_folder");
}

if (defined $get_files){
    print_OUT("Downloading files from [ ftp://ftp.ncbi.nih.gov/pub/taxonomy ]");
    foreach  my $file (@files) { 
        next if (defined $nucl_only and $file eq "gi_taxid_prot.dmp.gz");
        next if (defined $prot_only and $file eq "gi_taxid_nucl.dmp.gz");
        
        print_OUT("   '-> [ $file ]");
        getstore("$tax_ftp/$file", $file);
        if (($file =~ m/tgz$/) or ($file =~ m/tar.gz$/)) {
            my $tar = Archive::Tar->new;
            $tar->read($file);
            $tar->extract();
        } elsif (($file =~ m/\.gz$/) or ($file =~ m/.Z$/)) {
            my ($output) = $file;
            $output =~ s/\.gz$//;
            $output =~ s/\.Z$//;
            my $status = gunzip $file => $output or die "gunzip failed: $GunzipError\n";        
        }
    }
}

my $seq_id_files = 'both';
if ($nucl_only) {
    $seq_id_files = "nucl";        
} elsif ($prot_only) {
    $seq_id_files = "prot";        
}

&build_database("seq_to_gi",".",$seq_id_files);

unless($keep_files) {
    opendir(DIR, "./") or die $!;
    my @dir = readdir(DIR);
    
    foreach my $thing (@dir){
        next if (-d $thing);
        next if ($thing =~ m/\.db$/);
        unlink ($thing);
    }
}

print_OUT("Done");