############################ NOTE  ########
# This routines were obtained from the ebot software output as a way to facilitate
# the use of the NCBI DBs. They are not my original work but I have made some modification to them.
###################
#This script contains the routines of the NCBI_PowerScripting.pm module used in the course
# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
#
# Author:  Eric W. Sayers  sayers@ncbi.nlm.nih.gov
# http://www.ncbi.nlm.nih.gov/Class/PowerTools/eutils/course.html
#
#
# ---------------------------------------------------------------------------


#Contains the following subroutines:
#read_params
#egquery
#esearch
#esearch_links
#esummary
#esummary_links_by_id
#efetch
#efetch_batch
#efetch_links_by_id
#elink
#elink_history
#elink_batch
#elink_batch_to
#elink_by_id
#elink_by_id_to
#elink_out
#epost_uids
#epost_file
#epost_set
#print_summary
#print_links
#print_link_summaries
#get_uids
#read_index
#get_linknames
#get_link_report
#extract_links
#get_ftp_file


package NCBI_PowerScripting;
#use strict;
use warnings;
use Carp;
use Exporter qw (import);
use LWP::Simple;
use LWP::UserAgent;
use Net::FTP;


my $delay = 0;
my $maxdelay = 3;
my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";

####
our (@EXPORT, @EXPORT_OK, %EXPORT_TAGS);

@EXPORT = qw(read_params egquery esearch esearch_links esummary efetch efetch_batch elink elink_history elink_batch elink_batch_to elink_by_id elink_by_id_to elink_out epost_uids epost_file epost_set print_summary print_links print_link_summaries get_uids read_index get_linknames get_link_report extract_links get_ftp_file);				# symbols to export by default
@EXPORT_OK = qw(read_params egquery esearch esearch_links esummary efetch efetch_batch elink elink_history elink_batch elink_batch_to elink_by_id elink_by_id_to elink_out epost_uids epost_file epost_set print_summary print_links print_link_summaries get_uids read_index get_linknames get_link_report extract_links get_ftp_file);			# symbols to export on request

####
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
#************************************************************
#** BEGIN NCBI_PowerScripting MODULE ROUTINES ***************
#*************************************************************

sub read_params {
    
    # Reads input parameters from file supplied on command line
    # Input file must have lines of the following format:
    #   parameter|value
    # where parameter is the URL parameter name and value is the
    # value to be assigned to parameter
    # For ELink, the parameter "dbfrom" must be on a line before
    # the id parameters. This allows multiple &id parameters
    # Input: file named on command line
    # Output: %params; keys are parameter names, values are values
    # Example: $params{db} = 'nucleotide'
    # $params{id} is an array if "dbfrom" parameter is in input file
    
    my ($param, $value);
    my (@keys, @test);
    my %params;
    my %mark;
    my $dbfrom;
    
    #check for correct command line syntax
    if ($#ARGV != 0) { die "Usage: [eutil].pl input_file\n"; }
    
    #read input parameter file
    open(INPUT, "<$ARGV[0]") || die "Aborting. Can't open $ARGV[0]\n";
    
    while (<INPUT>) {
        
        chomp;
        ($param, $value) = split(/\|/);
        if ($param eq 'dbfrom') { $dbfrom = 1; }
        if (($param eq 'id') && ($dbfrom)) {
            push (@{$params{$param}}, $value);
        }
        else {
            $params{$param} = $value;
        }
    }
    
    close INPUT;
    
    return (%params);
    
}

#************************************************************************

sub egquery {
    
    # Performs EGQuery.
    # Input: %params:
    # $params{term} - Entrez query
    # $params{tool} - tool name
    # $params{email} - e-mail address
    # Output = %results; keys are databases, values are UID counts
    
    my %params = @_;
    my ($url, $raw);
    my @out;
    my $database;
    my %results;
    my ($begin, $end);
    
    sleep($delay);
    
    $url = $base . "egquery.fcgi?term=$params{term}";
    $url .= "&tool=$params{tool}&email=$params{email}";
    
    print "\n$url\n\n" if ($params{verbose});
    
    $begin = time;
    $raw = get($url);
    
    @out = split(/^/, $raw);
    
    foreach (@out) {
        
        if (/<DbName>(.*)<\/DbName>/) { $database = $1; }
        if (/<Count>(\d+)<\/Count>/) { $results{$database} = $1; }
        
    }
    
    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    return(%results);
    
}

#*********************************************************************

sub esearch {
    
    # Performs ESearch.
    # Input: %params
    # $params{db} - database
    # $params{term} - Entrez query
    # $params{usehistory} (y/n) - flag for using the Entrez history server, default = y
    # $params{retstart} - first item in results list to display (default = 0)
    # $params{retmax} - number of items in results list to display (default = 20)
    # $params{WebEnv} - Web Environment for accessing existing data sets
    # $params{reldate} - relative date, days preceding current date
    # $params{mindate} - begin date of range
    # $params{maxdate} - end date of range
    # $params{datetype} - type of date limited by relate, mindate, maxdate (ie edat, mdat, pdat, cdat)
    # $params{sort} - sort key
    # $params{tool} - tool name
    # $params{email} - e-mail address
    # $params{verbose} - (y/n) - causes messages to be sent to STDOUT; default = y
    # $params{http} - 'get' - uses HTTP Get; otherwise uses HTTP Post
    #
    # Output: %results: keys are 'db', 'count', 'query_key', 'WebEnv', 'uids'
    # $results{uids} is an array
    
    my %params = @_;
    my ($url, $url_params, $raw, $raw_cont);
    my @out;
    my %results;
    my ($begin, $end);
    my @options = qw(usehistory WebEnv retstart retmax reldate mindate maxdate datetype sort tool email);
    
    $params{verbose} = 'y' unless ($params{verbose});
    
    sleep($delay);
    
    $params{usehistory} = 'y';
    $params{http} = '' unless ($params{http});
    if ( ($params{db}) && ($params{term}) ) {
        $url_params = "db=$params{db}&term=$params{term}";
    }
    else {
        print "\nWARNING: ESearch requires both &db and &term!\n\n";
    }
    
    foreach my $opt (@options) {
        $url_params .= "&$opt=$params{$opt}" if ($params{$opt});
    }
    
    if ($params{http} eq 'get') {
        # use HTTP Get
        $url = $base . "esearch.fcgi?$url_params";
        $raw_cont = get($url);
    }
    else {
        # use HTTP Post
        $url = $base . "esearch.fcgi";
        
        #create user agent
        my $ua = new LWP::UserAgent;
        $ua->agent("esearch/1.0 " . $ua->agent);
        
        #create HTTP request object
        my $req = new HTTP::Request POST => "$url";
        $req->content_type('application/x-www-form-urlencoded');
        $req->content("$url_params");
        
        $begin = time;
        #post the HTTP request
        $raw = $ua->request($req);
        
        #print "\n$url?$url_params\n\n" if ($params{verbose} eq 'y');
        
        $raw_cont = $raw->content;
    }
    
    $raw_cont =~ /<Count>(\d+)<\/Count>/s;
    $results{count} = $1;
    $raw_cont =~ /<QueryKey>(\d+)<\/QueryKey>.*<WebEnv>(\S+)<\/WebEnv>/s;
    $results{query_key} = $1 if ($params{usehistory} eq 'y');
    $results{WebEnv} = $2;
    $results{db} = $params{db};
    $results{usehistory} = 'y' if ($params{usehistory} eq 'y');
    @out = split(/^/, $raw_cont);
    
    foreach (@out) {
        if (/<Id>(\d+)<\/Id>/) { push (@{$results{uids}}, $1); }
    }
    
    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    if ( ($results{count} == 0) && ($params{verbose} eq 'y') ) {
        print "ALERT: ESearch found no records for this query:\n";
        print "$params{term}\n";
    }
    
    $results{email} = $params{email};
    $results{tool} = $params{tool};
    
    return(%results);
    
}

#****************************************************************

sub esearch_links {
    
    # Performs ESearch on the output of elink_by_id ONLY
    # Input: %params:
    # $params{db} - database
    # $params{term} - Entrez query, where # represents the query key from elink
    # $params{WebEnv} - Web Environment for input data set
    # $params{query_key} - query key for input data set
    # $params{linkfile} - index file (.idx) produced by elink_by_id
    # $params{infile} - same as linkfile to provide backward compatibility
    # $params{outfile} - index file (.idx) containing results of esearch for each UID
    # 	default = infile_search.idx
    # $params{tool} - tool name
    # $params{email} - e-mail address
    #
    # Output: one hash and one file:
    # %results: 'query_key', 'WebEnv', 'linkfile'
    #   query_key and WEbEnv point to the whole set of limited UIDs
    #   $results{linkfile} = name of output index file
    # outfile.idx - index file containing lines of the form
    #   input UID in dbfrom:linked UIDs in db (comma-delimited list)
    
    my %params = @_;
    my (%results, %initial, %links, %output, %mark, %pparams);
    my (@uids, @init, @filt, @out, @diff);
    my ($uid, $final, $file, $uidlist);
    my @options = qw(tool email);
    
    # Run ESearch to get the set matching the limiting $term
    
    unless ( ($params{db}) && ($params{term}) && ($params{query_key}) && ($params{WebEnv}) ) {
        print "\nWARNING: esearch_links requires db, term, query_key and WebEnv in input hash!\n\n";
    }
    
    unless ( ($params{linkfile}) || ($params{infile}) ) {
        print "WARNING: esearch_links requires linkfile in input hash!\n\n";
    }
    
    foreach my $opt (@options) {
        $params{$opt} = '' unless ($params{$opt});
    }
    
    $params{term} =~ s/#/%23$params{query_key}/;
    
    $params{linkfile} = $params{infile} if ( ($params{infile}) && (!$params{linkfile}) );
    
    %results = esearch(%params);
    
    @uids = get_uids(%results);
    
    # Read input index file from elink_by_id
    
    %initial = read_index($params{linkfile});
    
    # Write new index file
    
    if ($params{outfile}) {
        $file = $params{outfile};
    }
    else {
        if ($params{linkfile} =~ /(\S+)\..*$/) { $file = $1; }
        else { $file = $params{linkfile} };
        
        $file .= '_search.idx';
    }
    
    open (OUTPUT, ">$file") || die "Can't open $file!\n";
    
    @links{@uids} = ();
    
    foreach $uid (keys %initial) {
        
        undef @filt;
        @init = split(/,/, $initial{$uid});
        
        foreach (@init) {
            push(@filt, $_) if exists $links{$_};
        }
        
        $final = join(',', @filt);
        print OUTPUT "$uid:$final\n";
        
        # make the final UID list non-redundant
        grep($mark{$_}++, @out);
        @diff = grep(!$mark{$_}, @filt);
        
        @out = (@out, @diff);
        
    }
    
    close OUTPUT;
    
    # post set of UIDs
    
    $uidlist = join(',', @out);
    
    $pparams{db} = $params{db};
    $pparams{id} = $uidlist;
    
    %results = epost_set(%pparams);
    
    print "Wrote index file to $file.\n";
    
    $results{linkfile} = $file;
    $results{email} = $params{email};
    $results{tool} = $params{tool};
    
    return %results;
    
}

#****************************************************************

sub esummary {
    
    # Performs ESummary.
    # Input: %params:
    # $params{db} - database
    # $params{id} - UID list (ignored if query_key exists)
    # $params{query_key} - query_key
    # $params{WebEnv} - web environment
    # $params{retstart} - first DocSum to retrieve
    # $params{retmax} - number of DocSums to retrieve
    # $params{outfile} - name of output file for XML (default = docsums)
    # $params{batch} - size of batch to retrieve (default = 50,000)
    # $params{tool} - tool name
    # $params{email} - e-mail address
    #
    # Output: XML written to $params{outfile}
    
    my %params = @_;
    my ($url, $raw, $done, $trial, $url_params);
    my @out;
    my $id;
    my %results;
    my ($begin, $end, $retmax, $tempfile, $limit, $expect, $c);
    my $batch = 10000;
    my @options = qw(retstart retmax tool email);
    
    sleep($delay);
    
    $params{outfile} = 'docsums' unless ($params{outfile});
    $batch = $params{batch} if ($params{batch});
    unless ( ($params{db}) && ( ($params{id}) || ( ($params{query_key}) && ($params{WebEnv}) ) ) ) {
        print "\nWARNING: ESummary requires &db and either &id or (&WebEnv, &query_key)!\n\n";
    }
    foreach my $opt (@options) {
        $params{$opt} = '' unless ($params{$opt});
    }
    
    #first use ESearch to determine the size of the dataset
    
    if ($params{retmax}) {
        $count = $params{retmax};
    }
    elsif ($params{num}) {
        $count = $params{num};
    }
    elsif ($params{query_key}) {
        $params{term} = "%23" . "$params{query_key}";
        $params{usehistory} = 'y';
        
        %results = esearch(%params);
        
        $count = $results{count};
    }
    else {
        @out = split(/,/, $params{id});
        $count = @out;
    }
    
    print "Retrieving ";
    print "no more than " if ($params{retmax});
    print "$count DocSums from $params{db}...\n";
    
    my $ua = new LWP::UserAgent;
    $ua->agent("esummary/1.0 " . $ua->agent);
    
    $tempfile = $params{outfile} . "_temp";
    
    $url = $base . "esummary.fcgi?db=$params{db}";
    
    if ($params{query_key}) {
        $url .= "&query_key=$params{query_key}&WebEnv=$params{WebEnv}&tool=$params{tool}&email=$params{email}";
    }
    else {
        $url = $base . "esummary.fcgi";
        $url_params = "db=$params{db}&id=$params{id}&tool=$params{tool}&email=$params{email}";
    }
    
    print "\n$url\n\n" if ($params{verbose});
    
    $begin = time;
    
    if ($params{query_key}) {
        # This is the default routine that returns the XML Esummary output
        # set retstart/retmax to input values if they exist
        # Otherwise loop in batches to get entire set
        if ( ($params{retstart}) || ($params{retmax}) ) {
            $url .= "&retstart=$params{retstart}&retmax=$params{retmax}";
            $raw = $ua->get($url, ':content_file' => $params{outfile});
        }
        else {
            open (OUTFILE, ">$params{outfile}");
            close OUTFILE;
            open (OUTFILE, ">>$params{outfile}");
            for (my $retstart=1; $retstart <= $count; $retstart+=$batch) {
                if ($retstart+$batch > $count) { $limit = $count; }
                else { $limit = $retstart + $batch - 1; }
                $expect = $limit - $retstart + 1;
                $done = 0; $trial = 0;
                until ($done) {
                    print "Retrieving records $retstart - $limit...";
                    
                    $url = $base . "esummary.fcgi?db=$params{db}";
                    $url .= "&query_key=$params{query_key}&WebEnv=$params{WebEnv}";
                    $url .= "&retstart=$retstart&retmax=$batch";
                    
                    $raw = $ua->get($url, ':content_file' => $tempfile);
                    
                    #count the number of DocSums retrieved
                    $c = 0;
                    open(TEMP, $tempfile);
                    while (<TEMP>) { $c++ if /<DocSum>/; }
                    close TEMP;
                    
                    print "$c records confirmed.";
                    if ($c >= $expect) {
                        $done = 1;
                        print "\n";
                    }
                    else {
                        $trial++;
                        if ($trial == 3) {
                            $done = 1;
                            print "Giving up!\n";
                        }
                        else {
                            print " Trying again...\n";
                        }
                    }
                }
                
                open (TEMP, $tempfile);
                while (<TEMP>) { print OUTFILE; }
                close TEMP;
                
            }
            close OUTFILE;
            
            unlink($tempfile);
            
        }
    }
    else {
        $ua = new LWP::UserAgent;
        $ua->agent("esummary/1.0 " . $ua->agent);
        
        #create HTTP request object
        my $req = new HTTP::Request POST => "$url";
        $req->content_type('application/x-www-form-urlencoded');
        $req->content("$url_params");
        
        #post the HTTP request
        $raw = $ua->request($req);
        
        #print "\n$url?$url_params\n\n" if ($params{verbose} eq 'y');
        
        open (OUTFILE, ">$params{outfile}");
        print OUTFILE $raw->content;
        close OUTFILE;
    }
    
    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    print "Document summaries written to $params{outfile}.\n";
    
}

#****************************************************************

sub efetch {
    
    # Performs EFetch.
    # Input: %params:
    # $params{db} - database
    # $params{id} - UID list (ignored if query_key exists)
    # $params{query_key} - query key
    # $params{WebEnv} - web environment
    # $params{retmode} - output data format
    # $params{rettype} - output data record type
    # $params{retstart} - first record in set to retrieve
    # $params{retmax} - number of records to retrieve
    # $params{seq_start} - retrieve sequence starting at this position
    # $params{seq_stop} - retrieve sequence until this position
    # $params{strand} - which DNA strand to retrieve (1=plus, 2=minus)
    # $params{complexity} - determines what data object to retrieve
    # $params{report} - report format for db=taxonomy and snp
    # $params{tool} - tool name
    # $params{email} - e-mail address
    #
    # Output: $raw; raw EFetch output
    
    my %params = @_;
    my ($url, $raw);
    my ($begin, $end);
    my @options = qw(retmode rettype retstart retmax seq_start seq_stop strand complexity report tool email verbose);
    
    sleep($delay);
    
    unless ( ($params{db}) && ( ($params{id}) || ( ($params{query_key}) && ($params{WebEnv}) ) ) ) {
        print_OUT("\nWARNING: EFetch requires &db and either &id or (&WebEnv, &query_key)!");
    }
    foreach my $opt (@options) {
        $params{$opt} = '' unless ($params{$opt});
    }
    
    $url = $base . "efetch.fcgi?db=$params{db}";
    
    if ($params{query_key}) {
        $url .= "&query_key=$params{query_key}&WebEnv=$params{WebEnv}";
    }
    else {
        $url .= "&id=$params{id}";
    }
    
    if (($params{report}) && ($params{db} eq 'taxonomy') ) {
        $url .= "&report=$params{report}&retmode=text";
    }
    elsif ( ($params{report}) && ($params{db} eq 'snp') ) {
        $url .= "&report=$params{report}&retmode=$params{retmode}";
    }
    else {
        $url .= "&retmode=$params{retmode}&rettype=$params{rettype}";
    }
    
    $url .= "&retstart=$params{retstart}&retmax=$params{retmax}";
    $url .= "&seq_start=$params{seq_start}&seq_stop=$params{seq_stop}";
    $url .= "&strand=$params{strand}&complexity=$params{complexity}";
    
    $url .= "&tool=$params{tool}&email=$params{email}";
    
    print "\n$url\n\n" if ($params{verbose});
    
    $begin = time;
    $raw = get($url);
    
    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    return($raw);
    
}

#****************************************************************

sub efetch_batch {
    
    # Uses efetch to download a large data set in 500 record batches
    # The data set must be stored on the History server
    # The output is sent to a file named $params{outfile}
    # Input: %params:
    # $params{db} - link to database
    # $params{query_key} - query key
    # $params{WebEnv} - web environment
    # $params{retmode} - output data format
    # $params{rettype} - output data record type
    # $params{seq_start} - retrieve sequence starting at this position
    # $params{seq_stop} - retrieve sequence until this position
    # $params{strand} - which DNA strand to retrieve (1=plus, 2=minus)
    # $params{complexity} - determines what data object to retrieve
    # $params{tool} - tool name
    # $params{email} - e-mail address
    # $params{outfile} - name of output file
    # $params{report} - report format for db=taxonomy
    # $params{batch} - number of records retreived with each efetch URL (default = 500)
    #     setting 'batch' = -1 sets $retmax to null
    #
    # Output: nothing returned; raw EFetch output sent to $params{outfile}
    #   default file name - fetch.out
    # Other output: periodic status messages sent to standard output
    
    my %params = @_;
    my ($url, $raw);
    my ($begin, $end);
    my %results;
    my ($count, $first, $last);
    my ($retstart, $retmax);
    my @options = qw(retmode rettype seq_start seq_stop strand complexity report tool email);
    
    unless ( ($params{db}) && ( ($params{id}) || ( ($params{query_key}) && ($params{WebEnv}) ) ) ) {
        print_OUT("\nWARNING: EFetch requires &db and either &id or (&WebEnv, &query_key)!");
    }
    foreach my $opt (@options) {
        $params{$opt} = '' unless ($params{$opt});
    }
    
    $params{outfile} = 'fetch.out' unless ($params{outfile});
    $params{batch} = 0 unless ($params{batch});
    
    if ($params{batch} == -1) {
        $retmax = '';
    }
    elsif ($params{batch}) {
        $retmax = $params{batch};
    }
    else {
        $retmax = 500;
    }
    
    if ($params{retmax}) {
        $count = $params{retmax};
        $retmax = $count if ($count < $retmax);
    }
    
    #first use ESearch to determine the size of the dataset
    
    if ($params{num}) {
        $count = $params{num};
    }
    else {
        $params{term} = "%23" . "$params{query_key}";
        $params{usehistory} = 'y';
        
        %results = esearch(%params);
        
        $count = $results{count} unless ($count);
    }
    
    $params{retmax} = $retmax;
    
    print_OUT("Retrieving $count records from $params{db}...");
    
    open (OUT, ">$params{outfile}") || die "Aborting. Can't open $params{outfile}\n";
    
    if ($retmax) {
        # retrieve the data set in batches
        for ($retstart = 0; $retstart < $count; $retstart += $retmax) {
            
            sleep($delay);
            $params{retstart} = $retstart;
            $begin = time;
            $raw = efetch(%params);
            
            print OUT $raw;
            
            if ($retstart + $retmax > $count) { $last = $count; }
            else { $last = $retstart + $retmax; }
            $first = $retstart + 1;
            
            print_OUT("Received records $first - $last.");
            $end = time;
            $delay = $maxdelay - ($end - $begin);
            if ($delay < 0) { $delay = 0; }
        }
    }
    else {
        # retrieve the data set in one URL
        sleep($delay);
        $begin = time;
        $raw = efetch(%params);
        print OUT $raw;
        print_OUT("Received records 1 - $count.");
        $end = time;
        $delay = $maxdelay - ($end - $begin);
        if ($delay < 0) { $delay = 0; }
        
    }
    
    close OUT;
    
    #print_OUT("Wrote data to $params{outfile}.");
    
}

#****************************************************************

sub elink {
    
    # Performs ELink.
    # Input: %params:
    # $params{dbfrom} - link from database
    # $params{db} - link to database
    # $params{id} - array of UID lists (ignored if query_key exists)
    # $params{query_key} - query key
    # $params{WebEnv} - web environment)
    # $params{term} - Entrez term used to limit link results
    # $params{linkname} - Linkname to be returned (optional)
    # $params{tool} - tool name
    # $params{email} - e-mail address
    #
    # Output: %links:
    # @{$links{from}{$set}} = array of input UIDs in set $set
    # @{$links{to}{$linkname}{$set}} = array of linked UIDs in $db in set $set
    # @{$links{score}{$linkname}{$set}} = array of similarity scores for linked UIDs in $db in set $set
    # $links{db}{$linkname} = db name for $linkname links
    # where $set = integer corresponding to one &id parameter
    # value in the ELink URL
    
    my %params = @_;
    my ($url, $url_params, $raw, $out);
    my ($line, $getdata, $getid, $link, $id, $set, $name);
    my @out;
    my @link_ids;
    my $ids;
    my %results;
    my $db;
    my ($begin, $end);
    my $giveup = 3;
    my $ua;
    my @options = qw(term linkname tool email);
    
    sleep($delay);
    
    $set = 0;
    
    $url = $base . "elink.fcgi";
    
    unless ( ($params{dbfrom}) && ($params{db}) ) {
        print "\nWARNING: ELink requires both &dbfrom and &db!\n\n";
    }
    unless ( ($params{id}) || ( ($params{WebEnv}) && ($params{query_key}) ) ) {
        print "\nWARNING: ELink requires either &id or (&WebEnv, &query_key)!\n\n";
    }
    foreach my $opt (@options) {
        $params{$opt} = '' unless ($params{$opt});
    }
    
    if ($params{query_key}) {
        # use HTTP Get
        $url .= "?dbfrom=$params{dbfrom}&db=$params{db}&term=$params{term}";
        $url .= "&query_key=$params{query_key}&WebEnv=$params{WebEnv}";
        $url .= "&linkname=$params{linkname}&tool=$params{tool}&email=$params{email}";
        
        print "\n$url\n\n" if ($params{verbose});
        
        $begin = time;
        $trial = 0;
        $failure = 1;
        while (($failure) && ($trial < $giveup)) {
            $raw = get($url);
            if ( ($raw =~ /ERROR/) || ($raw =~ /Error/) || ($raw !~ /<DbFrom>/) ) {
                #    print "Links failed. Trying again...\n";
            }
            else { $failure = 0; }
            $trial++;
        }
        $end = time;
        print "Links failed after $giveup trials. Giving up!\n" if ($failure);
        
    }
    else {
        # use HTTP Post
        $url_params = "dbfrom=$params{dbfrom}&db=$params{db}&term=$params{term}";
        $url_params .= "&linkname=$params{linkname}&tool=$params{tool}&email=$params{email}";
        foreach $ids (@{$params{id}}) {
            $url_params .= "&id=$ids";
        }
        $ua = new LWP::UserAgent;
        $ua->agent("elink/1.0 " . $ua->agent);
        
        #create HTTP request object
        my $req = new HTTP::Request POST => "$url";
        $req->content_type('application/x-www-form-urlencoded');
        $req->content("$url_params");
        
        #post the HTTP request
        $begin = time;
        $out = $ua->request($req);
        $end = time;
        $raw = $out->content;
    }
    
    
    # parse output XML
    @out = split(/^/,$raw);
    
    $getdata = 0;
    
    $set = 0;
    while ($raw =~ /<LinkSet>(.*?)<\/LinkSet>/sg) {
        
        $linkset = $1;
        if ($linkset =~ /<IdList>(.*?)<\/IdList>/sg) {
            $ids = $1;
            while ($ids =~ /<Id>(\d+)<\/Id>/sg) {
                push (@{$results{from}{$set}}, $1);
            }
        }
        
        while ($linkset =~ /<LinkSetDb>(.*?)<\/LinkSetDb>/sg) {
            $linksetdb = $1;
            if ($linksetdb =~ /<LinkName>(.*)<\/LinkName>/) {
                $linkname = $1;
                $results{db}{$linkname} = $1 if ($linksetdb =~ /<DbTo>(.*)<\/DbTo>/);
                while ($linksetdb =~ /<Link>(.*?)<\/Link>/sg) {
                    $link = $1;
                    push (@{$results{to}{$linkname}{$set}}, $1) if ($link =~ /<Id>(\d+)<\/Id>/);
                    push (@{$results{score}{$linkname}{$set}}, $1) if ($link =~ /<Score>(\d+)<\/Score>/);
                }
            }
        }
        $set++;
    }
    
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    return(%results);
    
}

#************************************************************

sub elink_history {
    
    # Uses ELink with &cmd=neighbor_history to post ELink results
    # on the History server
    
    # Input: %params:
    # $params{dbfrom} - link from database
    # $params{db} - link to database
    # $params{id} - array of UID lists (ignored if query_key exists)
    # $params{query_key} - query key
    # $params{WebEnv} - web environment
    # $params{term} - Entrez term used to limit link results
    # $params{linkname} - Linkname to be returned (optional)
    # $params{tool} - tool name
    # $params{email} - e-mail address
    #
    # Output: %links:
    # @{$links{from}{$set}} = array of input UIDs in set $set
    # $links{to}{$linkname}{$set}{query_key} = query_key of $linkname links in set $set
    # $links{db}{$linkname} = db name for $linkname links
    # $links{WebEnv} = Web Environment of linked UID sets
    # where $set = integer corresponding to one &id parameter
    # value in the ELink URL
    # NOTE: If no links are found, query_key will be set to -1
    
    my %params = @_;
    my ($url, $url_params, $raw);
    my ($line, $getdata, $getid, $link, $id, $set, $name);
    my @out;
    my @link_ids;
    my $ids;
    my %results;
    my $db;
    my ($begin, $end);
    my @options = qw(term linkname tool email);
    
    sleep($delay);
    
    unless ( ($params{dbfrom}) && ($params{db}) ) {
        print "\nWARNING: ELink requires both &dbfrom and &db!\n\n";
    }
    unless ( ($params{id}) || ( ($params{WebEnv}) && ($params{query_key}) ) ) {
        print "\nWARNING: ELink requires either &id or (&WebEnv, &query_key)!\n\n";
    }
    foreach my $opt (@options) {
        $params{$opt} = '' unless ($params{$opt});
    }
    
    $set = 0;
    
    $url = $base . "elink.fcgi";
    
    $url_params = "dbfrom=$params{dbfrom}&db=$params{db}";
    $url_params .= "&cmd=neighbor_history&term=$params{term}";
    
    if ($params{query_key}) {
        
        $url_params .= "&query_key=$params{query_key}&WebEnv=$params{WebEnv}";
        
    }
    else {
        
        foreach $ids (@{$params{id}}) {
            $url_params .= "&id=$ids";
        }
    }
    $url_params .= "&linkname=$params{linkname}&tool=$params{tool}&email=$params{email}";
    
    print "\n$url_params\n\n" if ($params{verbose});
    
    #create user agent
    my $ua = new LWP::UserAgent;
    $ua->agent("elink/1.0 " . $ua->agent);
    
    #create HTTP request object
    my $req = new HTTP::Request POST => "$url";
    $req->content_type('application/x-www-form-urlencoded');
    $req->content("$url_params");
    
    $begin = time;
    #post the HTTP request
    $raw = $ua->request($req);
    
    #print "\n$url?$url_params\n\n" if ($params{verbose} eq 'y');
    
    $raw_cont = $raw->content;
    #print $raw_cont;
    #$begin = time;
    #$raw = get($url);
    
    #parse XML output
    $set = 0;
    while ($raw_cont =~ /<LinkSet>(.*?)<\/LinkSet>/sg) {
        
        $linkset = $1;
        if ($linkset =~ /<IdList>(.*?)<\/IdList>/sg) {
            $ids = $1;
            while ($ids =~ /<Id>(\d+)<\/Id>/sg) {
                push (@{$results{from}{$set}}, $1);
            }
        }
        
        while ($linkset =~ /<LinkSetDbHistory>(.*?)<\/LinkSetDbHistory>/sg) {
            $linksetdb = $1;
            if ($linksetdb =~ /<LinkName>(.*)<\/LinkName>/) {
                $linkname = $1;
                $results{to}{$linkname}{$set}{query_key} = $1 if ($linksetdb =~ /<QueryKey>(.*)<\/QueryKey>/);
                $results{db}{$linkname} = $1 if ($linksetdb =~ /<DbTo>(.*)<\/DbTo>/);
            }
        }
        $results{WebEnv} = $1 if ($linkset =~ /<WebEnv>(.*)<\/WebEnv>/);
        $set++;
    }
    
    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    return(%results);
    
}

#********************************************************************

sub elink_batch {
    
    # Produces links for a single set of records posted on the history server
    # from dbfrom to db. The routine segments the set in batches of size $batch
    # and then produces a non-redundant set of links for the entire set.
    
    #input hash: {WebEnv} = web environment of input set
    #	     {query_key} = query key of input set
    #	     {id} = list of UIDs (ignored if query_key exists)
    #	     {dbfrom} = database of input set, source db
    #	     {db} = destination db for elink
    #	     {term} = term parameter for elink
    #	     {linkname} = name of desired link; if set, output hash is one-dimensional
    #	     {http} - 'get' - HTTP Get; otherwise HTTP Post
    # 	     {tool} - tool name
    # 	     {email} - e-mail address
    
    #output: if {linkname} is NOT set:
    #        %links{$linkname}{query_key} - query key for unique linkname links
    #              {$linkname}{WebEnv} - web environment for unique linkname links
    #	       {$linkname}{db} - database containing $linkname links
    #
    #	 if (linkname} is set:
    #	 %links{query_key} - query key for links
    #              {WebEnv} - web environment for links
    #	       {db} - database containing links
    
    my %params = @_;
    
    my $batch = 10000;
    my $giveup = 3;
    my ($retstart, $first, $last, $max, $trial, $failure, $name, $cur);
    my (%sparams, %lparams, %ct, %pparams, %posted, %iparams);
    my (%sresults, %lresults, %presults);
    my (%links, %count, %foundnames);
    my %output;
    my @new;
    my @options = qw(term linkname http tool email);
    
    unless ( ($params{dbfrom}) && ($params{db}) ) {
        print "\nWARNING: ELink requires both &dbfrom and &db!\n\n";
    }
    unless ( ($params{id}) || ( ($params{WebEnv}) && ($params{query_key}) ) ) {
        print "\nWARNING: ELink requires either &id or (&WebEnv, &query_key)!\n\n";
    }
    foreach my $opt (@options) {
        $params{$opt} = '' unless ($params{$opt});
    }
    
    if ( ($params{id}) && (!$params{query_key}) ) {
        # input data set is NOT on the history, so put it there with epost
        if ($params{id}[0]) {
            $params{id} = join(',',@{$params{id}});
        }
        
        $iparams{db} = $params{dbfrom};
        $iparams{id} = $params{id};
        $iparams{http} = $params{http};
        
        %iparams = epost_set(%iparams);
        
        $params{query_key} = $iparams{query_key};
        $params{WebEnv} = $iparams{WebEnv};
        
    }
    
    # use smaller batches for computational neighbors
    $batch = 50 if ($params{dbfrom} eq $params{db});
    
    # use esearch to determine the size of the input data set
    $sparams{db} = $params{dbfrom};
    $sparams{term} = "%23$params{query_key}";
    $sparams{retmax} = $batch;
    $sparams{WebEnv} = $params{WebEnv};
    $sparams{usehistory} = 'y';
    
    
    %sresults = esearch(%sparams);
    $max = $sresults{count};
    
    if ($params{linkname}) {
        print "Finding links named $params{linkname} for $max $params{dbfrom} records...\n";
    }
    else {
        print "Finding all links from $max $params{dbfrom} records to $params{db}...\n";
    }
    
    $lparams{dbfrom} = $params{dbfrom};
    $lparams{db} = $params{db};
    $lparams{term} = $params{term};
    $lparams{linkname} = $params{linkname};
    
    $hparams{usehistory} = 'y';
    $hparams{db} = $params{db};
    $pparams{db} = $params{db};
    
    #batch elink from dbfrom to db
    
    for ($retstart=0; $retstart < $max; $retstart += $batch) {
        
        if (($retstart + $batch) > $max) {
            $last = $max;
        }
        else {
            $last = $retstart + $batch;
        }
        
        $first = $retstart + 1;
        
        # use esearch to retrieve each batch of input UIDs, and then use elink_history
        $sparams{retstart} = $retstart;
        
        %sresults = esearch(%sparams);
        
        $lparams{id}[0] = join(',', @{$sresults{uids}});
        
        $trial = 0;
        $failure = 1;
        while (($failure) && ($trial < $giveup)) {
            %lresults = elink_history(%lparams);
            foreach $name (keys %{$lresults{to}}) {
                if ($lresults{to}{$name}{0}{query_key} > 0) {
                    $failure = 0;
                }
                elsif ($lresults{to}{$name}{0}{query_key} == -1) {
                    print "No links found with name $name.\n";
                    $failure = 0;
                }
            }
            if ($failure) {
                if ( $trial < $giveup - 1 ) {
                    #         print "Links failed. Trying again...\n";
                }
                else {
                    print "No links found for records $first - $last.\n";
                }
            }
            $trial++;
        }
        
        foreach $key (keys %{$lresults{to}}) {
            $foundnames{$key}{ct}++;
            $foundnames{$key}{db} = $lresults{db}{$key};
        }
        
        # for each linkname found, add the new links for the current batch to the links found
        # for previous batches and then non-redundify this set
        # store the resulting set in @{$links{$linkname}}
        unless ($failure) {
            foreach $name (keys %foundnames) {
                undef %ct;
                if (exists($links{$name})) { $cur = @{$links{$name}} }
                else { $cur = 0; }
                if ($cur > 0) {
                    foreach (@{$links{$name}}) {
                        $ct{$_}++;
                    }
                }
                
                if ($lresults{to}{$name}{0}{query_key}) {
                    $pparams{db} =  $lresults{db}{$name};
                    $pparams{query_key} = $lresults{to}{$name}{0}{query_key};
                    $pparams{WebEnv} = $lresults{WebEnv};
                    
                    @new = get_uids(%pparams);
                    
                    foreach (@new) {
                        $ct{$_}++;
                    }
                }
                
                @{$links{$name}} = keys %ct;
                $count{$name} = @{$links{$name}};
            }
            
            print "Links complete for records $first - $last.\n";
            foreach (keys %count) {
                print "So far, $count{$_} unique links for $_.\n";
            }
        }
    } #end of batch loop
    
    foreach $name (keys %links) {
        
        $pparams{id} = join(',', @{$links{$name}});
        $pparams{db} = $foundnames{$name}{db};
        $pparams{http} = $params{http};
        
        %{$output{$name}} = epost_set(%pparams);
        $output{$name}{email} = $params{email};
        $output{$name}{tool} = $params{tool};
        
    }
    
    %output = extract_links($params{linkname}, %output) if ($params{linkname});
    
    return %output;
    
}

#*********************************************************************

sub elink_batch_to {
    
    # Runs elink_batch using a universal hash, with the destination db defined
    # by the $dbto parameter
    
    my ($dbto, %params) = @_;
    
    $params{dbfrom} = $params{db};
    $params{db} = $dbto;
    
    %params = elink_batch(%params);
    
    return(%params);
    
}

#*********************************************************************

sub elink_by_id {
    
    # Produces links for each member of a set of records posted on the history server
    # from dbfrom to db. The routine segments the set in batches of size $batch
    # and then produces a set of links for each UID in the set and places these on the
    # history using elink_history.
    
    #input hash: {WebEnv} = web environment of input set
    #	     {query_key} = query key of input set
    #	     {id} = list of UIDs (ignored if query_key exists)
    #	     {dbfrom} = database of input set, source db
    #	     {db} = destination db for elink
    #	     {term} = term parameter for elink
    #	     {linkname} = name of desired link; if set, output hash is one-dimensional
    #	     {outfile} = name of index file if get_uids is not 'n'
    #			  default = dbfrom_db.idx
    #	     {scorefile} = name of index file containing similarity scores if db=dbfrom
    #			  default = db.sco
    #	     {http} - 'get' - use HTTP Get; otherwise use HTTP Post
    
    #output: one hash and one file
    # if {linkname} is NOT set:
    # %links:
    # $links{$linkname}{query_key}
    # $links{$linkname}{WebEnv}
    # $links{$linkname}{db}
    # $links{$linkname}{linkfile}
    # $links{$linkname}{scorefile}
    #
    # if {linkname} is set:
    # %links:
    # $links{query_key}
    # $links{WebEnv}
    # $links{db}
    # $links{linkfile}
    # $links{scorefile}
    # Output WebEnv and querykey point to a non-redundant list of UIDs in db linked to all UIDs in dbfrom
    # NOTE: If no links are found, query_key will be set to -1
    # outfile.idx - index file containing lines of the form
    #   input UID in dbfrom:linked UIDs in db (comma-delimited list)
    
    my %params = @_;
    
    my $batch = 5000;
    my $giveup = 3;
    my ($retstart, $max, $trial, $failure, $set, $first, $last);
    my ($uid, $file, $num, $expect, $scorefile, $getscore, $name);
    my (@input, @diff, @uids, @temp);
    my (%sparams, %lparams, %iparams, %iresults);
    my (%sresults, %links, %mark, %output, %foundnames);
    my @options = qw(term linkname http outfile scorefile tool email);
    
    unless ( ($params{dbfrom}) && ($params{db}) ) {
        print "\nWARNING: ELink requires both &dbfrom and &db!\n\n";
    }
    unless ( ($params{id}) || ( ($params{WebEnv}) && ($params{query_key}) ) ) {
        print "\nWARNING: ELink requires either &id or (&WebEnv, &query_key)!\n\n";
    }
    foreach my $opt (@options) {
        $params{$opt} = '' unless ($params{$opt});
    }
    
    if ( ($params{id}) && (!$params{query_key}) ) {
        
        if ($params{id}[0]) {
            $params{id} = join(',',@{$params{id}});
        }
        
        $iparams{db} = $params{dbfrom};
        $iparams{id} = $params{id};
        $iparams{http} = $params{http};
        
        %iresults = epost_set(%iparams);
        
        $params{query_key} = $iresults{query_key};
        $params{WebEnv} = $iresults{WebEnv};
        
    }
    
    if ($params{dbfrom} eq $params{db}) {
        $batch = 50;
        $getscore = 1;
    }
    
    $sparams{db} = $params{dbfrom};
    $sparams{query_key} = $params{query_key};
    $sparams{WebEnv} = $params{WebEnv};
    
    @input = get_uids(%sparams);
    
    $lparams{dbfrom} = $params{dbfrom};
    $lparams{db} = $params{db};
    $lparams{term} = $params{term};
    $lparams{linkname} = $params{linkname};
    
    $max = @input;
    
    if ($params{linkname}) {
        print "Finding links named $params{linkname} for $max records...\n";
    }
    else {
        print "Finding all links from $max $params{dbfrom} records to $params{db}...\n";
    }
    
    #batch elink from dbfrom to db
    
    for ($retstart=0; $retstart < $max; $retstart += $batch) {
        
        if (($retstart + $batch) > $max) {
            $last = $max;
            $expect = $last - $retstart;
        }
        else {
            $last = $retstart + $batch;
            $expect = $batch;
        }
        $first = $retstart + 1;
        
        @{$lparams{id}} = @input[$retstart..$last-1];
        
        
        
        #put UIDs into arrays
        $trial = 0;
        $num = 0;
        %lresults = elink(%lparams);
        
        foreach $key (keys %{$lresults{to}}) { $foundnames{$key}++; }
        
        foreach $key (sort keys %{$lresults{from}}) {
            
            $uid = $lresults{from}{$key}[0];
            # remove self-hit from list of computational links
            foreach $name (keys %{$lresults{to}}) {
                if ($getscore) {
                    
                    @temp = @{$lresults{to}{$name}{$key}};
                    $links{$uid}{$name} = join(',', @temp[1..$#temp] );
                    
                    @temp =  @{$lresults{score}{$name}{$key}};
                    $links{$uid}{score}{$name} = join(',', @temp[1..$#temp] );
                }
                elsif ($lresults{to}{$name}{$key}) {
                    $links{$uid}{$name} = join(',', @{$lresults{to}{$name}{$key}} );
                }
            }
        }
        
        @temp = keys %{$lresults{from}};
        $num = @temp;
        if ($num == $expect) {
            print "Links found for records $first - $last.\n";
        }
        else {
            print "WARNING: For records $first - $last, found links for $num out of $expect UIDs.\n";
        }
        
    }
    
    if ($params{get_uids} ne 'n') {
        foreach $name (keys %foundnames) {
            # write index file and combine UIDs
            if ($params{outfile}) { $file = $params{outfile} . "_$name.idx"; }
            else { $file = $name . '.idx'; }
            
            if ($params{scorefile}) { $scorefile = $params{scorefile} . "_$name.sco"; }
            else { $scorefile = $name . '.sco'; }
            
            open (OUTPUT, ">$file") || die "Can't open $file!\n";
            
            if ($getscore) {
                open (SCORES, ">$scorefile") || die "Can't open $scorefile!\n";
            }
            
            undef @uids;
            undef %mark;
            
            foreach $key (keys %links) {
                if ($links{$key}{$name}) {
                    print OUTPUT "$key:$links{$key}{$name}\n";
                    @temp = split(/,/, $links{$key}{$name});
                }
                else { @temp = (); }
                print SCORES "$key:$links{$key}{score}{$name}\n" if ( ($getscore) && ($links{$key}{score}{$name}) );
                
                # make the final UID list non-redundant
                grep($mark{$_}++, @uids);
                @diff = grep(!$mark{$_}, @temp);
                
                @uids = (@uids, @diff);
                
            }
            
            close OUTPUT;
            close SCORES if ($getscore);
            
            print "Wrote link index file to $file.\n";
            print "Wrote scores to $scorefile.\n" if ($getscore);
            
            
            # post set of UIDs
            
            $pparams{db} = $lresults{db}{$name};
            $uidlist = join(',', @uids);
            $pparams{id} = $uidlist;
            $pparams{http} = $params{http};
            
            %{$output{$name}} = epost_set(%pparams);
            $output{$name}{linkfile} = $file;
            $output{$name}{scorefile} = $scorefile if ($scorefile);
            $output{$name}{tool} = $params{tool};
            $output{$name}{email} = $params{email};
            
        }
    }
    
    %output = extract_links($params{linkname}, %output) if ($params{linkname});
    
    return %output;
    
}

#*********************************************************************

sub elink_by_id_to {
    
    # Runs elink_by_id using a universal hash, with the destination db defined
    # by the $dbto parameter
    # $params{get_uids} is forced to null
    #	 {http} - 'get' - use HTTP Get; otherwise use HTTP Post
    
    my ($dbto, %params) = @_;
    
    $params{dbfrom} = $params{db};
    $params{db} = $dbto;
    $params{get_uids} = '';
    
    %params = elink_by_id(%params);
    
    return(%params);
    
}

#*********************************************************************

sub elink_out {
    
    # Performs ELink to find LinkOut data for records on the input
    # history set. The routine uses only cmd=llinks, llinkslib, and prlinks.
    # The input UIDs are processed in batches of size $batch (default = 100)
    
    #input hash: {WebEnv} = web environment of input set
    #	     {query_key} = query key of input set
    #	     {id} = list of UIDs (ignored if query_key exists)
    #	     {db} = dbfrom for elink
    #	     {cmd} = cmd parameter for elink: llinks (default), llinkslib, or prlinks
    #	     {holding} = string to search for holding provider
    #	     {outfile} = name of output file to write (default = $db_linkout)
    #	     {html} = if 'y' will produce an HTML output file in addition to raw XML
    #No output is produced except the file
    
    my %params = @_;
    my (%iparams, %sparams, %lparams, %sresults);
    my ($linkouts, $linkurl, $initial, $final, $middle, $tempfile);
    my $batch = 500;
    
    my $ua = new LWP::UserAgent;
    $ua->agent("elink/1.0 " . $ua->agent);
    
    #sanity check input
    unless ( ($params{db}) && ( ($params{id}) || ( ($params{WebEnv}) && ($params{query_key}) ) ) ) {
        print "\nWARNING: ELink requires either &id or (&WebEnv, &query_key)!\n\n";
    }
    foreach my $opt (@options) {
        $params{$opt} = '' unless ($params{$opt});
    }
    
    unless ( ($params{cmd} eq 'llinks') || ($params{cmd} eq 'llinkslib') || ($params{cmd} eq 'prlinks') || ($params{cmd} eq '') ) {
        print "\nWARNING: Invalid Elink LinkOut command mode: $params{cmd}\n";
        print "Use either llinks, llinkslib, or prlinks.\n";
        print "For the moment, llinks will be used.\n\n";
        $params{cmd} = 'llinks';
    }
    
    if ( ($params{holding}) && ($params{cmd} eq 'prlinks') ) {
        print "\nWARNING: The holding parameter cannot be used with cmd=prlinks.\n";
        print "Continuing with cmd=llinks.\n\n";
        $params{cmd} = 'llinks';
    }
    
    $params{cmd} = 'llinks' unless ($params{cmd});
    $params{outfile} = "$params{db}_linkout" unless ($params{outfile});
    $tempfile = $params{outfile} . ".temp";
    $params{html} = '' unless ($params{html});
    
    if ( ($params{id}) && (!$params{query_key}) ) {
        # input data set is NOT on the history, so put it there with epost
        if ($params{id}[0]) {
            $params{id} = join(',',@{$params{id}});
        }
        
        $iparams{db} = $params{dbfrom};
        $iparams{id} = $params{id};
        
        %iparams = epost_set(%iparams);
        
        $params{query_key} = $iparams{query_key};
        $params{WebEnv} = $iparams{WebEnv};
        
    }
    
    # use esearch to determine the size of the input data set
    $sparams{db} = $params{db};
    $sparams{term} = "%23$params{query_key}";
    $sparams{retmax} = $batch;
    $sparams{WebEnv} = $params{WebEnv};
    $sparams{usehistory} = 'y';
    
    
    %sresults = esearch(%sparams);
    $max = $sresults{count};
    
    print "Finding LinkOut data for $max records in $params{db}...\n";
    
    $lparams{dbfrom} = $params{db};
    $lparams{cmd} = $params{cmd};
    $lparams{holding} = $params{holding};
    
    #batch elink from dbfrom
    
    open (OUTPUT, ">$params{outfile}") || die "Can't open $params{outfile}: $!\n";
    
    for ($retstart=0; $retstart < $max; $retstart += $batch) {
        
        $initial = $final = $middle = 0;
        if ($max > $batch) {
            if ($retstart == 0) { $initial = 1; }
            elsif ( ($retstart + $batch) >= $max ) { $final = 1; }
            else { $middle = 1; }
        }
        
        if (($retstart + $batch) > $max) {
            $last = $max;
        }
        else {
            $last = $retstart + $batch;
        }
        
        $first = $retstart + 1;
        
        # use esearch to retrieve each batch of input UIDs, and then use elink_history
        $sparams{retstart} = $retstart;
        
        %sresults = esearch(%sparams);
        
        $lparams{id} = join(',', @{$sresults{uids}});
        
        $linkurl = $base . "elink.fcgi?dbfrom=$lparams{dbfrom}&cmd=$lparams{cmd}&id=$lparams{id}&holding=$params{holding}";
        
        #post the HTTP request
        $raw = $ua->get($linkurl, ':content_file' => $tempfile);
        
        open (TEMP, "$tempfile");
        
        # merge the batched XML documents into one large XML document by removing opening/closing lines as appropriate
        while (<TEMP>) {
            if ( (/<\?xml version/) || (/<\!DOCTYPE/) || (/<eLinkResult>/) || (/<LinkSet>/) || (/<DbFrom>/) || (/<IdUrlList>/) ) {
                print OUTPUT unless ( ($middle) || ($final) );
            }
            elsif ( (/<\/LinkSet>/) || (/<\/eLinkResult>/) || (/IdUrlList>/) ) {
                print OUTPUT unless ( ($initial) || ($middle) );
            }
            else { print OUTPUT; }
        }
        
        print "Links complete for records $first - $last.\n";
        
        close TEMP;
        
    } #end of batch loop
    
    close OUTPUT;
    
    unlink $tempfile;
    
    # Optional HTML output
    if ($params{html} eq 'y') {
        my ($data, $id, $prov, $provname, $provurl, $provgif, $provset, $provlink, $cur);
        my @list;
        my $htmlfile = $params{outfile} . ".html";
        
        open (TEMP, ">$tempfile") || die "Can't open $tempfile: $!\n";
        
        $/ = "</IdUrlSet>";
        
        open (IN, "$params{outfile}") || die "Can't open $params{outfile}!: $!\n";
        
        while (<IN>) {
            $data = $_;
            $id = $1 if ($data =~ /<Id>(\d+)<\/Id>/);
            push (@list, $id) if ($data =~ /<ObjUrl>/);
            # print TEMP "$id\t";
            while ($data =~ /<ObjUrl>(.*?)<\/ObjUrl>/sg) {
                $provurl = $provgif = $provlink = '';
                $prov = $1;
                $provurl = $1 if ($prov =~ /<Url>(.*)<\/Url>/);
                $provgif = $1 if ($prov =~ /<IconUrl>(.*)<\/IconUrl>/);
                $provlink = $1 if ($prov =~ /<LinkName>(.*)<\/LinkName>/);
                $prov =~ /<Provider>(.*)<\/Provider>/s;
                $provset = $1;
                $provname = $1 if ($provset =~ /<Name>(.*)<\/Name>/);
                
                $provurl =~ s/\&amp;/\&/g;
                if ($provurl =~ /^\/entrez/) {
                    $provurl = 'http://www.ncbi.nlm.nih.gov' . $provurl;
                }
                
                print TEMP "$id\t$provname\t$provlink\t$provurl\t$provgif\n";
                
            }
        }
        close IN;
        close TEMP;
        
        $/ = "\n";
        
        #print HTML
        open (OUT, ">$htmlfile");
        
        print OUT "<html><head><title>LinkOuts</title></head>\n";
        print OUT "<body>\n";
        print OUT "<h3 id=\"top\">LinkOut Results</h3><p>\n";
        foreach (@list) {
            print OUT "<a href=\"#$_\">$_</a><br>\n";
        }
        
        open (IN, "$tempfile");
        while (<IN>) {
            chomp;
            ($id, $provname, $provlink, $provurl, $provgif) = split("\t", $_);
            if ($id ne $cur) {
                print OUT "<hr><h4 id=\"$id\">LinkOuts for UID ";
                print OUT "<a href=\"http://www.ncbi.nlm.nih.gov/sites/entrez?db=$params{db}&list_uids=$id&cmd=retrieve&dopt=DocSum\" target=\"entrez\">$id</a>&nbsp;&nbsp;";
                print OUT "<a href=\"#top\">Top</a></h4><p>\n";
                $cur = $id;
            }
            print OUT "<a target=\"linkout\" href=\"$provurl\">$provname ";
            print OUT "($provlink)" if ($provlink);
            print OUT "</a>&nbsp;";
            print OUT "<img src=\"$provgif\">" if ($provgif);
            print OUT "<br>\n";
        }
        
        print OUT "</body></html>";
        
        close OUT;
        close IN;
        unlink $tempfile;
        
    }
    
}

#*********************************************************************

sub epost_uids {
    
    # Performs EPost, placing UIDs in the URL.
    # Input: %params:
    # $params{db} - database
    # $params{id} - list of UIDs
    # $params{WebEnv} - Web environment for existing history sets
    # $params{tool} - tool name
    # $params{email} - e-mail address
    #
    #Output: %results: keys are 'WebEnv' and 'query_key'
    
    my %params = @_;
    my ($url, $raw);
    my ($begin, $end);
    
    sleep($delay);
    
    $url = $base . "epost.fcgi?db=$params{db}&id=$params{id}";
    $url .= "&WebEnv=$params{WebEnv}";
    $url .= "&tool=$params{tool}&email=$params{email}";
    
    print "\n$url\n\n" if ($params{verbose});
    
    $begin = time;
    $raw = get($url);
    
    $raw =~ /<QueryKey>(\d+)<\/QueryKey>.*<WebEnv>(\S+)<\/WebEnv>/s;
    $results{query_key} = $1;
    $results{WebEnv} = $2;
    # Section added by Inti Pedroso: 5 oct 2012
    # count number of results
    my %count_results = %params;
    $count_results{'WebEnv'}= $results{'WebEnv'};
    $count_results{'query_key'}= $results{'query_key'};
    $count_results{term} = "%23" . "$results{'query_key'}";
    $count_results{usehistory} = 'y';
    
    %count_results = esearch(%params);
    
    $results{count} = $count_results{count};
    # End section added by Inti Pedroso: 5 oct 2012

    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    return(%results);
    
}

#*********************************************************************

sub epost_file {
    
    # Performs EPost, accepts input from file.
    # Input file must have one UID per line.
    # Input: %params:
    # $params{db} - database
    # $params{id} - filename containing a list of UIDs
    # $params{WebEnv} - Web environment for existing history sets
    # $params{tool} - tool name
    # $params{email} - e-mail address
    #
    # Output: %results: keys are 'WebEnv' and 'query_key';
    #		    num - number of records in input file
    
    my %params = @_;
    my ($uids, $id);
    my @list;
    my ($begin, $end, $count);
    my (%results, %current);
    my @options = qw(WebEnv tool email);
    
    sleep($delay);
    
    unless ( ($params{db}) && ($params{id}) ) {
        print "\nWARNING: EPost requires both &db and &id!\n\n";
    }
    foreach my $opt (@options) {
        $params{$opt} = '' unless ($params{$opt});
    }
    
    #read input file of UIDs, one per line
    open (INPUT, "$params{id}") || die "Can't open $params{id}\n";
    
    while (<INPUT>) {
        
        if (/^(\d+)\r*\n*$/) {
            $id = $1;
            push (@list, $id);
        }
        else {
            print "ALERT: Found invalid uid in input file: $_\n";
        }
        
    }
    
    $params{id} = join (',', @list);
    
    $begin = time;
    
    %results = epost_set(%params);
    
    %current = %results;
    
    $count = @list;
    print "Posted $count records to $params{db}.\n";
    $results{num} = $count;
    
    $current{term} = "%23$results{query_key}";
    %current = esearch(%current);
    
    $current{count} = 0 unless ($current{count});
    
    print "ALERT: Only $current{count} records are current!\n" if ($count != $current{count});
    
    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    return(%results);
    
}

#***********************************************************

sub epost_set {
    
    # Uses EPost to post a set of UIDs using the POST method
    # Useful for large sets of UIDs not from a disk file
    # Accepts a comma-delimited list of UIDs in $params{id}
    # $params{WebEnv} - Web environment for existing history sets
    # $params{http} - 'get' - uses HTTP Get; otherwise uses HTTP Post
    # Output: $results{query_key}, $results{WebEnv}
    #	  $results{num} - number of records in input set
    
    my (%params) = @_;
    my ($url_params, $raw, $url, $raw_cont);
    my ($begin, $end);
    my %results;
    my @options = qw(WebEnv tool email http);
    
    unless ( ($params{db}) && ($params{id}) ) {
        print "\nWARNING: EPost requires both &db and &id!\n\n";
    }
    foreach my $opt (@options) {
        $params{$opt} = '' unless ($params{$opt});
    }
    
    $url_params = "db=$params{db}&id=$params{id}";
    $url_params .= "&WebEnv=$params{WebEnv}";
    $url_params .= "&tool=$params{tool}&email=$params{email}";
    
    $url = $base . "epost.fcgi";
    
    @list = split(/,/, $params{id});
    $len = @list;
    
    if ($params{http} eq 'get') {
        $url .= "?$url_params";
        $raw_cont = get($url);
    }
    else {
        #create user agent
        my $ua = new LWP::UserAgent;
        $ua->agent("epost_file/1.0 " . $ua->agent);
        
        #create HTTP request object
        my $req = new HTTP::Request POST => "$url";
        $req->content_type('application/x-www-form-urlencoded');
        $req->content("$url_params");
        
        $begin = time;
        #post the HTTP request
        $raw = $ua->request($req);
        
        $raw_cont = $raw->content;
    }
    
    $raw_cont =~ /<QueryKey>(\d+)<\/QueryKey>.*<WebEnv>(\S+)<\/WebEnv>/s;
    $results{query_key} = $1;
    $results{WebEnv} = $2;
    $results{db} = $params{db};
    $results{tool} = $params{tool};
    $results{email} = $params{email};
    $results{num} = $len;
    
    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    return (%results);
    
}

#***********************************************************

sub print_summary {
    
    # Input: %results output from sub esummary
    
    my %results = @_;
    my ($id, $count, $i);
    my (@a, @b, @c);
    
    $count = 0;
    
    if ($results{homologene} eq 'y') {
        
        @a = sort keys %results;
        @b = sort keys %{$results{$a[0]}};
        foreach (@b) {
            @c = @{$results{$a[0]}{$_}};
            $count = $#c if ($count < $#c);
        }
        
        foreach $id (sort keys %results) {
            
            unless ($id eq 'homologene') {
                print "\nID $id:\n";
                for ($i=0; $i <= $count; $i++) {
                    foreach (sort keys %{$results{$id}}) {
                        print "$_: $results{$id}{$_}[$i]\n" if ($results{$id}{$_}[$i]);
                    }
                    print "\n";
                }
            }
        }
    }
    
    else {
        
        foreach $id (sort keys %results) {
            
            print "\nID $id:\n";
            foreach (sort keys %{$results{$id}}) {
                print "$_: $results{$id}{$_}\n";
            }
        }
    }
}

#***********************************************************

sub print_links {
    
    # Input: %results output from sub elink
    
    my %results = @_;
    my ($key, $db);
    
    foreach $key (sort keys %{$results{from}}) {
        print "Links from: ";
        foreach (@{$results{from}{$key}}) {
            print "$_ ";
        }
        foreach $db (keys %{$results{to}}) {
            print "\nto $db:";
            foreach (@{$results{to}{$db}{$key}}) {
                print "$_ ";
            }
        }
        print "\n***\n";
    }
    
}

#**********************************************************

sub print_link_summaries {
    
    # Input: %results output from sub link_history
    # Output: Docsums for linked records arranged by input UID
    # set and linked database
    
    my %results = @_;
    my (%params,%docsums);
    my ($db, $set);
    
    foreach $set ( sort keys %{$results{to}} ) {
        
        print "Links from set $set\n";
        foreach $db (keys %{$results{to}{$set}} ) {
            
            $params{db} = $db;
            $params{WebEnv} = $results{WebEnv};
            $params{query_key} = $results{to}{$set}{$db}{query_key};
            %docsums = esummary(%params);
            print "$db\n\n";
            print_summary(%docsums);
            print "\n";
        }
    }
    
}

#**********************************************************

sub get_uids {
    
    # Retrieves all UIDs from an Entrez history set
    # Input: %params:
    # $params{WebEnv} - web environment
    # $params{query_key} - query_key
    # $params{db} - database
    # $params{verbose} - prints output message
    # Output: array containing UIDs
    
    my %params = @_;
    my %results;
    my $num;
    
    unless ( ($params{db}) && ($params{WebEnv}) && ($params{query_key}) ) {
        print "\nWARNING: get_uids requires db, query_key and WebEnv keys in input hash!\n";
    }
    $params{verbose} = '' unless ($params{verbose});
    
    $params{usehistory} = 'y';
    $params{term} = "%23$params{query_key}";
    $params{retmax} = 100000000;
    %results = esearch(%params);
    
    $num = @{$results{uids}};
    print "Retrieved $num UIDs from query key $params{query_key}.\n" if ($params{verbose} eq 'y');
    
    return @{$results{uids}};
    
}

#*********************************************************

sub read_index {
    
    # reads index file (.idx) or score file (.sco) produced by elink_by_id or search_links
    # Output: hash %index: $index{id} = comma-delimited list of linked UIDs or scores
    
    my $file = $_[0];
    my %index;
    my ($key, $list);
    
    open (INPUT, "$file") || die "Can't open $file!\n";
    
    while (<INPUT>) {
        
        chomp;
        ($key, $list) = split(/:/, $_);
        $index{$key} = $list;
        
    }
    
    close INPUT;
    
    return %index;
    
}

#*********************************************************

sub get_linknames {
    
    # Using EInfo, collects available link names for given initial and destination databases
    # Input: $dbfrom, $db
    # Output: @linknames - array of link names
    
    my ($dbfrom, $db) = @_;
    my ($url, $out, $rec);
    my @linknames;
    
    $url = $base . "einfo.fcgi?db=$dbfrom";
    
    $out = get($url);
    
    while ($out =~ /<Link>(.*?)<\/Link>/sg) {
        
        $rec = $1;
        $rec =~ /<Name>(.*)<\/Name>/;
        $rec = $1;
        push (@linknames, $rec) if ($rec =~ /$db/);
        
    }
    
    return @linknames;
    
}

#*********************************************************

sub get_link_report {
    
    # Prints a summary report of links found, given an output hash from elink_batch_to or elink_by_id_to
    # Input: %links
    # Output: @linknames - array of link names found
    
    my %links = @_;
    my %sparams;
    my @linknames;
    my ($link, $word);
    
    foreach $link (keys %links) {
        
        $sparams{db} = $links{$link}{db};
        $sparams{term} = "%23$links{$link}{query_key}";
        $sparams{WebEnv} = $links{$link}{WebEnv};
        
        %sparams = esearch(%sparams);
        if ($sparams{count} == 1) { $word = 'link'; }
        else { $word = 'links'; }
        print "$link: $sparams{count} $word\n";
        
        push (@linknames, $link);
    }
    
    return (@linknames);
    
}

#********************************************************

sub extract_links {
    
    # Creates a universal hash (db,query_key,WebEnv) from %links hash (output of elink_batch_to or
    # elink_by_id_to) for linkname
    # Input: $linkname, %links
    # Output: %dblinks
    
    my ($linkname, %links) = @_;
    my %dblinks;
    
    if ($links{$linkname}) {
        
        %dblinks = %{$links{$linkname}};
        
    }
    
    else {
        print "WARNING: No link $linkname in input links hash!\n";
    }
    
    return %dblinks;
    
}

#*********************************************************

sub get_ftp_file {
    
    # Retrieves an ftp url and writes file to $outfile
    # Any file extensions in the original ftp file and added to
    # $outfile to complete the file name
    
    # Input: $url, $outfile;
    # Output: file $outfile (plus file extensions)
    
    my ($url, $outfile) = @_;
    my ($server, $path, $ftp);
    my (@dirs, @parts);
    
    $url =~ /ftp:\/\/(.*?)\/(.*)/;
    
    $server = $1;
    $path = $2;
    
    @dirs = split(/\//, $path);
    
    $file = pop(@dirs);
    
    @parts = split(/\./, $file);
    
    $outfile = 'pubchem' unless ($outfile);
    
    shift(@parts);
    
    $file = join('.', @parts);
    
    $file = $outfile . '.' . $file;
    
    $ftp = Net::FTP->new("$server") || die "Can't connect: $@\n";
    
    $ftp->login("anonymous", "powerscripting") || die "Can't authenticate!\n";
    
    $ftp->binary;
    
    $ftp->get($path, $file) || die "Can't fetch $path: $!\n";
    
    $ftp->quit();
    
    print "Data written to $file\n";
    
}
1;

