#!/usr/bin/perl

# Argument one: plain text, GI file (one per line)
# Argumento two: NCBI api (from your NCBI account)

use strict;
use warnings;
use LWP::Simple;

my $post_url;
my @epost_result_array;
my $line;
my $query_key;
my $web_environment;

    $post_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db=snp&id=242,28853987"; 
    @epost_result_array = get( "$post_url" ); 
        foreach $line ( @epost_result_array ) # Search each returned line for... 
            { if( $line =~ m/<QueryKey>(.*)<\/QueryKey>/ ) # ... the query key, and ... 
                { $query_key = $1; } 
              if( $line =~ m/<WebEnv>(.*)<\/WebEnv>/ ) # ... the Web environment. 
                { $web_environment = $1; } }

# Open GI file
open my $gi_list_fh, '<', "$ARGV[0]"
    or die("Missing GI list");

# Variable declaration
my $counter_input = 0;
my @gi_list;
my @gi_query;
my $input;
my $start1 = 0;
my $start2 = 1;
my $iteration = 1;

# Store GIs in array
while (<$gi_list_fh>) {
    chomp $_;
    my $gi = $_;
    push (@gi_list, $gi);
    }


# Capture 250 gis
while ($iteration < ((scalar @gi_list)/250)) {
    @gi_query = @gi_list[(250*$start1)..($start2*250)];

# URL for retrieving query
$input = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" . join(",", @gi_query) . "&rettype=acc" . "&query_key=$query_key&WebEnv=$web_environment";

# Retrieve 250 accession numbers
print "\nStarting query...\n";
my $queried = get("$input");
print $queried;

$start1++;
$start2++;
$iteration++;

}
