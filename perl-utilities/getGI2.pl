#!/usr/bin/perl

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

# Abre o arquivo de GIs
open my $gi_list_fh, '<', "$ARGV[0]"
    or die("Missing GI list");

# Declaração de variáveis
my $counter_input = 0;
my @gi_list;
my @gi_query;
my $input;
my $start1 = 0;
my $start2 = 1;
my $iteration = 1;

# Estruturar os GIs em array
while (<$gi_list_fh>) {
    chomp $_;
    my $gi = $_;
    push (@gi_list, $gi);
    }


# Loop for a cada 250
while ($iteration < ((scalar @gi_list)/250)) {
    @gi_query = @gi_list[(250*$start1)..($start2*250)];

# URL for retrieving query
$input = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" . join(",", @gi_query) . "&rettype=acc&api_key=6649a6c44ab0d97545450984407d7c2ace08" . "&query_key=$query_key&WebEnv=$web_environment";

#print $input;
# Consulta dos 1000
print "\nStarting query...\n";
my $queried = get("$input");
print $queried;

$start1++;
$start2++;
$iteration++;

}
