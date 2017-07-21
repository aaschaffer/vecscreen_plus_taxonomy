#!/usr/bin/perl -w
# the first line of perl code has to be as above
#
# Author: Alejandro Schaffer
# Documntation help and code review: Eric Nawrocki 
#
#
# Code to evaluate vecscreen matches using taxonomy information about
# the (query) sequence and the (subject) vector interval.
# 
# The sequence is assumed to come from a single taxon and that should
# be specified as either a genus-level taxid or 1, if the genus is
# unknown. In contrast, vectors may be composite sequences that have
# different pieces from different genera, so the matching interval is
# needed.
#
# Vectors are assumed to be entries in UniVec. What is a single vector
# in the wet lab may be split into multiple UniVec entries. In any
# comments below this line, the term "vector" is equivalent to "UniVec
# entry", despite the discrepancy noted above.
#
# Taxonomic source information about various intervals of vectors has
# been pre-computed and is provided as one set of inputs to this
# program.
#
# Some vectors and vector segments are artificial, meaning that they
# were synthesized in a labratory and do not come from a natural
# source.
# 
#
# Usage: compare_vector_matches_wtaxa.pl \
#        --input_summary <input summary file> \
#        --input_taxa <input file with NCBI's taxonomy tree as a list of edges>  \
#        --input_artificial_vectors <vectors that are entirely artificial> \
#        --input_artificial_segments \
#        --input_univec_sources <list of input files for UniVec sources> \
#        --input_amr <known vector intervals that confer AMR> \
#        --input_sequences <FASTA file with sequences being considered> \
#        --outfile <combined output summary file> 

use strict;
use warnings;
use Getopt::Long;

require "epn-options.pm";


# input/output file names
my $input_summary_file;        #input file that summarizes vecscreen matches 
my $input_taxa_file;           #input file with taxonomy in a tab-delimited four-column format 
                               #(taxid, parent taxid, rank, depth where the depth of the root is 1)
my $input_univec_genome_file;  #input file with list of files having matches of UniVec sequences to known genomes
my $input_artificial_vectors;  #input file with entire vectors that are artificial
my $input_artificial_segments; #input file with vector segments that are artificial
my $input_amr_segments;        #input file with vector segments that are known to confer anti-microbial resistance (AMR)
my $input_sequences;           #input file with candidate sequences
my $output_file;               #output summary file with three columns added showing a) evaluation of the match, 
                               #b) taxid of the vector interval if known, and c) lowest common ancestor of the query 
                               #and the subject interval, if applicable or 1 otherwise
my $microsat_file = "Microsatellite_vectors.txt"; #fixed file with names of vectors that contain microsatellites

# variables used in processing input
my $nextline;  #one line of the input summary file
my $m;         #loop index
my @fields;    #entries in one row
my $accession; #accession part of an identifier  

# The following five column indices are for the $input_summary_file with 
# vecscreen matches passed in as the argument --input_summary
my $ACC_COLUMN          = 0; #the column with the identifier
my $GENUS_COLUMN        = 1; #the column with the genus of the sequence represented by the identifier
my $TAXID_COLUMN        = 2; #the column with the genus of the sequence represented by the identifier
my $UNIVEC_COLUMN       = 5; #the column with the identifier of the matching UniVec sequence
my $UNIVEC_START_COLUMN = 6; #the column with the start position in the UniVec sequence
my $UNIVEC_END_COLUMN   = 7; #the column with the end position in the UniVec sequence

# The following four column indices are for the files with 
# UniVec source information passed in as the arguments 
# --input_artificial_segments ($input_artificial_segments) 
#  and --input_univec_sources ($input_artificial_vectors)
my $SOURCE_VECTOR_COLUMN = 0;     #an entry in UniVec
my $SOURCE_START_COLUMN  = 1;     #start of a vector interval
my $SOURCE_END_COLUMN    = 2;     #end of a vector interval
my $SOURCE_TAXON_COLUMN  = 3;     #the taxid (which is the source) for this vector interval 

# special taxon constants
my $ARTIFICIAL_TAXON     = 81077; #taxon to use for artificial sequences
my $UNCULTURED_BACTERIA_TAXON = 77133; #taxid of uncultured bacteria
my $TAXONOMY_ROOT        = 1;     #taxon to use as a default when no other taxon fits the situation
my $BACTERIA_TAXON       = 2;     #taxon that is the root of all Bacteria

# variables used in constructing output
my $output_result;       # decision about a single vector match
my $output_vector_taxon; # taxonomic information about the vector segment to support the decision
my $output_lca;          # lowest common ancestor, if any, between the sequence and the vector

# data structures for keeping track of information on matches
my %genome_match_taxon;        # hash of arrays, what is the taxon of the ith genome match, first layer by taxonomic rank
my %genome_match_univec;       # hash of arrays, what is the UniVec identifier of the ith genome match, first layer by taxonomic rank
my %genome_match_start;        # hash of arrays, what is the start position of the ith genome match, first layer by taxonomic rank
my %genome_match_end;          # hash of arrays, what is the end position of the ith genome match, first layer by taxonomic rank
my @artificial_segment_univec; # array, what is the vector of the ith artificial segment match
my @artificial_segment_start;  # array, what is the start position of the ith artificial segment match
my @artificial_segment_end;    # array, what is the end position of the ith artificial segment match
my @amr_segment_univec;        # array, what is the vector of the ith AMR segment match
my @amr_segment_start;         # array, what is the start position of the ith AMR segment match
my @amr_segment_end;           # array, what is the end position of the ith AMR segment match

# hashes for match indices
my %artificial_first_match_index; # hash, for each vector with an artificial segment, what is the first index of its matches 
                                  # in the artificial_segment_start and artificial_segment_end arrays
my %artificial_last_match_index;  # hash, for each vector with an artificial segment, what is the last index of its matches 
                                  # in the artificial_segment_start and artificial_segment_end arrays
my %biological_first_match_index; # for each taxonomic rank and for each vector with a biological segment, what is the first 
                                  # index of its matches in the genome_match_start and genome_match_end hashes of arrays
my %biological_last_match_index;  # for each taxonomic rank and for each vector with an artificial segment, what is the last 
                                  # index of its matches in the genome_match_start and genome_match_end hashes of arrays
my %amr_first_match_index;        # for each vector with an AMR segment, what is the first index of its matches in the 
                                  # amr_segment_start and amr_segment_end arrays
my %amr_last_match_index;         # for each vector with an AMR segment, what is the last index of its matches in the 
                                  # amr_segment_start and amr_segment_end arrays

# hashes mapping taxon to other information
my %taxonomy_parent;              # hash, maps a taxon to its parent
my %taxonomy_level;               # hash, maps a taxon to its level (a.k.a. depth) in the tree, where the level of the root is 1
my %rank_hash;                    # hash, maps each taxon to its taxonomic rank (e.g., phylum)
my %belongs_in_bacteria;          # hash, maps each taxon to 0 or 1 depending on whether it belongs in bacteria (1)

# variables used for processing the main input file, the vecscreen summary file
my $i;                      # loop index
my $num_columns_input = 11;  # number of columns in input vecscreen summary
my $start_index;            # first index to check for matches (loop lower bound)
my $end_index;              # last index to check for matches (loop upper bound)
my $this_vector;            # the vector for one match to UniVec
my $one_start;              # start of match in a UniVec sequence
my $one_end;                # end of match in a UniVec sequence
my $one_genus;              # genus (in the taxonomy tree) of one query sequence in vecscreen summary
my $one_taxid;              # taxid (in the taxonomy tree) of one query sequence in vecscreen summary

# hard-coded overlap thresholds
my $overlap_threshold     = 16; #amount of overlap between UniVec match and genome match to be considered a biological overlap
my $amr_overlap_threshold = 12; #amount of overlap between UniVec match and AMR segment to be considered FALSE_AMR

# variables used when parsing artificial vectors files  and candidates file
my %artificial_whole_vectors; # hash to store the identifiers of whole vectors that are artificial, for quick look up
my %microsat_containing_vectors; # hash to store the identifiers of whole vectors that contain a microsatellite, for quick look up
my %deflines; #maps an accession to its defline
my @taxonomy_ranks = ("superkingdom","kingdom","phylum","class","order","family","tribe","genus"); #these are the ranks for which matches to UniVec entries have been collected

my $num_taxonomy_ranks = scalar(@taxonomy_ranks); # size of @taxonomy_ranks
my $rank_index;        #loop index for iteration in @taxonomy_ranks
my $one_rank;          #entry in @taxonomy_ranks
my $this_lca;          #one return value from the subroutine find_lca
my $this_vector_taxid; #the taxid for which this vector matches as a source

# variables used in special debugging mode
my $DEBUG             = 0;  # set to 1 to enter debug mode
my $genus_to_test     = 10373;     
my $accession_to_test = "M93360.1";
my $vector_to_test    = "uv|DQ649431.1:3092-4131-49";

# variables for dealing with options (epn-options.pm)
my %opt_HH           = ();
my @opt_order_A      = ();
my %opt_group_desc_H = ();

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
# The opt_Add() function is the way we add options to %opt_HH.
# It takes values of for each of the 2nd dim keys listed above.
$opt_group_desc_H{"1"} = "basic options";
#       option                         type      default group   requires incompat preamble-outfile                                                                help-outfile
opt_Add("-h",                          "boolean",0,          0,    undef, undef,   undef,                                                                          "display this help",                                       \%opt_HH, \@opt_order_A);
opt_Add("--input_summary",             "string", undef,      1,    undef, undef,   "input file of vecscreen summary",                                              "File name <s> of vecscreen summary",                      \%opt_HH, \@opt_order_A);
opt_Add("--input_taxa",                "string", undef,      1,    undef, undef,   "input file mapping vescreen matches to taxa",                                  "File name <s> mapping vecscreen matches to taxa",         \%opt_HH, \@opt_order_A);
opt_Add("--input_artificial_vectors",  "string", undef,      1,    undef, undef,   "input file listing vectors that are entirely artificial",                      "File name <s> with vectors that are entirely artificial", \%opt_HH, \@opt_order_A);
opt_Add("--input_artificial_segments", "string", undef,      1,    undef, undef,   "input file listing vector segments that are artificial but not whole vectors", "File name <s> with vector segments that are artificial but not whole vectors", \%opt_HH, \@opt_order_A);
opt_Add("--input_univec_sources",      "string", undef,      1,    undef, undef,   "input file listing multiple files with univec sources",                        "File name <s> listing multiple files with univec sources", \%opt_HH, \@opt_order_A);
opt_Add("--input_amr",                 "string", undef,      1,    undef, undef,   "input file listing vector segments that confer anti-microbial resistance",     "File name <s> with vector segments that confer AMR",      \%opt_HH, \@opt_order_A);
opt_Add("--input_sequences",           "string", undef,      1,    undef, undef,   "input file listing candidate sequences in FASTA format",                       "File name <s> with candidate sequences",                  \%opt_HH, \@opt_order_A);
opt_Add("--outfile",                   "string", undef,      1,    undef, undef,   "output file to create",                                                        "Name <s> of output file to create",                       \%opt_HH, \@opt_order_A);


my $synopsis = "compare_vector_matches_wtaxa.pl: evaluate vecscreen matches taking taxonomy of query and match into account\n";
my $usage    = "Usage:\n\n"; 
$usage      .= "\tcompare_vector_matches_wtaxa.pl \\\n";
$usage      .= "\t--input_summary <input vecscreen summary file> \\\n"; 
$usage      .= "\t--input_taxa <input file with NCBI's taxonomy> \\\n";
$usage      .= "\t--input_artificial_vectors <vectors that are entirely artificial> \\\n";
$usage      .= "\t--input_artificial_segments <vector segments that are artificial> \\\n";
$usage      .= "\t--input_univec_sources <list of input files for UniVec sources> \\\n";
$usage      .= "\t--input_amr <known vector intervals that confer AMR> \\\n";
$usage      .= "\t--outfile <output file>\n\n";
$usage      .= "\nFor example:\n";
$usage      .= "compare_vector_matches_wtaxa.pl --input_summary summary_combined.txt "; 
$usage      .= "--input_taxa taxonomy.txt ";
$usage      .= "--input_artificial_vectors artificial_whole_vectors.txt ";
$usage      .= "--input_artificial_segments artificial_segments.txt ";
$usage      .= "--input_univec_matches univec_sources.txt ";
$usage      .= "--input_amr amr_vector_segments.txt ";
$usage      .= "--input_sequences candidate_sequences.txt ";
$usage      .= "--outfile updated_summary.txt\n";

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $options_okay =
    &GetOptions('h'                            => \$GetOptions_H{"-h"},
                'input_summary=s'              => \$GetOptions_H{"--input_summary"},
                'input_taxa=s'                 => \$GetOptions_H{"--input_taxa"},
                'input_artificial_vectors=s'   => \$GetOptions_H{"--input_artificial_vectors"},
                'input_artificial_segments=s'  => \$GetOptions_H{"--input_artificial_segments"},
                'input_univec_sources=s'       => \$GetOptions_H{"--input_univec_sources"},
                'input_amr=s'                  => \$GetOptions_H{"--input_amr"},
                'input_sequences=s'            => \$GetOptions_H{"--input_sequences"},
                'outfile=s'                    => \$GetOptions_H{"--outfile"});

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) {
    opt_OutputHelp(*STDOUT, $synopsis . $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
    if(! $options_okay) { die "ERROR, unrecognized option;"; }
    else                { exit 0; } # -h, exit with 0 status
}

# set options in %opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

# define file names
$input_summary_file        = opt_Get("--input_summary", \%opt_HH);
$input_taxa_file           = opt_Get("--input_taxa", \%opt_HH);
$input_artificial_vectors  = opt_Get("--input_artificial_vectors", \%opt_HH);
$input_artificial_segments = opt_Get("--input_artificial_segments", \%opt_HH);
$input_univec_genome_file  = opt_Get("--input_univec_sources", \%opt_HH);
$input_amr_segments        = opt_Get("--input_amr", \%opt_HH);
$input_sequences           = opt_Get("--input_sequences", \%opt_HH);
$output_file               = opt_Get("--outfile", \%opt_HH);

# die if any of the required options were not used
my $errmsg = undef;
if(! defined $input_summary_file)        { $errmsg .= "ERROR, --input_summary option not used.\n"; }
if(! defined $input_taxa_file)           { $errmsg .= "ERROR, --input_taxa option not used.\n"; }
if(! defined $input_artificial_vectors)  { $errmsg .= "ERROR, --input_artificial_vectors option not used.\n"; }
if(! defined $input_artificial_segments) { $errmsg .= "ERROR, --input_artificial_segments option not used.\n"; }
if(! defined $input_univec_genome_file)  { $errmsg .= "ERROR, --input_univec_sources option not used.\n"; }
if(! defined $input_amr_segments)        { $errmsg .= "ERROR, --input_amr option not used.\n"; }
if(! defined $input_sequences)           { $errmsg .= "ERROR, --input_sequences option not used.\n"; }
if(! defined $output_file)               { $errmsg .= "ERROR, --outfile option not used.\n"; }
if(defined $errmsg) { 
  die $errmsg . "\n$usage\n";
}


# open output files
open(SUMMARY, "<", $input_summary_file) or die "Cannot open 1 $input_summary_file for input\n"; 
open(OUTPUT,  ">", $output_file)        or die "Cannot open 2 $output_file for output\n"; 

# process input files and store the relevant information 
process_taxonomy_tree($input_taxa_file);
process_microsat_vectors($microsat_file);
process_candidate_sequences($input_sequences);
process_artificial_vectors($input_artificial_vectors);
process_artificial_segments($input_artificial_segments);
process_amr_segments($input_amr_segments);
process_genome_match_info($input_univec_genome_file);

######################################################################
#
# Process the summary input file, one line at a time. Each line
# includes a query (the sequence in which the vector match was found
# by vecscreen) and a subject (the vector match). For each line 
# (query/subject pair) we will classify it into one of 5 possible 
# evaluations:
#
# 1. NO_DATA:           there is no data about the source of the vector segment 
#                       in the match.
# 
# 2. TRUE_ARTIFICIAL:   the vector segment matched is an ARTIFICIAL sequence and 
#                       hence the match is TRUE contamination.
# 3. TRUE_ARTIFICIAL_MICROSAT:   the vector segment matched is an ARTIFICIAL sequence and 
#                       hence the match is TRUE contamination and the vector contains a microsatellite.
#
# 4. FALSE_AMR:         The query sequence is bacterial; the vector segment matches 
#                       a known sequence that confers anti-microbial resistance 
#                       and these can often be transferred horizontally between
#                       bacteria that may be taxonomically distant.
#
# 5. FALSE_BIOLOGICAL:  The subject's biological origin is known and its taxid is 
#                       deemed close enough to
#                       that of the query, so that the match is not contamination.
#
# 6. TRUE_BIOLOGICAL :  The query and the matching subject originate from
#                       taxa that are too far apart for the vector to occur
#                       plausibly in the query.
# 7. TRUE_MICROSAT :  The query and the matching subject originate from
#                       taxa that are too far apart for the vector to occur
#                       plausibly in the query, and the vector contains a microsatellite.
#
# 8. LIKE_FALSE_BACTERIAL: The query is from uncultered bacteria and the matching subject is
#                        from bacteria
#
# For each possibility, we further determine taxonomic information about the 
# query/subject pair, namely the $output_vector_taxon and $output_lca, as follows:
#
# $output_result     $output_vector_taxon  $output_lca   notes
# ----------------   --------------------  -----------   -------------------------
# TRUE_ARTIFICIAL    81077                 81077         81077 is hard-coded taxonomy value for artificical sequences ($ARTIFICIAL_TAXON)
# TRUE_ARTIFICIAL_MICROSAT    81077                 81077         81077 is hard-coded taxonomy value for artificical sequences ($ARTIFICIAL_TAXON)
# FALSE_AMR          2                     2             2     is hard-coded taxonomy value for the root of all bacteria ($BACTERIA_TAXON)
# FALSE_BIOLOGICAL   $this_vector_taxid    $this_lca     determined for each vector match using NCBI taxonomy 
# TRUE_BIOLOGICAL    $this_vector_taxid    $this_lca     determined for each vector match using NCBI taxonomy 
# TRUE_MICROSAT    $this_vector_taxid    $this_lca     determined for each vector match using NCBI taxonomy 
#
while(defined($nextline = <SUMMARY>)) {
    chomp($nextline);
    # default initialization of three output fields
    $output_result       = "NO_DATA";
    $output_vector_taxon = $TAXONOMY_ROOT;
    $output_lca          = $TAXONOMY_ROOT;
    @fields = split /\t/, $nextline;
    $accession   = $fields[$ACC_COLUMN];
    $one_genus   = $fields[$GENUS_COLUMN];
    $one_taxid   = $fields[$TAXID_COLUMN];
    $this_vector = $fields[$UNIVEC_COLUMN];
    if ($fields[$UNIVEC_START_COLUMN] < $fields[$UNIVEC_END_COLUMN]) {
	$one_start = $fields[$UNIVEC_START_COLUMN];
	$one_end   = $fields[$UNIVEC_END_COLUMN];
    }
    else {
	$one_start = $fields[$UNIVEC_END_COLUMN];
	$one_end   = $fields[$UNIVEC_START_COLUMN];
    }
    if (defined($artificial_whole_vectors{$this_vector})) {
        # TRUE_ARTIFICIAL:   the vector segment matched is an ARTIFICIAL sequence 
        #                    (whole artificial sequence) and hence the match is TRUE 
        #                    contamination.
	if (defined($microsat_containing_vectors{$this_vector}) && defined($deflines{$accession}) && ($deflines{$accession} =~m/satellite/ )) {
	    $output_result       = "TRUE_ARTIFICIAL_MICROSAT";
	}
	else {
	    $output_result       = "TRUE_ARTIFICIAL";
	}
	$output_vector_taxon = $ARTIFICIAL_TAXON;
	$output_lca          = $ARTIFICIAL_TAXON;
    }
    else { 
	if (defined($artificial_first_match_index{$this_vector})) {
	    $start_index = $artificial_first_match_index{$this_vector};
	    $end_index   = $artificial_last_match_index{$this_vector};
	    for($m = $start_index; $m <= $end_index; $m++) {
		#things are set up so that in each pair of positions start <= end 
		if (find_overlap($one_start, $one_end, $artificial_segment_start[$m], $artificial_segment_end[$m]) >= $overlap_threshold) {
                    # TRUE_ARTIFICIAL:   the vector segment matched is an ARTIFICIAL sequence 
                    #                    (segment of an artificial sequence) and hence the match
                    #                    is TRUE contamination.
		    if (defined($microsat_containing_vectors{$this_vector}) && defined($deflines{$accession}) && ($deflines{$accession} =~m/satellite/ )) {
			$output_result       = "TRUE_ARTIFICIAL_MICROSAT";
		    }
		    else {
			$output_result       = "TRUE_ARTIFICIAL";
		    }
		    $output_vector_taxon = $ARTIFICIAL_TAXON;
		    $output_lca          = $ARTIFICIAL_TAXON;
		    last;
		}
	    }
	}
	if ("NO_DATA" eq $output_result) {
	    if ((defined($amr_first_match_index{$this_vector})) && 
		(($UNCULTURED_BACTERIA_TAXON == $one_taxid) || ($BACTERIA_TAXON == find_ancestor($one_genus,"superkingdom")))) {
		$start_index = $amr_first_match_index{$this_vector};
		$end_index   = $amr_last_match_index{$this_vector};
		for($m = $start_index; $m <= $end_index; $m++) {
		    #things are set up so that in each pair of positions start <= end 
		    if (find_overlap($one_start, $one_end, $amr_segment_start[$m], $amr_segment_end[$m]) >= $amr_overlap_threshold) {
                        # FALSE_AMR:  The query sequence is bacterial; the vector segment matches 
                        #             a known sequence that confers anti-microbial resistance 
                        #             and these can often be transferred horizontally between
                        #             bacteria that may be taxonomically distant.
			$output_result       = "FALSE_AMR";
			$output_vector_taxon = $BACTERIA_TAXON;
			$output_lca          = $BACTERIA_TAXON;
			last;
		    }
		}
	    }
	}
    }

    # If we didn't define $this_vector as TRUE_ARTIFICIAL or FALSE_AMR  
    # we now try to match from most specific to least specific biological taxa.
    # We go in this particular order because a match can be changed from TRUE_BIOLOGICAL 
    # at a lower rank to FALSE_BIOLOGICAL at a higher rank, but not the other way around
    for($rank_index = $num_taxonomy_ranks -1; $rank_index >= 0 ; $rank_index--) {
        $one_rank = $taxonomy_ranks[$rank_index];
	if (("NO_DATA" eq $output_result) || ("TRUE_BIOLOGICAL" eq $output_result)) {
	    if (defined($biological_first_match_index{$one_rank}{$this_vector})) {
		$start_index = $biological_first_match_index{$one_rank}{$this_vector};
		$end_index   = $biological_last_match_index{$one_rank}{$this_vector};
		for($m = $start_index; $m <= $end_index; $m++) {
		    if (("NO_DATA" eq $output_result) || ("TRUE_BIOLOGICAL" eq $output_result)) {
			if ($DEBUG) {
			    if ($this_vector eq $vector_to_test) {
				print "Testing rank $one_rank for vector $this_vector for taxid $genome_match_taxon{$one_rank}[$m] with ends $one_start $one_end $genome_match_start{$one_rank}[$m] $genome_match_end{$one_rank}[$m]\n";
			    }
			}
			if (find_overlap($one_start, $one_end, $genome_match_start{$one_rank}[$m], $genome_match_end{$one_rank}[$m]) >= $overlap_threshold) {
			    $this_vector_taxid = $genome_match_taxon{$one_rank}[$m];
			    if ((1 == $one_genus) && ($UNCULTURED_BACTERIA_TAXON != $one_taxid)) {
				$this_lca          = find_lca($one_taxid, $this_vector_taxid);
			    }
			    else {
				$this_lca          = find_lca($one_genus, $this_vector_taxid);
			    }
			    if ($DEBUG) {
				if ($this_vector eq $vector_to_test) {
				    print "Testing rank $one_rank lca is $this_lca and taxid is $this_vector_taxid\n";
				}
			    }
			    if ($this_lca == $this_vector_taxid) {
				#  FALSE_BIOLOGICAL:  The subject's biological origin is known and its taxid is 
				#                     deemed close enough to that of the query, so that the 
				#                     match is not contamination.
				$output_result       = "FALSE_BIOLOGICAL";
				$output_vector_taxon = $this_vector_taxid;
				$output_lca          = $this_lca;
			    }
			    else {
				if (($UNCULTURED_BACTERIA_TAXON == $fields[$TAXID_COLUMN]) && ($belongs_in_bacteria{$this_vector_taxid})) {
				    $output_result = "LIKELY_FALSE_BACTERIAL";
				    $output_vector_taxon = $this_vector_taxid;
				    $output_lca = $UNCULTURED_BACTERIA_TAXON;
				}
				else {
                                # TRUE_BIOLOGICAL:   The query and the matching subject originate from
                                #                    taxa that are too far apart for the vector to occur
                                #                    plausibly in the query (at least so far, this 
                                #                    may be updated to FALSE_BIOLOGICAL as $m gets larger
                                #                    (taxonomic rank gets higher).
				    $output_result = "TRUE_BIOLOGICAL";
				    $output_vector_taxon = $this_vector_taxid;
				    $output_lca = $this_lca;     
				}
			    }
			}
		    }
		}
	    }
	}
    }
    # output original information from input summary file ($input_summary_file)
    for($i = 0; $i < $num_columns_input; $i++) {
	print OUTPUT "$fields[$i]\t";
    }
    # output additional 3 columns we just determined
    if (("TRUE_BIOLOGICAL" eq $output_result) && (defined($microsat_containing_vectors{$this_vector})) && defined($deflines{$accession}) && ($deflines{$accession} =~m/satellite/ )) {
	$output_result = "TRUE_MICROSAT";
    }
    print OUTPUT "$output_result\t$output_vector_taxon\t$output_lca";  
    print OUTPUT "\n";
}
close (SUMMARY);
close (OUTPUT);

################################################
# SUBROUTINES
################################################
# 
#
# List of subroutines:
# 
# Subroutines that parse input files: 
# process_taxonomy_tree():       parse taxonomy tree file ($input_taxa_file)
# process_microsat_vectors();    record which vectors contain microsatellites
# process_artificial_vectors():  parse artificial vectors file ($input_artificial_vectors)
# process_artificial_segments(): parse artificial segments file ($input_artificial_segments)
# process_amr_segments():        parse AMR segments file ($input_amr_segments)
# process_genome_match_info():   parse genome match info file ($input_univec_genome_file)
# recognize_rank():              helper function for process_genome_match_info() 
#
# Other subroutines:
# find_overlap():  return number of overlapping positions between two intervals
# find_ancestor(): return ancestor at a specified rank for a taxon
# find_lca():      return lowest common ancestor of two integer taxids
#
################################################
# Subroutine: process_taxonomy_tree()
# Synopsis: reads a file that includes NCBI's taxonomy information in four columns (taxon, parent taxon, rank, level)
#
# Args: $taxonomy_information_file
#
# Returns: nothing

sub process_taxonomy_tree {

    my $sub_name = "process_taxonomy_tree()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_taxonomy_file) = @_;
    my $local_nextline; #one line of taxonomy information
    my @local_fields;   #split of one line of taxonomy information
    my $local_TAXID_COLUMN       = 0;
    my $local_PARENT_COLUMN      = 1;
    my $local_FORMAL_RANK_COLUMN = 2;
    my $local_LEVEL_COLUMN       = 3;
    my $local_BACTERIA_COLUMN       = 4;

    open(TAXONOMY, "<", $local_taxonomy_file) or die "Cannot open $local_taxonomy_file for input in $sub_name\n"; 
    
    while(defined($local_nextline = <TAXONOMY>)) {
	chomp($local_nextline);
	@local_fields = split /\t/, $local_nextline;
	$taxonomy_parent{$local_fields[$local_TAXID_COLUMN]} = $local_fields[$local_PARENT_COLUMN];
	$taxonomy_level{$local_fields[$local_TAXID_COLUMN]}  = $local_fields[$local_LEVEL_COLUMN];
	$belongs_in_bacteria{$local_fields[$local_TAXID_COLUMN]}  = $local_fields[$local_BACTERIA_COLUMN];
	if ($DEBUG && ($genus_to_test == $local_fields[0])) {
	    print "Testing taxon $genus_to_test\n"; 
	}
        $rank_hash{$local_fields[$local_TAXID_COLUMN]} = $local_fields[$local_FORMAL_RANK_COLUMN];
    }
    close(TAXONOMY);
}

# Subroutine: process_candidate_sequences()
# Synopsis: reads a file that contains a list of input sequences in FASTA and stores the deflines in a hash
#
# Args: $candidate_sequences_file
#
# Returns: nothing

sub process_candidate_sequences {

    my $sub_name = "process_candidate_sequences()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_sequences_file) = @_;
    my $local_nextline; #one line with one vector
    my $local_accession;

    open(SEQUENCES, "<", $local_sequences_file) or die "Cannot open $local_sequences_file for input in $sub_name\n"; 
    
    while(defined($local_nextline = <SEQUENCES>)) {
	chomp($local_nextline);
	if ($local_nextline =~m/^>/) {
	    ($local_accession) = ($local_nextline =~m/^>(\S+)/);
	    $deflines{$local_accession} = $local_nextline;
	}
    }
    close(SEQUENCES);
}

# Subroutine: process_microsat_vectors()
# Synopsis: reads a file that contains a list of vectors that contain a microsatellite and stores them in a hash
#
# Args: $microsat_vector_file
#
# Returns: nothing

sub process_microsat_vectors {

    my $sub_name = "process_microsat_vectors()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_vector_file) = @_;
    my $local_nextline; #one line with one vector

    open(VECTORS, "<", $local_vector_file) or die "Cannot open $local_vector_file for input in $sub_name\n"; 
    
    while(defined($local_nextline = <VECTORS>)) {
	chomp($local_nextline);
	$microsat_containing_vectors{$local_nextline} = 1;
    }
    close(VECTORS);
}

# Subroutine: process_artificial_vectors()
# Synopsis: reads a file that contains a list of artificial vectors and stores them in a hash
#
# Args: $artificial_vector_file
#
# Returns: nothing

sub process_artificial_vectors {

    my $sub_name = "process_artificial_vectors()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_vector_file) = @_;
    my $local_nextline; #one line with one vector

    open(VECTORS, "<", $local_vector_file) or die "Cannot open $local_vector_file for input in $sub_name\n"; 
    
    while(defined($local_nextline = <VECTORS>)) {
	chomp($local_nextline);
	$artificial_whole_vectors{$local_nextline} = 1;
    }
    close(VECTORS);
}

# Subroutine: process_artificial_segments()
# Synopsis: reads a file that contains UniVec segments (but not whole vectors) that are artificial
#           The file should be sorted by the first column, if not we will detect it and die.
#
# Args: $artificial_segment_file
#       
# Returns: nothing

sub process_artificial_segments {

    my $sub_name = "process_artificial_segments()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_artificial_segment_file) = @_;
    my $local_nextline; #one line of genome match
    my @local_fields; #split of one line of genome match information
    my $local_line_index = 0; #index for entries in file
    my $prv_source_vector = undef; # source vector ($local_fields[$SOURCE_VECTOR_COLUMN]) of previous line
                                   # used to make sure the file is sorted appropriately (by first column)

    open(SEGMENTS, "<", $local_artificial_segment_file) or die "Cannot open $local_artificial_segment_file for input in $sub_name\n"; 
    
    while(defined($local_nextline = <SEGMENTS>)) {
	chomp($local_nextline);
	@local_fields = split /\t/, $local_nextline;
        # We take advantage of fact that file is sorted by first column. If it is not,
        # we will detect it and die.
	if ($ARTIFICIAL_TAXON == $local_fields[$SOURCE_TAXON_COLUMN]) {
            # EPN  added code here to check that file is sorted appropriately and die if not
            #              added code here to check if start <= end and die if not
            # is this a new source vector? 
            if((defined $prv_source_vector) && ($prv_source_vector ne $local_fields[$SOURCE_VECTOR_COLUMN])) { 
                # yes it is new, $artificial_first_match_index{$local_fields[$SOURCE_VECTOR_COLUMN]} should be undef
                # or else the file wasn't sorted by the first column as we expected it to be
                if(defined($artificial_first_match_index{$local_fields[$SOURCE_VECTOR_COLUMN]})) { 
                    die "ERROR in $sub_name, $local_artificial_segment_file not properly sorted by first column, we saw $local_fields[$SOURCE_VECTOR_COLUMN] as the first token in two non-contiguous groups of >= 1 line"; 
                }
            }
            if($local_fields[$SOURCE_START_COLUMN] > $local_fields[$SOURCE_END_COLUMN]) { 
              die "ERROR in $sub_name, $local_artificial_segment_file includes line with start > end:\n$local_nextline\n";
            }
	    $artificial_segment_univec[$local_line_index] = $local_fields[$SOURCE_VECTOR_COLUMN];
	    $artificial_segment_start[$local_line_index]  = $local_fields[$SOURCE_START_COLUMN];
	    $artificial_segment_end[$local_line_index]    = $local_fields[$SOURCE_END_COLUMN];
	    if (!defined($artificial_first_match_index{$local_fields[$SOURCE_VECTOR_COLUMN]})) {
		$artificial_first_match_index{$local_fields[$SOURCE_VECTOR_COLUMN]} = $local_line_index;
	    }
	    $artificial_last_match_index{$local_fields[$SOURCE_VECTOR_COLUMN]} = $local_line_index;
            $prv_source_vector = ($local_fields[$SOURCE_VECTOR_COLUMN]);
	}
	$local_line_index++;
    }
    close(SEGMENTS);
}

# Subroutine: process_amr_segments()
# Synopsis: reads a file that contains UniVec segments (but not whole vectors) that confer anti-microbial resistance.
#           The file should be sorted by the first column, if not we will detect it and die.
#
# Args: $amr_segment_file
#       
# Returns: nothing

sub process_amr_segments {

    my $sub_name = "process_amr_segments()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_amr_segment_file) = @_;
    my $local_nextline; #one line of genome match
    my @local_fields; #split of one line of genome match information
    my $local_line_index = 0; #index for entries in file
    my $prv_source_vector = undef; # source vector ($local_fields[$SOURCE_VECTOR_COLUMN]) of previous line
                                   # used to make sure the file is sorted appropriately (by first column)

    open(SEGMENTS, "<", $local_amr_segment_file) or die "Cannot open $local_amr_segment_file for input in $sub_name\n"; 
    
    while(defined($local_nextline = <SEGMENTS>)) {
	chomp($local_nextline);
	@local_fields = split /\t/, $local_nextline;
        # We take advantage of fact that file is sorted by first column. If it is not,
        # we will detect it and die.
        # EPN  added code here to check that file is sorted appropriately and die if not
        #              added code here to check if start <= end and die if not
        # is this a new source vector? 
        if((defined $prv_source_vector) && ($prv_source_vector ne $local_fields[$SOURCE_VECTOR_COLUMN])) { 
            # yes it is new, $amr_first_match_index{$local_fields[$SOURCE_VECTOR_COLUMN]} should be undef
            # or else the file wasn't sorted by the first column as we expected it to be
            if(defined($amr_first_match_index{$local_fields[$SOURCE_VECTOR_COLUMN]})) { 
                      die "ERROR in $sub_name, $local_amr_segment_file not properly sorted by first column, we saw $local_fields[$SOURCE_VECTOR_COLUMN] as the first token in two non-contiguous groups of >= 1 line"; 
            }
        }
        if($local_fields[$SOURCE_START_COLUMN] > $local_fields[$SOURCE_END_COLUMN]) { 
          die "ERROR in $sub_name, $local_amr_segment_file includes line with start > end:\n$local_nextline\n";
        }
	$amr_segment_univec[$local_line_index] = $local_fields[$SOURCE_VECTOR_COLUMN];
	$amr_segment_start[$local_line_index]   = $local_fields[$SOURCE_START_COLUMN];
	$amr_segment_end[$local_line_index]     = $local_fields[$SOURCE_END_COLUMN];
	if (!defined($amr_first_match_index{$local_fields[$SOURCE_VECTOR_COLUMN]})) {
	    $amr_first_match_index{$local_fields[$SOURCE_VECTOR_COLUMN]} = $local_line_index;
	}
	$amr_last_match_index{$local_fields[$SOURCE_VECTOR_COLUMN]} = $local_line_index;
        $prv_source_vector = ($local_fields[$SOURCE_VECTOR_COLUMN]);
	$local_line_index++;
    }
    close(SEGMENTS);
}

# Subroutine: process_genome_match_info()
# Synopsis: reads a file that contains a list of files (one per taxonomic rank) 
#           such that each file contains information on matches between UniVec sequences and genomes.
#           In the top level file, there are two values on each row, the file name and the taxonomic rank.
#           There is expected to be at most one file per rank.
#           
# Args: $genome_match_file
#       
# Returns: nothing

sub process_genome_match_info {

    my $sub_name = "process_genome_match_info()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_genome_match_file) = @_;
    my $local_file_line; #one line with a file name and a rank
    my %local_ranks_seen; #taxonomic ranks seen
    my $local_one_rank_file; #file of matches for one taxonomic rank
    my $local_rank; #the rank for this file;
    my $local_nextline; #one line of genome match
    my @local_file_fields; #split of line in the upper level file
    my @local_fields; #split of one line of genome match information
    my $local_line_index = 0; #index for entries in file
    my $prv_source_vector = undef; # source vector ($local_fields[$SOURCE_VECTOR_COLUMN]) of previous line
                                   # used to make sure the file is sorted appropriately (by first column)
    my %read_rank = ();    # hash, key: rank $r, value is set to 1 once we open a file for rank $r
                           # used only to make sure we don't read two files for a single rank
                           # if we do, we die in error

    open(FILES_BY_RANK, "<", $local_genome_match_file) or die "Cannot open $local_genome_match_file for input in $sub_name\n"; 
    while(defined($local_file_line = <FILES_BY_RANK>)) {
	chomp($local_file_line);
	if ($DEBUG) {
	    printf "Reading file with rank: $local_file_line\n";
	}
	@local_file_fields = split /\t/, $local_file_line;
	$local_one_rank_file = $local_file_fields[0];
	if ($DEBUG) {
	    printf "Reading file: $local_one_rank_file\n";
	}
        $local_rank = $local_file_fields[1];
	if ($DEBUG) {
	    printf "Working on rank $local_rank\n";
	}
	if (recognize_rank($local_rank)) {
	    open(GENOME, "<", $local_one_rank_file) or die "Cannot open $local_one_rank_file in $sub_name\n"; 
            # EPN added this block to check that multiple files are not read for same rank, we die if
            #              we do
            if((defined $read_rank{$local_rank}) && ($read_rank{$local_rank} == 1)) { 
              die "ERROR in $sub_name, tried to read multiple files for rank $local_rank in file $local_genome_match_file, should only be one file per rank";
            }
            $read_rank{$local_rank} = 1;
            $prv_source_vector = undef;
	    $local_line_index = 0;
	    while(defined($local_nextline = <GENOME>)) {
		chomp($local_nextline);
		@local_fields = split /\t/, $local_nextline;
                # We take advantage of fact that file is sorted by first column. If it is not,
                # we will detect it and die.
                # EPN added code here to check that file is sorted appropriately and die if not
                #              added code here to check if start <= end and die if not
                # is this a new source vector? 
                if((defined $prv_source_vector) && ($prv_source_vector ne $local_fields[$SOURCE_VECTOR_COLUMN])) { 
                    # yes it is new, $amr_first_match_index{$local_fields[$SOURCE_VECTOR_COLUMN]} should be undef
                    # or else the file wasn't sorted by the first column as we expected it to be
                    if(defined($biological_first_match_index{$local_rank}{$local_fields[$SOURCE_VECTOR_COLUMN]})) { 
                        die "ERROR in $sub_name, $local_one_rank_file read from $local_genome_match_file not properly sorted by first column, we saw $local_fields[$SOURCE_VECTOR_COLUMN] as the first token in two non-contiguous groups of >= 1 line"; 
                    }
                }
                if($local_fields[$SOURCE_START_COLUMN] > $local_fields[$SOURCE_END_COLUMN]) { 
                  die "ERROR in $sub_name, $local_one_rank_file read from $local_genome_match_file includes line with start > end:\n$local_nextline\n";
                }
		$genome_match_univec{$local_rank}[$local_line_index] = $local_fields[$SOURCE_VECTOR_COLUMN];
		$genome_match_start{$local_rank}[$local_line_index]  = $local_fields[$SOURCE_START_COLUMN];
		$genome_match_end{$local_rank}[$local_line_index]    = $local_fields[$SOURCE_END_COLUMN];
		$genome_match_taxon{$local_rank}[$local_line_index]  = $local_fields[$SOURCE_TAXON_COLUMN];
		if ($DEBUG) {
		    if ($vector_to_test eq $local_fields[$SOURCE_VECTOR_COLUMN]) {
			print "$local_nextline\n";
		    }
		}
		if (!defined($biological_first_match_index{$local_rank}{$local_fields[$SOURCE_VECTOR_COLUMN]})) {
		    $biological_first_match_index{$local_rank}{$local_fields[$SOURCE_VECTOR_COLUMN]} = $local_line_index;
		}
		$biological_last_match_index{$local_rank}{$local_fields[$SOURCE_VECTOR_COLUMN]} = $local_line_index;
                $prv_source_vector = ($local_fields[$SOURCE_VECTOR_COLUMN]);
		$local_line_index++;
	    }
	    close(GENOME);
	}
	else {
	    printf STDERR "Warning: The rank $local_rank is not recognized\n";
	}
    }
    close(FILES_BY_RANK);
}

# Subroutine: recognize_rank()
# Synopsis: takes a candidate rank as its only argument;
#           returns 1 if the rank is recognized and 0 otherwise.
#
# Args: $rank
#
# Returns: 1 if the rank is recognized and 0 otherwise

sub recognize_rank {

    my $sub_name = "recognize_rank()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_rank) = @_;
    my $local_r;
    
    for($local_r = 0; $local_r < $num_taxonomy_ranks; $local_r++) {
	if ($local_rank eq $taxonomy_ranks[$local_r]) { 
	    return 1;
	}
    }
    return 0;
}

# Subroutine: find_overlap()
# Synopsis: find the amount of overlap between two intervals, specified by four endpoints.
#
# Args: $local_left_end1, $local_right_end1, $local_left_end2, $local_right_end2
#
# Returns: amount of overlap

sub find_overlap {
    my $sub_name = "find_overlap()";
    my $nargs_exp = 4;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_left_end1, $local_right_end1, $local_left_end2, $local_right_end2) = @_;

    my $local_overlap_amount;  #amount of overlap between the two intervals
    my $local_intersect_left;  #the left endpoint that is further to the right
    my $local_intersect_right; #the right endpoint that is further to the left

    if ($local_right_end1 < $local_left_end2) {
        return 0;
    }
    if ($local_right_end2 < $local_left_end1) {
        return 0;
    }
    if ($local_left_end1 <= $local_left_end2) {
        $local_intersect_left = $local_left_end2;
    }
    else {
        $local_intersect_left = $local_left_end1;
    }
    if ($local_right_end1 <= $local_right_end2) {
        $local_intersect_right = $local_right_end1;
    }
    else {
        $local_intersect_right = $local_right_end2;
    }
    $local_overlap_amount = $local_intersect_right - $local_intersect_left + 1;
    return $local_overlap_amount;
}

# Subroutine: find_ancestor()
# Synopsis: given a taxon, returns its ancestor that is at a specific rank such as "genus";
#           it returns 1 (for root) if there is no ancestor at the specified target rank
#
# Args: $local_taxon
#       $local_target_rank
#
# Returns: the ancestor taxon that is at the specified rank or 1 (for root)

sub find_ancestor {

    my $sub_name = "find_ancestor()";
    my $nargs_exp = 2;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_taxon, $local_target_rank) = @_;
    my $local_ancestor; #one line of taxonomy information

    if (!defined($taxonomy_parent{$local_taxon})) {
        print STDERR "Taxid $local_taxon is not in the taxonomy tree, the tree likely needs to be updated\n";
	exit;
    }
    $local_ancestor = $local_taxon;
    while ((1 != $local_ancestor) && ((!(defined($rank_hash{$local_ancestor}))) || (!($rank_hash{$local_ancestor} eq $local_target_rank)))) {
        $local_ancestor = $taxonomy_parent{$local_ancestor};
    }            
    return($local_ancestor);
}        

# Subroutine: find_lca ()
# Synopsis: finds the lowest common ancestor of two integer taxids.
#
# Args: $local_taxid1, $local_taxid2
#
# Returns: lowest common ancestor of the two taxids

sub find_lca {
    my $sub_name = "find_lca()";
    my $nargs_exp = 2;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_taxid1, $local_taxid2) = @_;
    my $local_ancestor1 = $local_taxid1;    #ancestor of local_taxid1
    my $local_ancestor2 = $local_taxid2;
    while ($taxonomy_level{$local_ancestor1} < $taxonomy_level{$local_ancestor2}) {
	$local_ancestor2 = $taxonomy_parent{$local_ancestor2};
    }
    while ($taxonomy_level{$local_ancestor1} > $taxonomy_level{$local_ancestor2}) {
	$local_ancestor1 = $taxonomy_parent{$local_ancestor1};
    }
    while ($local_ancestor1 != $local_ancestor2) {
	$local_ancestor1 = $taxonomy_parent{$local_ancestor1};
	$local_ancestor2 = $taxonomy_parent{$local_ancestor2};
    }
    return($local_ancestor1);
}
