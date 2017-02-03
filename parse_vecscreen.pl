#!/usr/bin/perl -w
# the first line of perl code has to be above
#
# Author: Alejandro Schaffer with help from Eric Nawrocki
# Code to parse vecscreen output for vector screening.
#
# Usage: parse_vecscreen.pl --input <input file>                 \ [REQUIRED]
#                           --outfile_internal <output internal> \ [REQUIRED]
#                           --outfile_terminal <output terminal> \ [REQUIRED]
#                           --verbose                            \ [OPTIONAL]
#                           --debug <query-for-debugging>          [OPTIONAL]
#
use strict;
use warnings;
use Getopt::Long;

require "epn-options.pm";
  
# variables related to command line options
my $input_file;           # input file of vecscreen output
my $output_internal_file; # output file for sequences from Genbank with only internal matches
my $output_terminal_file; # output file for sequences from Genbank with only terminal matches
my $verbose_mode;         # did user specify verbose mode, in which extra columns are printed
my $debug_mode;           # did user specify debugging mode
my $debug_query = undef;  # query accession for debugging

# hard-coded vecscreen constants, these are treated as global variables
my $DISTANCE_FOR_TERMINAL         = 25; # number of nucleotides to the end for a match to be considered terminal
my $STRONG_SCORE_LOWER_TERMINAL   = 24; # lower score threshold for a terminal match to be a strong match
my $STRONG_SCORE_LOWER_INTERNAL   = 30; # lower score threshold for an internal match to be a strong match
my $MODERATE_SCORE_UPPER_TERMINAL = 23; # upper score threshold for a terminal match to be a moderate match
my $MODERATE_SCORE_UPPER_INTERNAL = 29; # upper score threshold for an internal match to be a moderate match
my $MODERATE_SCORE_LOWER_TERMINAL = 19; # lower score threshold for a terminal match to be a moderate match
my $MODERATE_SCORE_LOWER_INTERNAL = 25; # lower score threshold for an internal match to be a moderate match
my $WEAK_SCORE_LOWER_TERMINAL     = 16; # lower score threshold for a terminal match to be a weak match
my $WEAK_SCORE_LOWER_INTERNAL     = 23; # lower score threshold for an internal match to be a weak match

# different possible states for the finite automaton used to help parse the input file
my $State_Naive           = 0;
my $State_FoundBLASTN     = 1;
my $State_FoundQuery      = 2;
my $State_FoundAtLeastOne = 3;
my $State_FoundMatch      = 4;

# different line types in the input file
my $Linetype_BLASTN      = 0;
my $Linetype_match       = 1;
my $Linetype_position    = 2;
my $Linetype_query_name  = 3;
my $Linetype_length      = 4;
my $Linetype_vector_name = 5;
my $Linetype_score       = 6;
my $Linetype_query_aln   = 7;
my $Linetype_vector_aln  = 8;
my $Linetype_other       = 9;

# variables used while parsing input file
my $i;                     # loop index
my $state;                 # current state of the finite automaton
my $linetype;              # type of current line
my $next_score;            # score of the next match
my $line;                  # one line of input file
my $match_adjective;       # Strong, Moderate, Weak, or Suspect
my $query_name;            # the identifier for the query
my $vector_name;           # string identifying a vector that matches a query
my $final_score;           # score of the current match
my $suspect_output_string; # string for suspect origins
my $print_flag;            # '0' or '1' valued variable that determines whether we print a match or not
my $startPos;              # starting position of a match to vector
my $endPos;                # ending position of a match to vector
my $query_length;          # length of current query
my $seen_vector_name;      # have we seen the vector identifier for this match or not
my $vector_start;          # first aligned position in one block of the vector (subject) sequence
my $vector_end;            # last aligned position in one block of the vector (subject) sequence
my $query_start;           # first aligned position in one block of the query sequence
my $query_end;             # last aligned position in one block of the query  sequence
my $processing_alignment;  # are we possibly in the middle of processing an alignment
my $overall_vector_start;  # vector start for an alignment that can have multiple blocks
my $overall_vector_end;    # vector end for an alignment that can have multiple blocks
my $overall_query_start;   # query start for an alignment that can have multiple blocks
my $overall_query_end;     # query end for an alignment that can have multiple blocks
my $has_strong_match;      # does the query have a strong match
my $has_moderate_match;    # does the query have a moderate match
my $has_weak_match;        # does the query have a weak match
my $output_line;           # a line of output
my $internal_match;        # is this match an internal match

my $alignments_one_query;        # how many alignments have we seen for one query
my $alignments_one_matched_pair; # how many alignments have we seen for one query-vector pair

# file handles
my $terminal_FH;  # output filehandle for terminal matches
my $internal_FH;  # output filehandle for internal matches

# variables related to command line options, see epn-options.pm
my %opt_HH = ();
my @opt_order_A = ();
my %opt_group_desc_H = ();

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
# The opt_Add() function is the way we add options to %opt_HH.
# It takes values of for each of the 2nd dim keys listed above.
#       option                  type       default group   requires incompat  preamble-outfile                             help-outfile
opt_Add("-h",                   "boolean", 0,          0,    undef, undef,    undef,                                       "display this help",  \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"1"} = "required options";
opt_Add("--input",              "string",  undef,      1,    undef, undef,   "input fasta file",                           "REQUIRED: input file name <s> with vecscreen output",           \%opt_HH, \@opt_order_A);
opt_Add("--outfile_internal",   "string",  undef,      1,    undef, undef,   "output of sequences with internal matches",  "REQUIRED: output file <s> for sequences with internal matches", \%opt_HH, \@opt_order_A);
opt_Add("--outfile_terminal",   "string",  undef,      1,    undef, undef,   "output of sequences with terminal matches",  "REQUIRED: output file <s> for sequences with terminal matches", \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"2"} = "other options (not required)";
opt_Add("--verbose",            "boolean", 0,          2,    undef, undef,   "be verbose",                                 "be verbose; output commands to stdout as they're run",          \%opt_HH, \@opt_order_A);
opt_Add("--debug",              "string",  0,          2,    undef, undef,   "turn debugging mode on for query <s>",       "turn debugging mode on for query <s>",                          \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $all_options_recognized =
    &GetOptions(
    'h'                  => \$GetOptions_H{"-h"},
    'input=s'            => \$GetOptions_H{"--input"},
    'outfile_internal=s' => \$GetOptions_H{"--outfile_internal"},
    'outfile_terminal=s' => \$GetOptions_H{"--outfile_terminal"},
    'verbose'            => \$GetOptions_H{"--verbose"},
    'debug'              => \$GetOptions_H{"--debug"});


my $synopsis = "parse_vecscreen.pl :: convert vecscreen output file to one or more tab-delimited output files\n";
my $usage    = "Usage: perl parse_vecscreen.pl ";

my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.01";
my $releasedate   = "Jan 2017";

# set options in %opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

# define file names
$input_file           = opt_Get("--input",           \%opt_HH);
$output_internal_file = opt_Get("--outfile_internal", \%opt_HH);
$output_terminal_file = opt_Get("--outfile_terminal", \%opt_HH);
$verbose_mode         = opt_Get("--verbose",         \%opt_HH);
$debug_mode           = opt_IsUsed("--debug",        \%opt_HH);
if($debug_mode) { 
  $debug_query = opt_Get("--debug", \%opt_HH);
}

# Die if any of:
# - non-existent option is used
# - any of the required options are not used. 
# - -h is used
my $reqopts_errmsg = "";
if(! defined $input_file)           { $reqopts_errmsg .= "ERROR, --input_file not used. It is required.\n"; }
if(! defined $output_internal_file) { $reqopts_errmsg .= "ERROR, --output_internal option not used. It is required.\n"; }
if(! defined $output_terminal_file) { $reqopts_errmsg .= "ERROR, --output_terminal option not used. It is required.\n"; }
if(($GetOptions_H{"-h"}) || ($reqopts_errmsg ne "") || (! $all_options_recognized)) { 
    opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
    if($GetOptions_H{"-h"})          { exit 0; } # -h, exit with 0 status
    elsif($reqopts_errmsg ne "")     { die $reqopts_errmsg; }
    elsif(! $all_options_recognized) { die "ERROR, unrecognized option;"; } 
}

# open the input and output files
open(INPUT,        "<", $input_file)           or die "Cannot open $input_file for reading\n"; 
open($internal_FH, ">", $output_internal_file) or die "Cannot open $output_internal_file for writing\n"; 
open($terminal_FH, ">", $output_terminal_file) or die "Cannot open $output_terminal_file for writing\n"; 

###################################################################
# Parse the input file (the vecscreen output).
#
# The following line types need to be recognized and distinguished:
#
#   A. BLASTN lines:            lines that begin with "BLASTN", indicating the start of output for a new query 
#                               isBLASTNLine() returns '1' 
#                               1 per query
#
#   B. Match lines:             lines that begin with "Strong", "Moderate", "Weak" or "Suspect" indicating a type of match
#                               getMatch() returns name of match (!0)
#                               >= 0 per query
#
#   C. Position lines:          lines that determine whether a query has a terminal match or has an internal match
#                               isPositionLine() returns '1'
#                               getPositions() returns $start and $stop
#                               >= 0 per query
#
#   D. Query name lines:        lines that begin with "Query=", including the name of a query
#                               isQueryNameLine() returns '1'
#                               1 per query
#
#   E. Length lines:            lines that begin with "Length", including the length of the query or vector (subject)
#                               getLength() returns length (!0)
#                               These lines are *only* relevant when they pertain to the query
#                               which is only true when $seen_vector_name is 0
#                               >=1 per query, only one of which is relevant (length of query)
#
#   F. Vector name lines:       lines that begin with "> gnl"   including the name of a vector (subject)
#                               isVectorNameLine() returns '1'
#                               >= 0 per query
#
#   G. Score lines:             lines that contain "Score" including score for a match
#                               getScore() returns score (!0)
#                               >= 0 per query
#
#   H. Query alignment lines:   lines that begin with "Query "  including a segment of the query in an alignment
#                               isQueryAlnmentLine() returns '1'
#                               >= 0 per query
#
#   I. Vector alignment lines:  lines that begin with "Sbjct "  including a segment of the vector (subject) in an alignment
#                               isVectorAlnmentLine() returns '1'
#                               >= 0 per query
#
# - these types were named based on order they appear in the file, A
#   is first, B second, etc.  
#
# - there is not exactly 1 line of each type per each query, some have
#    0, others > 1, and many line type counts are query-dependent.
#
# The determineLineType() subroutine determines the line type of 
# each line and updates the $linetype variable to one of 10 values
# (A through I is 9, and $Linetype_other is the 10th).
#
# Below is a table with 5 columns that shows an example input file
# (column 5) and what type each line is (column 2). Note that some
# lines have no type (are ignored). Column 3 shows the subroutine that
# is most relevant to recognizing each line type.  The value of the
# $state variable is shown as what it will be after each line in
# column 4.  The file line numbers (column 1) are sometimes referred
# to in the comments within the code of the large while() loop below
# that loops over lines of the input file.
# 
# Of course, this simple example does not cover all possibilities.
# For example, there can be more than one match between the 
# same query/vector pair, which is not demonstrated in this
# example. 
#
#----------------------------------------------------------------------------------------------------------------------------------------
#     |      | value of $linetype      |                          |
#     |      | after each line, set by |                          |
#file |      | determineLineType()     |                          |
#line | line | (blank=>$Linetype_other | value of                 | example file input line 
# num | type |  which are not parsed)  | $state after each line   | value of $line, read from <INPUT> at top of loop
#-----|--------------------------------|--------------------------|-----------------------------------------------------------------------
#   1 |    A | ->$Linetype_vector_aln  | ->$State_FoundBLASTN     |BLASTN 2.6.0+
#   2 |      |                         | "                        |
#   3 |      |                         | "                        |
#   4 |      |                         | "                        | Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.
#   5 |      |                         | "                        | Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.
#   6 |      |                         | "                        | Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of
#   7 |      |                         | "                        | protein database search programs", Nucleic Acids Res. 25:3389-3402.
#   8 |      |                         | "                        |
#   9 |      |                         | "                        |
#  10 |      |                         | "                        |
#  11 |      |                         | "                        |Database: UniVec (build 9.0)
#  12 |      |                         | "                        |           5,456 sequences; 1,049,913 total letters
#  13 |      |                         | "                        |
#  14 |    B | ->$Linetype_match       | ->$State_FoundMatch      |Moderate match
#  15 |    C | ->$Linetype_position    | ->$State_FoundAtLeastOne |1285	1333
#  16 |    C | ->$Linetype_position    | "                        |1408	1436
#  17 |    B | ->$Linetype_match       | "                        |Suspect origin
#  18 |    C | ->$Linetype_position    | "                        |1437	1467
#  19 |      |                         | "                        |
#  20 |      |                         | "                        |
#  21 |    D | ->$Linetype_query_aln   | ->$State_FoundQuery      |Query= XM_715279.1 Candida albicans SC5314 Spl1p (SPL1), partial mRNA
#  22 |      |                         | "                        |
#  23 |    E | ->$Linetype_length      | "                        |Length=1467
#  24 |      |                         | "                        |
#  25 |      |                         | "                        |
#  26 |    F | ->$Linetype_vector_name | "                        |> gnl|uv|U29899.1:1847-4463 Cloning vector pACT2 MatchmakerII
#  27 |    E | ->$Linetype_length      | "                        |Length=2617
#  28 |      |                         | "                        |
#  29 |    G | ->$Linetype_score       | "                        | Score = 58.6 bits (29),  Expect = 4e-06
#  30 |      |                         | "                        | Identities = 29/29 (100%), Gaps = 0/29 (0%)
#  31 |      |                         | "                        | Strand=Plus/Minus
#  32 |      |                         | "                        |
#  33 |    H | ->$Linetype_query_aln   | "                        | Query  1408  TTATGGGAAATGGTTCAAGAAGGTATTGA  1436
#  34 |      |                         | "                        |              |||||||||||||||||||||||||||||
#  35 |    I | ->$Linetype_vector_aln  | "                        | Sbjct  2140  TTATGGGAAATGGTTCAAGAAGGTATTGA  2112
#  36 |      |                         | "                        |
#  37 |      |                         | "                        |
#  38 |    F | ->$Linetype_vector_name | "                        |> gnl|uv|U03498.1:7366-8059 Yeast episomal vector YEp13
#  39 |    E | ->$Linetype_length      | "                        |Length=694
#  40 |      |                         | "                        |
#  41 |    G | ->$Linetype_score       | "                        | Score = 50.6 bits (25),  Expect = 0.001
#  42 |      |                         | "                        | Identities = 45/49 (92%), Gaps = 0/49 (0%)
#  43 |      |                         | "                        | Strand=Plus/Minus
#  44 |      |                         | "                        | 
#  45 |    H | ->$Linetype_query_aln   | "                        |Query  1285  GATGATGCCTTAGCTCATTCTTCAATTAGATTTGGTATTGGTAGATTTA  1333
#  46 |      |                         | "                        |             |||||||| ||||| |||||||| || ||||||||||||||||||||||
#  47 |    I | ->$Linetype_vector_aln  | "                        |Sbjct  113   GATGATGCATTAGCCCATTCTTCCATCAGATTTGGTATTGGTAGATTTA  65
#  48 |      |                         | "                        |
#  49 |      |                         | "                        |
#  50 |      |                         | "                        |
#  51 |      |                         | "                        |  Lambda      K        H
#  52 |      |                         | "                        |    1.39    0.747     1.38 
#  53 |      |                         | "                        | 
#  54 |      |                         | "                        | Gapped
#  55 |      |                         | "                        | Lambda      K        H
#  56 |      |                         | "                        |   1.39    0.747     1.38 
#  57 |      |                         | "                        |
#  58 |      |                         | "                        | Effective search space used: 1750000000000
#  59 |      |                         | "                        |
#  60 |      |                         | "                        |
#  61 |      |                         | "                        | Database: UniVec (build 9.0)
#  62 |      |                         | "                        |  Posted date:  Mar 23, 2015  3:44 PM
#  63 |      |                         | "                        |  Number of letters in database: 1,049,913
#  64 |      |                         | "                        |  Number of sequences in database:  5,456
#  66 |      |                         | "                        |
#  67 |      |                         | "                        |
#  68 |      |                         | "                        |
#  69 |      |                         | "                        | Matrix: blastn matrix 1 -5
#  70 |      |                         | "                        | Gap Penalties: Existence: 3, Extension: 3
############################################################################################
#
# Output: 
# Each query/target match is transformed into a single line of tabular output
# and output when that match is finished being processed from the input.
# There are 3 different situations in which we output a single tabular 
# line describing the most recent parsed match:
#
# output situation 1: output the final query/vector match for a query. 
#                     Occurs only when (i) we see the 'next' BLASTN line (line type A)
#                     or (ii) we run out of input lines because only then do we 
#                     know we are finished with a query.
#
# output situation 2: output a query/vector match n for a query/vector pair with 
#                     X (X>1) matches, where n < X. Occurs when we see the next 'Score' line
#                     (line type G) indicating a new match for the current query/vector 
#                     pair is coming up, so we output the previous one.
#
# output situation 3: output a query/vector match n for a query/vector pair with 
#                     X (X>=1) matches, where n == X when this vector is not the 
#                     final vector to have >=1 matches to this query (if it was, it 
#                     would be covered by situation 1). Occurs when we see the next VectorName
#                     line (line type G) indicating a new vector with at least one match to 
#                     the current query.
#
#

# Initializations
$state                = $State_Naive;
$seen_vector_name       = 0;
$processing_alignment = 0;

# For each line what we do depends strictly on the line type.
# The following while() loop has the following structure:
# if   (line type BLASTN)      { }
# elsif(line type match)       { }
# elsif(line type position)    { }
# elsif(line type query_name)  { }
# elsif(line type length)      { }
# elsif(line type vector_name) { }
# elsif(line type score)       { }
# elsif(line type query_aln)   { }
# elsif(line type vector_aln)  { }
# 
while(defined($line = <INPUT>)) { 
  chomp($line);
  $linetype = determineLineType($line, 
                                $Linetype_BLASTN,
                                $Linetype_match,
                                $Linetype_position,
                                $Linetype_query_name,
                                $Linetype_length,
                                $Linetype_vector_name,
                                $Linetype_score,
                                $Linetype_query_aln,
                                $Linetype_vector_aln,
                                $Linetype_other);

  ######################################################################################
  if($linetype == $Linetype_BLASTN) { 
    # A. BLASTN line: line that begins with "BLASTN", indicating the start of output 
    #                 for a new query 
    if ($print_flag) {
      # output situation 1: output the final query/vector match for a query. 
      #                     (sub-situation (i): we've seen the 'next' BLASTN line
      ($internal_match, $output_line) = 
          createOutputLine($verbose_mode, $query_name, $overall_query_start, $overall_query_end, $query_length, 
                           $vector_name, $overall_vector_start, $overall_vector_end, $has_strong_match,
                           $has_moderate_match, $has_weak_match, $final_score, $suspect_output_string);
      printOutputLine($output_line, $internal_FH, $terminal_FH, $internal_match, $debug_mode, $debug_query, $query_name, 1, $line);
      # 1 is calling 'situation' only used if in debug mode
    }
    # reset variables for the next query 
    $suspect_output_string       = "";
    $alignments_one_query        = 0;
    $alignments_one_matched_pair = 0;
    $print_flag                  = 0; 
    $state                       = $State_FoundBLASTN;
    $seen_vector_name            = 0;	
    $processing_alignment        = 0;
    $has_strong_match = $has_moderate_match = $has_weak_match = 0;
    if ($debug_mode && defined($query_name) && ($query_name =~m/$debug_query/)) {
      print "#DEBUG: set print_flag to $print_flag\n";
    }
  }
  ######################################################################################
  elsif($linetype == $Linetype_match) { 
    # B. Match lines: lines that begin with "Strong", "Moderate", "Weak" or "Suspect" indicating a type of match
    $match_adjective = getMatch($line);
    updateMatchAdjectives($match_adjective, \$has_strong_match, \$has_moderate_match, \$has_weak_match);
    if($state == $State_FoundBLASTN) { 
      # First match line (because $state == $State_FoundBLASTN) for this query
      # *** line num 14 in example file above ***
      $state = $State_FoundMatch;
    }
    elsif(($state == $State_FoundMatch) || ($state == $State_FoundAtLeastOne)) { 
      # NOT first match line for this query (because $state == $State_FoundMatch or $state == $State_FoundAtLeastOne)
      # *** line num 17 in example file above ***
      $state = $State_FoundAtLeastOne;
      $print_flag = 1;
      if ($debug_mode && defined($query_name) && ($query_name =~m/$debug_query/)) {
        print "#DEBUG: setting print_flag to 1 at location B; triggering line is $line\n"; 
      }
    }
  }
  ######################################################################################
  elsif($linetype == $Linetype_position) { 
    # C. Position line: line that determines whether a query has a terminal match or has an internal match
    if(($state == $State_FoundMatch) || ($state == $State_FoundAtLeastOne)) { 
      # *** line num 15 in above example file *** 
      ($startPos, $endPos) = getPositions($line);
      $state = $State_FoundAtLeastOne;
      $print_flag = 1; # we have a match, so we need to print it eventually
      if ($match_adjective eq "Suspect") { 
        # *** line num 18 in above example file ***
        # last seen match line was of type 'suspect'
        # update suspect_output_string which holds suspect regions
        $suspect_output_string .= $match_adjective . '[' . $startPos . ',' . $endPos . ']' . ';' ;
      }
      if ($debug_mode && defined($query_name) && ($query_name =~m/$debug_query/)) {
        print "#DEBUG: setting print_flag to 1 at location A; triggering line is $line\n"; 
      }
    }
  }
  ######################################################################################
  elsif($linetype == $Linetype_query_name) { 
    if($state == $State_FoundAtLeastOne) { 
      # D. Query name line: line that begins with "Query=", including the name of a query
      # *** line num 21 in example file above ***
      # by requiring that $state == $State_FoundAtLeastOne, we skip
      # queries that have zero matches
      $query_name = getQueryName($line);
      $state      = $State_FoundQuery;
    }
  }
  ######################################################################################
  elsif($linetype == $Linetype_length) { 
    if(! $seen_vector_name) { 
      # E. Length line: line that begins with "Length", including the length of the query or vector (subject)
      # *** line num 23 in above example file ***
      # the (!$seen_vector_name) part ensures this is the query length and not vector(subject) length
      $query_length = getLength($line);
    }
  }
  ######################################################################################
  elsif($linetype == $Linetype_vector_name) { 
    #F. Vector name line: line that begins with "> gnl" including the name of a vector (subject)
    # *** line num 38 in above example file ***
    if($print_flag) { 
      if ($alignments_one_query > 0) {
        # if we've seen an alignment of this query to a previous vector (subject)
        # (not the same as this vector) then output info on that previous match
        # output situation 3: output a query/vector match n for a query/vector pair with 
        #                     X (X>=1) matches, where n == X when this vector is not the 
        #                     final vector to have >=1 matches to this query (if it was, it 
        #                     would be covered by situation 1). Occurs when we see the next VectorName
        #                     line (line type G) indicating a new vector with at least one match to 
        #                     the current query.
        ($internal_match, $output_line) = 
            createOutputLine($verbose_mode, $query_name, $overall_query_start, $overall_query_end, $query_length, 
                             $vector_name, $overall_vector_start, $overall_vector_end, $has_strong_match,
                             $has_moderate_match, $has_weak_match, $final_score, $suspect_output_string);
        printOutputLine($output_line, $internal_FH, $terminal_FH, $internal_match, $debug_mode, $debug_query, $query_name, 3, $line);
        # reset variables
        $alignments_one_matched_pair = 0;
        $alignments_one_query        = 0;
      }
      # store new vector/vector name, very important not to do this before calling 
      # createOutputLine because we want to output match info for *previous* vector
      $seen_vector_name     = 1;
      $vector_name          = getVectorName($line);
      $processing_alignment = 0;
    }
  }
  ######################################################################################
  elsif($linetype == $Linetype_score) { 
    # G. Score line: line that contains "Score" including score for a match
    # *** line num 29 or 41 in above example file ***
    if ($state == $State_FoundQuery) { 
      $next_score = getScore($line);
      if ($alignments_one_matched_pair > 0) {
        # output situation 2: output a query/vector match n for a query/vector pair with 
        #                     X (X>1) matches, where n < X. Occurs when we see the next 'Score' line
        #                     (line type G) indicating a new match for the current query/vector 
        #                     pair is coming up, so we output the previous one.
        ($internal_match, $output_line) = 
            createOutputLine($verbose_mode, $query_name, $overall_query_start, $overall_query_end, $query_length, 
                             $vector_name, $overall_vector_start, $overall_vector_end, $has_strong_match,
                             $has_moderate_match, $has_weak_match, $final_score, $suspect_output_string);
        printOutputLine($output_line, $internal_FH, $terminal_FH, $internal_match, $debug_mode, $debug_query, $query_name, 2, $line);
        $processing_alignment = 0;
      }
      $final_score = $next_score; # important to wait until after the createOutputLine
                                  # call to update this because the output is on the *previous* match
                                  # whose score is in $final_score until we redefine it here.
    }
  }
  ######################################################################################
  elsif($linetype == $Linetype_query_aln) { 
    # H. Query alignment line:line that begins with "Query "including a segment of the query in an alignment
    ($query_start, $query_end) = getAlignmentPositions($line); 
    updateOverallPositions($processing_alignment, $query_start, $query_end, \$overall_query_start, \$overall_query_end);
  } 
  ######################################################################################
  elsif($linetype == $Linetype_vector_aln) { 
    # I. Vector alignment line: line that begins with "Sbjct "  including a segment of the vector (subject) in an alignment
    ($vector_start, $vector_end) = getAlignmentPositions($line); 
    updateOverallPositions($processing_alignment, $vector_start, $vector_end, \$overall_vector_start, \$overall_vector_end);
    if (! $processing_alignment) {
      # the first Vector Alignment line we've seen for this match
      $processing_alignment = 1;
      $alignments_one_query++;
      $alignments_one_matched_pair++;
    }
  } 
  ######################################################################################
} # end of 'while(defined($line = <INPUT>))'
######################################################################################
# final step: print the final match, if necessary
if ($print_flag) {
  # output situation 1: output the final query/vector match for a query. 
  #                     (sub-situation (ii): we've run out of input lines so we know
  #                     we're now finished with final query.
  ($internal_match, $output_line) = 
      createOutputLine($verbose_mode, $query_name, $overall_query_start, $overall_query_end, $query_length, 
                       $vector_name, $overall_vector_start, $overall_vector_end, $has_strong_match,
                       $has_moderate_match, $has_weak_match, $final_score, $suspect_output_string);
  printOutputLine($output_line, $internal_FH, $terminal_FH, $internal_match, $debug_mode, $debug_query, 4, $line);
}

# close files and exit
close(INPUT);
close($internal_FH);
close($terminal_FH);

exit 0;

################################################
# List of subroutines:
#
# Subroutines to extract needed info from an input file line:
# getMatch():              get type of match from a match line
# getPositions():          get start and end positions from a positions line
# getQueryName():          get name of query from a query name line
# getLength():             get length from a length line
# getVectorName():         get name of vector (subject) from a vector name line 
# getScore():              get score from a score line
# getAlignmentPositions(): get start and end positions from an alignment line 
# 
# Subroutines related to match adjectives:
# determineStrongestAdjective():
# determineMatchAdjective():
# updateMatchAdjectives():
#
# Subroutines related to output:
# createOutputLine():
# printOutputLine():
# 
# Miscellaneous subroutines: 
# determineLineType(): returns line type for input file line
# isInternalMatch():
# updateOverallPositions
#
######################################################################################

##########################################################################################
# Subroutine: getMatch()
# Synopsis:   Extracts the type of match from a match line.
#
# Args: $line: match line of input file
#
# Returns: the type of match, either "Strong", "Moderate", "Weak", or "Suspect"
# 
# Dies: if $line is not a match line;
#       does not start with "Strong", "Moderate", 
#       "Weak" or "Strong".
#
##########################################################################################
sub getMatch {
  my $sub_name = "getMatch()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($line) = @_;
  
  if   ($line =~m/^Strong/)   { return "Strong";   }
  elsif($line =~m/^Moderate/) { return "Moderate"; }
  elsif($line =~m/^Weak/)     { return "Weak";     }
  elsif($line =~m/^Suspect/)  { return "Suspect";  }
  else { die "ERROR in $sub_name, unexpected format of line: $line\n"; }
}

##########################################################################################
# Subroutine: getPositions() 
# Synopsis:   Extracts two numerical positions from
#             a position line.
#
# Args: $line: positions line of input file 
#
# Returns: TWO values: 
#          $start: the start position
#          $end:   the end position
# 
# Dies: if $line is not a position line;
#       does not match /(\d+)\s+(\d+)/
#
##########################################################################################
sub getPositions {
  my $sub_name = "getPositions()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($line) = @_;

  if($line =~m/(\d+)\s+(\d+)/) { 
    return ($1, $2);
  }
  else { 
    die "ERROR in $sub_name, unexpected format of line: $line\n";
  }
}

##########################################################################################
# Subroutine: getQueryName() 
# Synopsis:   Extracts the name of a query from 
#             a query name line.
#
# Args: $line: query name line of input file
#
# Returns: name of the query
# 
# Dies: if $line is not a query name;
#       does not match
#       /^Query=\s+\S+/;
#
##########################################################################################
sub getQueryName {
  my $sub_name = "getQueryName()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($line) = @_;

  if($line =~m/^Query=\s+(\S+)/) { 
    return $1;
  }
  else { 
    die "ERROR in $sub_name, unexpected format of line: $line\n";
  }
}

##########################################################################################
# Subroutine: getLength()
# Synopsis:   Extracts the length from a 
#             length line.
#
# Args: $line: length line of input file
#
# Returns: length
# 
# Dies: if $line is not a length name;
#       does not match
#       /^Length=(\d+)/.
#
##########################################################################################
sub getLength { 
  my $sub_name = "getLength()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($line) = @_;

  if($line =~m/^Length=(\d+)/) { 
    return $1;
  }
  else { 
    die "ERROR in $sub_name, unexpected format of line: $line\n";
  }
}

##########################################################################################
# Subroutine: getVectorName()
# Synopsis:   Extracts the vector name from 
#             a vector name line.
#
# Args: $line: vector name line of input file
#
# Returns: name of the vector (subject)
# 
# Dies: if $line is not a vector name;
#       does not match
#       /^>\s*gnl\|(\S+)/.
#
##########################################################################################
sub getVectorName { 
  my $sub_name = "getVectorName()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($line) = @_;

  if($line =~m/^>\s*gnl\|(\S+)/) { 
    return $1;
  }
  else { 
    die "ERROR in $sub_name, unexpected format of line: $line\n";
  }
}

##########################################################################################
# Subroutine: getScore()
# Synopsis:   Extracts the score from a 
#             score line.
#
# Args: $line: score line of input file
#
# Returns: score
# 
# Dies: if $line is not a score name;
#       does not match
#       /^Length=(\d+)/.
#
##########################################################################################
sub getScore { 
  my $sub_name = "getScore()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($line) = @_;

  if($line =~ m/^\s*Score\s*=\s+\S+\s+bits\s+\((\d+)\)/) { 
    return $1;
  }
  else { 
    die "ERROR in $sub_name, unexpected format of line: $line\n";
  }
}

##########################################################################################
# Subroutine: getAlignmentPositions() 
# Synopsis:   Extracts two numerical positions from
#             a query or vector alignment line.
#
# Args: $line: query alignment or vector alignment line of 
#              input file 
#
# Returns: TWO values: 
#          $start: the start position
#          $end:   the end position
# 
# Dies: if $line is not a query or vector alignment line;
#       does not match either of 
#       /Query\s+(\d+)\s+\D+(\d+)/ or
#       /Sbjct\s+(\d+)\s+\D+(\d+)/ 
#
##########################################################################################
sub getAlignmentPositions {
  my $sub_name = "getAlignmentPositions()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($line) = @_;
  
  if($line =~ m/^Query\s+(\d+)\s+\D+(\d+)/) { 
    return ($1, $2);
  }
  elsif($line =~ m/^Sbjct\s+(\d+)\s+\D+(\d+)/) { 
    return ($1, $2);
  }
  else { 
    die "in $sub_name(), unexpected format for input line: $line\n";
  }
}

##########################################################################################
# Subroutine: determineStrongestAdjective() 
# Synopsis:   Determine the strongest match adjective out 
#             of Strong, Moderate, Weak and return it.
#
# Args:       
# $has_strong:   '1' if query has at least one 'Strong' match
# $has_moderate: '1' if query has at least one 'Moderate' match
# $has_weak:     '1' if query has at least one 'Weak' match
#
# Returns: The adjective describing the strongest vector
#          match for this query. "None" if all 3 input
#          arguments are '0'.
#
# Dies:    Never.
##########################################################################################
sub determineStrongestAdjective {
  my $sub_name = "determineStrengthAdjective()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($has_strong, $has_moderate, $has_weak) = @_;

  if    ($has_strong)   { return "Strong";   }
  elsif ($has_moderate) { return "Moderate"; }
  elsif ($has_weak)     { return "Weak";     }
  else                  { return "None";     }
}

##########################################################################################
# Subroutine: determineMatchAdjective() 
# Synopsis:   Determine the adjective for one match out of 
#             Strong, Moderate, Weak and return it.
#
# Args: 
# $score:        score of current match
# $is_internal:  '1' if match is an internal match, else '0'
#             
# Returns: the adjective ("Strong", "Moderate" or "Weak") describing 
#          the quality of this match. If score is below the "Seak"
#          threshold, then return "None".
# 
# Dies:  Never.
##########################################################################################
sub determineMatchAdjective {
  my $sub_name = "determineMatchAdjective()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($score, $is_internal) = @_;

  if ($is_internal) {
    if   ($score >= $STRONG_SCORE_LOWER_INTERNAL)   { return "Strong";   }
    elsif($score >= $MODERATE_SCORE_LOWER_INTERNAL) { return "Moderate"; }
    elsif($score >= $WEAK_SCORE_LOWER_INTERNAL)     { return "Weak";     }
    else                                            { return "None";     }
  }
  else { # terminal match
    if   ($score >= $STRONG_SCORE_LOWER_TERMINAL)   { return "Strong";   }
    elsif($score >= $MODERATE_SCORE_LOWER_TERMINAL) { return "Moderate"; }
    elsif($score >= $WEAK_SCORE_LOWER_TERMINAL)     { return "Weak";     }
    else                                            { return "None";     }
  }
}

##########################################################################################
# Subroutine: updateMatchAdjectives() 
# Synopsis:   Update variables which indicate the severities 
#             of matches one query has.
#
# Args: 
# $adjective:      adjective describing current match
# $has_strong_R:   reference to variable storing whether query has a strong match or not
# $has_moderate_R: reference to variable storing whether query has a moderate match or not
# $has_weak_R:     reference to variable storing whether query has a weak match or not
#             
# Returns: Void. Updates $$has_strong_R, $$has_moderate_R, $$has_weak_R.
# 
# Dies:  If $adjective is none of "Strong", "Moderate", "Weak" or "Suspect".
#
##########################################################################################
sub updateMatchAdjectives {
  my $sub_name = "updateMatchAdjectives()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($adjective, $has_strong_R, $has_moderate_R, $has_weak_R) = @_;
  if    ($adjective eq "Strong")   { $$has_strong_R   = 1; }
  elsif ($adjective eq "Moderate") { $$has_moderate_R = 1; }
  elsif ($adjective eq "Weak")     { $$has_weak_R     = 1; }
  elsif ($adjective eq "Suspect")  { ; } # do nothing 
  else  { die "ERROR in $sub_name, unexpected adjective $adjective"; }
}

##########################################################################################
# Subroutine: createOutputLine()
#
# Synopsis: Format a line of output given all required 
#           information.
#
# Arguments: 
# $verbose_mode:       '1' if we are in verbose mode, else '0'
# $query_name:          name of the query we are printing
# $query_start:        start position of full alignment in query
# $query_end:          end position of full alignment in query
# $query_length:       total length of query
# $vector:             name of the vector (vector) we are printing
# $vector_start:       start position of full alignment in vector (vector)
# $vector_end:         end position of full alignment in vector (vector)
# $has_strong_match:   '1' if query has a strong match, else '0'
# $has_moderate_match: '1' if query has a moderate match, else '0'
# $has_weak_match:     '1' if query has a weak match, else '0'
# $final_score:        score of match 
# $suspect_string:     string for suspect regions, can be "" for none
#
# Returns: TWO values: 
#          value 1: '1' if match is 'internal', '0' if match is 'terminal'
#          value 2: string to output for this match
#
# Dies: Never.
##########################################################################################
sub createOutputLine { 
  my $sub_name = "createOutputLine()";
  my $nargs_exp = 13;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($verbose_mode, $query_name, $query_start, $query_end, $query_length, 
      $vector, $vector_start, $vector_end, 
      $has_strong_match, $has_moderate_match, $has_weak_match,
      $final_score, $suspect_string) = @_;
  
  my $overall_adjective;    # what is the strongest adjective for which this query has a match of that adjective
  my $one_match_adjective;  # what is the adjective for a single alignment
  my $internal_match;       # is this match an internal match

  $internal_match      = (isInternalMatch($overall_query_start, $overall_query_end, $query_length)) ? 1 : 0;
  $overall_adjective   = determineStrongestAdjective($has_strong_match, $has_moderate_match, $has_weak_match); 
  $one_match_adjective = determineMatchAdjective($final_score, $internal_match);

  if($suspect_string eq "") { $suspect_string = "None"; }

  if ($verbose_mode) {
    return ($internal_match, "$query_name\t$query_start\t$query_end\t$vector\t$vector_start\t$vector_end\t$one_match_adjective\t$overall_adjective\t$suspect_string\n");
  }
  else {
    return ($internal_match, "$query_name\t$vector\t$vector_start\t$vector_end\n");
  }
}

##########################################################################################
# Subroutine: printOutputLine()
#
# Synopsis: Print a line of output to either terminal 
#           or internal output file.
#
# Arguments: 
# $output_line:       line of output to print
# $internal_FH:       file handle for outputting internal match lines
# $terminal_FH:       file handle for outputting terminal match lines
# $internal_match:    '1' if this is an internal match we are printing
#                     '0' if it is a terminal match
# $debug_mode:        '1' if we're in debugging mode
# $debug_query:       query to output information for if $debug_mode is '1'
# $query:             name of query, only relevant if $debug_mode is 1
# $calling_situation: '1', '2', '3', or '4', indicating which of the four
#                     situations this subroutine was called for             
# $line:              triggering line (only used if $debug_mode)
#             
# Returns: Void.
#
# Dies: Never.
##########################################################################################
sub printOutputLine { 
  my $sub_name = "printOutputLine()";
  my $nargs_exp = 9;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($output_line, $internal_FH, $terminal_FH, $internal_match, $debug_mode, $debug_query, $query, $calling_situation, $line) = @_;

  if ($internal_match) { print $internal_FH $output_line; }
  else                 { print $terminal_FH $output_line; }

  if ($debug_mode && $query =~m/$debug_query/) {
    print "#DEBUG: in $sub_name, location $calling_situation; internal match is $internal_match;\n";
    print "#DEBUG: print_flag is $print_flag; triggering line is: $line\n";
    print "#DEBUG: output line is: $output_line";
  }

  return;
}

##########################################################################################
# Subroutine: determineLineType()
#
# Synopsis: Returns type of line $line.
#
# Args: $line:                 the line of input file
#       $Linetype_BLASTN:      '0' (hard-coded 'line type' value for BLASTN lines
#       $Linetype_match:       '1' (hard-coded 'line type' value for match lines
#       $Linetype_position:    '2' (hard-coded 'line type' value for position lines
#       $Linetype_query_name:   '3' (hard-coded 'line type' value for query name lines
#       $Linetype_length:      '4' (hard-coded 'line type' value for length lines
#       $Linetype_vector_name:  '5' (hard-coded 'line type' value for vector name lines
#       $Linetype_score:       '6' (hard-coded 'line type' value for score lines
#       $Linetype_query_aln:    '7' (hard-coded 'line type' value for query alignment lines
#       $Linetype_vector_aln:   '8' (hard-coded 'line type' value for vector alignment lines
#       $Linetype_other:       '9' (hard-coded 'line type' value for other (non-parsed) lines
#
# Returns: line type of $line
# 
# Dies: Never
#
##########################################################################################
sub determineLineType { 
    my $sub_name = "determineLineType()";
    my $nargs_exp = 11;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($line, $Linetype_BLASTN, $Linetype_match, $Linetype_position, $Linetype_query_name,
        $Linetype_length, $Linetype_vector_name, $Linetype_score, $Linetype_query_aln, 
        $Linetype_vector_aln, $Linetype_other) = @_;

    if($line =~ m/^BLASTN/)                              { return $Linetype_BLASTN;      }
    if($line =~ m/^Strong/)                              { return $Linetype_match;       }
    if($line =~ m/^Moderate/)                            { return $Linetype_match;       }
    if($line =~ m/^Weak/)                                { return $Linetype_match;       }
    if($line =~ m/^Suspect/)                             { return $Linetype_match;       } 
    if($line =~ m/^\d/)                                  { return $Linetype_position;    }
    if($line =~ m/^Query=/)                              { return $Linetype_query_name;  }
    if($line =~ m/^Length=\d+/)                          { return $Linetype_length;      }
    if($line =~ m/^>\s*gnl/)                             { return $Linetype_vector_name; }
    if($line =~ m/^\s*Score\s*=\s+\S+\s+bits\s+\(\d+\)/) { return $Linetype_score;       }
    if($line =~ m/^Query\s+\d+\s+\D+\d+/)                { return $Linetype_query_aln;   }
    if($line =~ m/^Sbjct\s+\d+\s+\D+\d+/)                { return $Linetype_vector_aln;  }
    return $Linetype_other; 
}


##########################################################################################
# Subroutine: isInternalMatch() 
# Synopsis:   Return '1' if a query/vector match is deemed to be in 
#             the internal part of the sequence.
#
# Args: 
# $query_start:  start of alignment of match inquery 
# $query_end:    end   of alignment of match inquery 
# $query_length: length of query            
#
# Returns: 1 of the match is internal and 0 otherwise
##########################################################################################
sub isInternalMatch {
  my $sub_name = "isInternalMatch()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($query_start, $query_end, $query_length) = @_;
  if (( $query_start                    > $DISTANCE_FOR_TERMINAL) && 
      (($query_length - $query_end + 1) > $DISTANCE_FOR_TERMINAL)) {
    return 1;
  }
  else {
    return 0;
  }
}

##########################################################################################
# Subroutine:  updateOverallPositions()
#
# Purpose:     Update the $$overall_start_R and $$overall_end_R 
#              which keep track of overall alignment positions
#              (start of full alignment and end of full alignment)
#              based on $start and $end (positions of current
#              alignment block). If ($processing_alignment == 1)
#              then this is not the first block, else it is
#              the first block.
#
# Arguments:
#   $processing_alignment: '1' if not first alignment block
#   $start:                start for current block we are processing
#   $end:                  end  for current block we are processing
#   $overall_start_R:      reference, overall start position
#   $overall_end_R:        reference, overall end position
#
# Returns:    void, but updates $$overall_start_R and $$overall_end_R
#
# Dies:       Never
#
##########################################################################################
sub updateOverallPositions { 
  my $sub_name = "updateOverallPositions()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($processing_alignment, $start, $end, $overall_start_R, $overall_end_R) = @_;

  if ($processing_alignment) { 
    # not the first Query Alignment line we've seen for this match, 
    # update $overall_end, but not $overall_start
    if (1 == abs($start - $$overall_end_R)) { 
      $$overall_end_R = $end;
    }
    # TODO: update this, should we die if this fails?
  }
  else {
    # the first alignment line we've seen for this match
    $$overall_start_R = $start;
    $$overall_end_R   = $end;
  }
  return; 
}
