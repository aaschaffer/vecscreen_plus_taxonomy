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
my $output_file_internal; # output file for sequences from Genbank with only internal matches
my $output_file_terminal; # output file for sequences from Genbank with only terminal matches
my $verbose_mode;         # did user specify verbose mode, in which extra columns are printed
my $debug_mode;           # did user specify debugging mode
my $debug_query = undef;  # query accession for debugging

# hard-coded vecscreen constants
my $DISTANCE_FOR_TERMINAL         = 25; # number of nucleotides to the end for a match to be considered terminal
my $STRONG_SCORE_LOWER_TERMINAL   = 24; # lower score threshold for a terminal match to be a strong match
my $STRONG_SCORE_LOWER_INTERNAL   = 30; # lower score threshold for an internal match to be a strong match
my $MODERATE_SCORE_UPPER_TERMINAL = 23; # upper score threshold for a terminal match to be a moderate match
my $MODERATE_SCORE_UPPER_INTERNAL = 29; # upper score threshold for an internal match to be a moderate match
my $MODERATE_SCORE_LOWER_TERMINAL = 19; # lower score threshold for a terminal match to be a moderate match
my $MODERATE_SCORE_LOWER_INTERNAL = 25; # lower score threshold for an internal match to be a moderate match
my $WEAK_SCORE_LOWER_TERMINAL     = 16; # lower score threshold for a terminal match to be a weak match
my $WEAK_SCORE_LOWER_INTERNAL     = 23; # lower score threshold for an internal match to be a weak match

my $State_Naive           = 0;
my $State_FoundBLASTN     = 1;
my $State_FoundQuery      = 2;
my $State_FoundAtLeastOne = 3;
my $State_FoundMatch      = 4;

my $state;           #state of the finite automaton
my $candidate_score; #test return value for whether a line has a score on it
my $final_score;     #score of the match
my $match_adjective; #Strong, Moderate, or Weak
my $num_queries;     #number of queries for which output has been processed
my $nextline;        #one line of input file
my $query;           #the identifier for the query
my $matchType;            #type of one match
my $suspectOutputString; #string for suspect origins
my $needToPrintMatch;    #0-1 valued variable that determines whether we print a match or not
my $startPos;            #starting position of a match to vector
my $endPos;              #ending position of a match to vector
my $vectorId;                    # string identifying a vector that matches a query
my $queryLength;                 # length of current query
my $testLength;                  # test value that will be either 0 or the length of the current query
my $internalMatch;               # is this match an internal match
my $seenVectorName = 0;          # have we seen the vector identifier for this match or not
my $subject_start;               # first aligned position in one block of the subject (database) sequence
my $subject_end;                 # last aligned position in one block of the subject (database) sequence
my $query_start;                 # first aligned position in one block of the query sequence
my $query_end;                   # last aligned position in one block of the query  sequence
my $i;                           # loop index
my $numMatches;                  # number of matches for this query
my $processing_alignment = 0;    # are we possibly in the middle of processing an alignment
my $overall_subject_start;       # subject start for an alignment that can have multiple blocks
my $overall_subject_end;         # subject end for an alignment that can have multiple blocks
my $overall_query_start;         # query start for an alignment that can have multiple blocks
my $overall_query_end;           # query end for an alignment that can have multiple blocks
my $has_strong_match;            # does the query have a strong match
my $has_moderate_match;          # does the query have a moderate match
my $has_weak_match;              # does the query have a weak match
my $overall_adjective;           # what is the strongest adjective for which this query has a match of that adjective
my $one_match_adjective;         # what is the adjective for a single alignment
my $alignments_one_query;        # how many alignments have we seen for one query
my $alignments_one_matched_pair; # how many alignments have we seen for one query-subject pair
my $identifier_filehandle;       # filehandle for file of core identifiers
my $terminal_filehandle;         # filehandle for terminal matches
my $internal_filehandle;         # filehandle for terminal matches

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
    'keep'               => \$GetOptions_H{"--keep"});


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
$output_internal_file = opt_Get("--output_internal", \%opt_HH);
$output_terminal_file = opt_Get("--output_terminal", \%opt_HH);
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
open($internal_FH, ">", $output_file_internal) or die "Cannot open $output_file_internal for writing\n"; 
open($terminal_FH, ">", $output_file_terminal) or die "Cannot open $output_file_terminal for writing\n"; 

# Parse input file (the vecscreen output):
#
# The following line types need to be recognized and distinguished:
#   A. BLASTN lines:            lines that begin with "BLASTN", indicating the start of output for a new query 
#                               isBLASTNLine() returns '1' 
#   B. Length lines:            lines that begin with "Length", including the length of the query or subject (vector)
#                               getLength() returns length (!0)
#   C. Query name lines:        lines that begin with "Query=", including the name of a query
#                               isQueryNameLine() returns '1'
#   D. Subject name lines:      lines that begin with "> gnl"   including the name of a subject (vector)
#                               isSubjectNameLine() returns '1'
#   E. Query alignment lines:   lines that begin with "Query "  including a segment of the query in an alignment
#                               isQueryAlignmentLine() returns '1'
#   F. Subject alignment lines: lines that begin with "Sbjct "  including a segment of the subject (vector) in an alignment
#                               isSubjectAlignmentLine() returns '1'
#   G. Match lines:             lines that begin with "Strong", "Moderate", "Weak" or "Suspect" indicating a type of match
#                               getMatch() returns name of match (!0)
#   H. Score lines:             lines that contain "Score" including score for a match
#                               getScore() returns score (!0)
#   I. Position lines:          lines that determine whether a query has a terminal match or has an internal match
#                               isPositionLine() returns '1'
#                               getPositions() returns $start and $stop
#line # line # relevant
# num # type # subroutine             # $state                 # <INPUT> line 
#-----#-------------------------------#------------------------#-----------------------------------------------------------------------
#   1 #    1 # isBLASTN               # ?                      #BLASTN 2.6.0+
#   2 #                               # $State_FoundBLASTN     #
#   3 #                               # "                      #
#   4 #                               # "                      # Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.
#   5 #                               # "                      # Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.
#   6 #                               # "                      # Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of
#   7 #                               # "                      # protein database search programs", Nucleic Acids Res. 25:3389-3402.
#   8 #                               # "                      #
#   9 #                               # "                      #
#  10 #                               # "                      #
#  11 #                               # "                      #Database: UniVec (build 9.0)
#  12 #                               # "                      #           5,456 sequences; 1,049,913 total letters
#  13 #                               # "                      #
#  14 # getMatch() => !0 ("Moderate") # "                      #Moderate match
#  15 # isPositionLine() =>1          # $State_FoundMatch      #1285	1333
#  16 # isPositionLine() =>1          # $State_FoundAtLeastOne #1408	1436
#  17 # getMatch() => !0 ("Suspect")  # "                      #Suspect origin
#  18 # isPositionLine() =>1          # "                      #1437	1467
#  19 #                               # "                      #
#  20 #                               # "                      #
#  21 # Type (3) isQueryNameLine()=>1 # "                      #Query= XM_715279.1 Candida albicans SC5314 Spl1p (SPL1), partial mRNA
#  22 #                               # $State_FoundQuery      #
#  23 # Line type (2)                 # "                      #Length=1467
#  24 #                               # "                      #
#  25 #                               # "                      #
#  26 # Type ? IsSubjectNameLine()=>1 # "                      #> gnl|uv|U29899.1:1847-4463 Cloning vector pACT2 MatchmakerII
#  27 # ignored b/c $seenVectorName=1 # "                      #Length=2617
#  28 #                               # "                      #
#  29 # getScore()=>!0             # "                      # Score = 58.6 bits (29),  Expect = 4e-06
#  30 #                               # "                      # Identities = 29/29 (100%), Gaps = 0/29 (0%)
#  31 #                               # "                      # Strand=Plus/Minus
#  32 #                               # "                      #
#  33 # isQueryAlignmentLine()=>1     # "                      # Query  1408  TTATGGGAAATGGTTCAAGAAGGTATTGA  1436
#  34 #                               # "                      #              |||||||||||||||||||||||||||||
#  35 # isSubjectAlignmentLine()=>1   # "                      # Sbjct  2140  TTATGGGAAATGGTTCAAGAAGGTATTGA  2112
#  36 #                               # "                      #
#  37 #                               # "                      #
#  38 # IsSubjectNameLine()=>1        # "                      #> gnl|uv|U03498.1:7366-8059 Yeast episomal vector YEp13
#  39 #                               #                        #Length=694
#  40 #                               #                        #
#  41 # getScore()=>!0             #                        # Score = 50.6 bits (25),  Expect = 0.001
#  42 #                               #                        # Identities = 45/49 (92%), Gaps = 0/49 (0%)
#  43 #                               #                        # Strand=Plus/Minus
#  44 #                               #                        # 
#  45 #                               #                        #Query  1285  GATGATGCCTTAGCTCATTCTTCAATTAGATTTGGTATTGGTAGATTTA  1333
#  46 #                               #                        #             |||||||| ||||| |||||||| || ||||||||||||||||||||||
#  47 #                               #                        #Sbjct  113   GATGATGCATTAGCCCATTCTTCCATCAGATTTGGTATTGGTAGATTTA  65
#  48 #                               #                        #
#  49 #                               #                        #
#  50 #                               #                        #
#  51 #                               #                        #  Lambda      K        H
#  52 #                               #                        #    1.39    0.747     1.38 
#  53 #                               #                        # 
#  54 #                               #                        # Gapped
#  55 #                               #                        # Lambda      K        H
#  56 #                               #                        #   1.39    0.747     1.38 
#  57 #                               #                        #
#  58 #                               #                        # Effective search space used: 1750000000000
#  59 #                               #                        #
#  60 #                               #                        #
#  61 #                               #                        # Database: UniVec (build 9.0)
#  62 #                               #                        #  Posted date:  Mar 23, 2015  3:44 PM
#  63 #                               #                        #  Number of letters in database: 1,049,913
#  64 #                               #                        #  Number of sequences in database:  5,456
#  66 #                               #                        #
#  67 #                               #                        #
#  68 #                               #                        #
#  69 #                               #                        # Matrix: blastn matrix 1 -5
#  70 #                               #                        # Gap Penalties: Existence: 3, Extension: 3
############################################################################################

# Initializations
$state       = $State_Naive;
$num_queries = 0;

while(defined($nextline = <INPUT>)) { 
  chomp($nextline);

  # determine type of line

  ######################################################################################
  $testLength = getLength($nextline);
  if ($testLength > 0 && (!$seenVectorName)) {
    # B. Length line: line that begins with "Length", including the length of the query or subject (vector)
    # *** line num 23 in above example file ***
    # the (!$seenVectorName) part ensures this is the query length and not vector(subject) length
    $queryLength = $testLength; # just stores $testLength (which may change in future loop iterations) in $queryLength
  }
  ######################################################################################

  ######################################################################################
  if (isQueryAlignmentLine($nextline)) {
    # E. Query alignment line:line that begins with "Query "including a segment of the query in an alignment
    ($query_start, $query_end) = getQueryAlignmentPositions($nextline); 
    if ($processing_alignment) { 
      # not the first Query Alignment line we've seen for this match, 
      # update $overall_query_end, but not $overall_query_start
      if (1 == abs($query_start - $overall_query_end)) { 
        $overall_query_end = $query_end;
      }
    }
    else {
      # the first Query Alignment line we've seen for this match
      $overall_query_start = $query_start;
      $overall_query_end   = $query_end;
    }
  } 
  ######################################################################################

  ######################################################################################
  if (isSubjectAlignmentLine($nextline)) {
    # F. Subject alignment line: line that begins with "Sbjct "  including a segment of the subject (vector) in an alignment
    ($subject_start, $subject_end) = getSubjectAlignmentPositions($nextline); 
    if ($processing_alignment) {
      # not the first Subject Alignment line we've seen for this match, 
      # update $overall_query_end, but not $overall_query_start
      if (1 == abs($subject_start - $overall_subject_end)) {
        $overall_subject_end = $subject_end;
      }
    }
    else {
      # the first Subject Alignment line we've seen for this match
      $overall_subject_start = $subject_start;
      $overall_subject_end   = $subject_end;
      $processing_alignment  = 1;
      $alignments_one_query++;
      $alignments_one_matched_pair++;
    }
  } 
  ######################################################################################

  ######################################################################################
  if (isBLASTNLine($nextline)) {
    # A. BLASTN line: line that begins with "BLASTN", indicating the start of output 
    #                 for a new query 
    # reset some state variables accordingly
    $seenVector  = 0;	
    $processing_alignment = 0;
    $state = $State_FoundBLASTN;

    if ($needToPrintMatch) {
      # printing situation 1: Print the final hit for a query to any
      #                       subject, all other hits will have been
      #                       printed in situations 2 or 3, as we see
      #                       either the next hit for the same
      #                       query/subject pair (situation 2) or the
      #                       final hit for a query/subject pair *for
      #                       which that same query has another hit to
      #                       a different subject (situation 3)*
      $internalMatch       = (isInternalMatch($overall_query_start, $overall_query_end, $queryLength)) ? 1 : 0;
      $overall_adjective   = getStrongestAdjective($has_strong_match, $has_moderate_match, $has_weak_match); 
      $one_match_adjective = getMatchAdjective($final_score,$internalMatch);
      selectivePrinting(1);
    }
    $needToPrintMatch = 0; # reset flag that we need to print a match
    if ($debug_mode && defined($query) && ($query =~m/$query_debug/)) {
      print "set needToPrintMatch to $needToPrintMatch\n";
    }
    # reset variables for next hit
    $numMatches                  = 0;
    $suspectOutputString         = "None";
    $internalMatch               = 0;
    $alignments_one_query        = 0;
    $alignments_one_matched_pair = 0;
    $has_strong_match = $has_moderate_match = $has_weak_match = 0;
    $num_queries++;
  }
  ######################################################################################

  ######################################################################################
  # Currently processing a sequence that has at least one vector hit in it
  if (($State_FoundMatch == $state) || ($State_FoundAtLeastOne == $state)) {
    if (isPositionLine($nextline)) {
      # *** line num 15 in above example file *** 
      ($startPos, $endPos) = getPositions($nextline);
      $state = $State_FoundAtLeastOne;
      if ($debug_mode && defined($query) && ($query =~m/$query_debug/)) {
        print "setting needToPrintMatch to 1 at location A; triggering line is $nextline\n"; 
      }
      $needToPrintMatch = 1; # we have a match, so we need to print it eventually
      if ("Suspect" eq $matchType) { 
        # line num 18 in above example file
        # update suspectOutputString which holds suspect regions
        if ("None" eq $suspectOutputString) {
          $suspectOutputString = $matchType . '[' . $startPos . ',' . $endPos . ']' . ';' ;
        }
        else {
          $suspectOutputString = $suspectOutputString . $matchType . '[' . $startPos . ',' . $endPos . ']' . ';' ;
        }
      }
    }
    else { 
      # not a position line (isPositionLine($nextline)) returned 0
      if ($matchType = getMatch($nextline)) {
        ($has_strong_match, $has_moderate_match, $has_weak_match) = updateMatchAdjectives($has_strong_match, $has_moderate_match, $has_weak_match, $matchType);
        $state = $State_FoundAtLeastOne;
        if ($debug_mode && defined($query) && ($query =~m/$query_debug/)) {
          print "setting needToPrintMatch to 1 at location B; triggering line is $nextline\n"; 
        }
        $needToPrintMatch = 1;
      }
      # print "$matchType\n";
    }
  }
  if (($State_FoundBLASTN == $state)) {
    if ($matchType = getMatch($nextline)) {
      # line num 14 in example file above
      ($has_strong_match, $has_moderate_match, $has_weak_match) = updateMatchAdjectives($has_strong_match, $has_moderate_match, $has_weak_match, $matchType);
      $state = $State_FoundMatch;
      # print "$matchType\n";
    }
    $numMatches++;
  }
  if (($State_FoundAtLeastOne == $state)) {
    if (isQueryNameLine($nextline)) {
      # line num 21 in example file above
      $query = getQuery($nextline);
      # print "query line: $query\n";
      $state = $State_FoundQuery;
    }
  }
  if ($State_FoundQuery == $state) {
    $candidate_score = getScore($nextline);
    if ($candidate_score > 0) {
      if ($alignments_one_matched_pair > 0) {
        # printing situation 2: X>1 alignments/hits between same query/subject pair
        #                       here we are printing alignment n where n < X.
        if (isInternalMatch($overall_query_start, $overall_query_end, $queryLength)) {
          $internalMatch = 1;
        }
        else {
          $internalMatch = 0;
        }
        $overall_adjective   = getStrongestAdjective($has_strong_match, $has_moderate_match, $has_weak_match); 
        $one_match_adjective = getMatchAdjective($final_score,$internalMatch);
        selectivePrinting(2);
        $processing_alignment = 0;
      }
      $final_score = $candidate_score;
    }
  }
  if (($State_FoundQuery == $state) && ($needToPrintMatch)) {
    if (isSubjectNameLine($nextline)) {
      # same query, but a new vector
      # *** line num 38 in example file above ***
      $seenVectorName = 1;
      if ($alignments_one_query > 0) {
        # if we've seen an alignment of this query to a previous subject (vector)
        # then output info on that previous alignment
        # printing situation 3: X>=1 alignments/hits between same query/subject pair
        #                       here we are printing alignment n == X (final alignment for that query/subject pair)
        $internalMatch       = isInternalMatch($overall_query_start, $overall_query_end, $queryLength);
        $overall_adjective   = getStrongestAdjective($has_strong_match, $has_moderate_match, $has_weak_match); 
        $one_match_adjective = getMatchAdjective($final_score,$internalMatch);
        selectivePrinting(3);
        # reset variables
        $alignments_one_matched_pair = 0;
        $alignments_one_query = 0;
      }
      # store new vector/subject name
      $vectorId = getVector($nextline);
      # print "$vectorId\n";
      $processing_alignment = 0;
    }
  }
}
#to handle last match

if ($needToPrintMatch) {
    if (isInternalMatch($overall_query_start, $overall_query_end, $queryLength)) {
	$internalMatch = 1;
    }
    else {
	$internalMatch = 0;
    }
    $overall_adjective = getStrongestAdjective($has_strong_match, $has_moderate_match, $has_weak_match); 
    $one_match_adjective = getMatchAdjective($final_score,$internalMatch);
    selectivePrinting(4);
}

# Removed this output line so output of from_vecscreen_to_summary_EPNEDITED.pl 
# is cleaner. Uncomment it if you want it back.
# print "Processed $num_queries queries\n";

close(INPUT);
close($internal_filehandle);
close($terminal_filehandle);

# Subroutine: isBLASTNLine() 
# Synopsis: Determines whether a line starts with the string BLASTN or not  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: 1 if yes, 0 if no
sub isBLASTNLine {
    my $sub_name = "isBLASTNLine()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;

    if (($local_one_line =~m/^BLASTN/)) {
	return 1;
    }
    else {
	return 0;
    }
}

# Subroutine: isPositionLine() 
# Synopsis: Determines whether a line starts with a digit
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: 1 if yes, 0 if no
sub isPositionLine {
    my $sub_name = "isPositionLine()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;

    if (($local_one_line =~m/^\d/)) {
	return 1;
    }
    else {
	return 0;
    }
}

# Subroutine: isQueryNameLine() 
# Synopsis: Determines whether a line starts with the string Query= or not  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: 1 if yes, 0 if no
sub isQueryNameLine {
    my $sub_name = "isQueryNameLine()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;

    if (($local_one_line =~m/^Query=/)) {
	return 1;
    }
    else {
	return 0;
    }
}

# Subroutine: getQuery() 
# Synopsis: Extracts the query identifier from a line that starts  Query=  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: the query identifier
sub getQuery {
    my $sub_name = "getQuery()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;

    my $local_query_id;

    ($local_query_id) = ($local_one_line =~m/^Query=\s+(\S+)/);
    return $local_query_id;
}

# Subroutine: isSubjectAlignmentLine() 
# Synopsis: Determines whether a line starts with the string Sbjct or not  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: 1 if yes, 0 if no
sub isSubjectAlignmentLine {
    my $sub_name = "isSubjectAlignmentLine()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;

    if (($local_one_line =~m/^Sbjct/)) {
	return 1;
    }
    else {
	return 0;
    }
}

# Subroutine: isQueryAlignmentLine() 
# Synopsis: Determines whether a line starts with the string Query or not  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: 1 if yes, 0 if no
sub isQueryAlignmentLine {
    my $sub_name = "isQueryAlignmentLine()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;

    if (($local_one_line =~m/^Query /)) {
	return 1;
    }
    else {
	return 0;
    }
}

# Subroutine: isSubjectNameLine() 
# Synopsis: Determines whether a line starts with the string > gnl| or not  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: 1 if yes, 0 if no
sub isSubjectNameLine {
    my $sub_name = "isSubjectNameLine()";
    my $nargs_exp = 1;


    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;

    if (($local_one_line =~m/^>\s*gnl/)) {
	return 1;
    }
    else {
	return 0;
    }
}

# Subroutine: getScore() 
# Synopsis: Determines whether a line includes the word Score  or not  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: score if yes, 0 if no
sub getScore {
    my $sub_name = "getScore()";
    my $nargs_exp = 1;


    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;
    my $local_score_return; #the score to return

    if (($local_one_line =~m/Score/)) {
	($local_score_return) = ($local_one_line =~m/\((\d+)\)/); 
	return $local_score_return;
    }
    else {
	return 0;
    }
}

# Subroutine: getVector() 
# Synopsis: Extracts the vector identifier from a line that starts  > gnl  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: the query identifier
sub getVector {
    my $sub_name = "getVector()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;

    my $local_vector_id;

    ($local_vector_id) = ($local_one_line =~m/^>\s*gnl\|(\S+)/);
    return $local_vector_id;
}

# Subroutine: getQuery() 
# Synopsis: Extracts the type of match from a line that may have a match  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: the query identifier
sub getMatch {
    my $sub_name = "getMatch()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;

    #print "Called getMatch $local_one_line\n";
    if ($local_one_line =~m/^Strong/) {
	return "Strong";
    }
    if ($local_one_line =~m/^Moderate/) {
	return "Moderate";
    }
    if ($local_one_line =~m/^Weak/) {
	return "Weak";
    }
    if ($local_one_line =~m/^Suspect/) {
	return "Suspect";
    }
    return 0;
}

# Subroutine: getPositions() 
# Synopsis: Extracts two numerical positions for a line that is supposed to contain them  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: the start position and the end position
sub getPositions {
    my $sub_name = "getPositions()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;

    my $local_start_pos; #stores the start position
    my $local_end_pos; #stores the end position
    #print "called getPositions $local_one_line\n";

    ($local_start_pos, $local_end_pos) = ($local_one_line =~m/(\d+)\s+(\d+)/);

    # print "from getPositions returning $local_start_pos\t$local_end_pos\n";
    return($local_start_pos, $local_end_pos);
}

# Subroutine: getSubjectAlignmentPositions() 
# Synopsis: Extracts two numerical positions for a line that is supposed to contain them  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: the start position and the end position
sub getSubjectAlignmentPositions {
    my $sub_name = "getSubjectAlignmentPositions()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;

    my $local_start_pos; #stores the start position
    my $local_end_pos; #stores the end position
    #print "called getSubjectAlignmentPositions $local_one_line\n";

    ($local_start_pos, $local_end_pos) = ($local_one_line =~m/Sbjct\s+(\d+)\s+\D+(\d+)/);

    # print "from getSubjectAlignmentPositions returning $local_start_pos\t$local_end_pos\n";
    return($local_start_pos, $local_end_pos);
}

# Subroutine: getQueryAlignmentPositions() 
# Synopsis: Extracts two numerical positions for a line that is supposed to contain them  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: the start position and the end position
sub getQueryAlignmentPositions {
    my $sub_name = "getQueryAlignmentPositions()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;

    my $local_start_pos; #stores the start position
    my $local_end_pos; #stores the end position
    #print "called getQueryAlignmentPositions $local_one_line\n";

    ($local_start_pos, $local_end_pos) = ($local_one_line =~m/Query\s+(\d+)\s+\D+(\d+)/);

    # print "from getQueryAlignmentPositions returning $local_start_pos\t$local_end_pos\n";
    return($local_start_pos, $local_end_pos);
}

# Subroutine: getLength() 
# Synopsis: If a line starts with Length=, returns the length; otherwise returns 0
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: the length or 0
sub getLength {
    my $sub_name = "getLength()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_one_line) = @_;

    my $local_length; #stores the length
    if ($local_one_line =~m/^Length/) {
	($local_length) = ($local_one_line =~m/Length=(\d+)/);
	return($local_length);
    }
    else {
	return 0;
    }
}

# Subroutine: getStrongestAdjective() 
# Synopsis: return the strongest match adjective out of Strong, Moderate, Weak
#
# Args: $local_has_strong, $local_has_moderate, $local_has_weak: whether each strength of match is present
#             
#             
#
# Returns: the adjective describing the strongest vector match for this sequence
sub getStrongestAdjective {
    my $sub_name = "getStrengthAdjective()";
    my $nargs_exp = 3;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_has_strong, $local_has_moderate, $local_has_weak) = @_;

    if ($local_has_strong) {
	return "Strong";
    }
    else {
	if ($local_has_moderate) {
	    return "Moderate";
	}
	else {
	    if ($local_has_weak) {
		return "Weak";
	    }
	    else {
		return "None";
	    }
	}
    }
}

# Subroutine: getMatchAdjective() 
# Synopsis: return the adjective for one match out of Strong, Moderate, Weak
#
# Args: $local_score, $local_is_internal
#             
#             
#
# Returns: the adjective (Strong, Moderate, or Weak) describing the quality of this match
sub getMatchAdjective {
    my $sub_name = "getMatchAdjective()";
    my $nargs_exp = 2;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_score, $local_is_internal) = @_;

    if ($local_is_internal) {
	if ($local_score >= $STRONG_SCORE_LOWER_INTERNAL) {
	    return "Strong";
	}
	else {
	    if ($local_score >= $MODERATE_SCORE_LOWER_INTERNAL) {
		return "Moderate";
	    }
	    else {
		if ($local_score >= $WEAK_SCORE_LOWER_INTERNAL) {
		    return "Weak";
		}
		else {
		    return "None"
		}
	    }
	}
    }
    else {
	if ($local_score >= $STRONG_SCORE_LOWER_TERMINAL) {
	    return "Strong";
	}
	else {
	    if ($local_score >= $MODERATE_SCORE_LOWER_TERMINAL) {
		return "Moderate";
	    }
	    else {
		if ($local_score >= $WEAK_SCORE_LOWER_TERMINAL) {
		    return "Weak";
		}
		else {
		    return "None"
		}
	    }
	}
    }
}

# Subroutine: updateMatchAdjectives() 
# Synopsis: update the local variables to indicate which severities of matches one query has
#
# Args: $input_has_strong_mtch, $input_has_moderate_match, $input_has_weak_match, $match_type
#             
#             
#
# Returns: updated values of $has_strong_match, $has_moderate_match, $has_weak_match
sub updateMatchAdjectives {
    my $sub_name = "updateMatchAdjectives()";
    my $nargs_exp = 4;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_has_strong_match, $local_has_moderate_match, $local_has_weak_match, $local_match_type) = @_;
    if ("Strong" eq $matchType) {
	$local_has_strong_match = 1;
    }
    if ("Moderate" eq $matchType) {
	$local_has_moderate_match = 1;
    }
    if ("Weak" eq $matchType) {
	$local_has_weak_match = 1;
    }
    return($local_has_strong_match, $local_has_moderate_match, $local_has_weak_match);
}

# Subroutine: selectivePrinting() 
# Synopsis: decides where to rpint a line of output
#
# Args: calling location
#             
#             
#
# Returns: nothing
sub selectivePrinting {
    my $sub_name = "selectivePrinting()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_calling_location) = @_;
    if ($debug_mode && $query =~m/$query_debug/) {
	print "printing query $query at location $local_calling_location; internal match is $internalMatch $overall_query_start $overall_query_end\n";
	print "needToPrintMatch is $needToPrintMatch; triggering line is: $nextline\n";
    }
    if ($internalMatch) {
	printOutputLine($internal_filehandle);
    }
    else {
	printOutputLine($terminal_filehandle);
    }
}

# Subroutine: printOutputLine() 
# Synopsis: return the adjective for one match out of Strong, Moderate, Weak
#
# Args: $local_filehandle
#             
#             
#
# Returns: nothing
sub printOutputLine {
    my $sub_name = "printOutputLine()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_filehandle) = @_;
    if ($verbose_mode) {
	print $local_filehandle "$query\t$overall_query_start\t$overall_query_end\t$vectorId\t$overall_subject_start\t$overall_subject_end\t$one_match_adjective\t$overall_adjective\t$suspectOutputString\n";
    }
    else {
	print $local_filehandle "$query\t$vectorId\t$overall_subject_start\t$overall_subject_end\n";
    }
}

# Subroutine: isInternalMatch() 
# Synopsis: i the match between vector sequence and quey deemed to be in the internal part of the sequence
#
# Args: $query_start, $query_end, $query_length
#             
#             
#
# Returns: 1 of the match is internal and 0 otherwise
sub isInternalMatch {
    my $sub_name = "isInternalMatch()";
    my $nargs_exp = 3;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
    my ($local_query_start, $local_query_end, $local_query_length) = @_;
    if (($local_query_start > $DISTANCE_FOR_TERMINAL) && (($local_query_length -$local_query_end +1) > $DISTANCE_FOR_TERMINAL)) {
	return 1;
    }
    else {
	return 0;
    }
}
