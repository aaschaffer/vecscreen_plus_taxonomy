#!/usr/bin/perl -w
# the first line of perl code has to be above
#
# Author: Alejandro Schaffer with help from Eric Nawrocki
# Code to parse vecscreen output for vector screening.
#
# Usage: parse_vecscreen.pl --input <input file> (optional) --verbose --outfile_internal <output internal> --outfile_terminal <output terminal>

use strict;
use warnings;
use Getopt::Long;

require "epn-options.pm";


my $DISTANCE_FOR_TERMINAL = 25; #number of nucleotides to the end for a match to be considered terminal
my $STRONG_SCORE_LOWER_TERMINAL = 24; #lower score threshold for a terminal match to be a strong match
my $STRONG_SCORE_LOWER_INTERNAL = 30; #lower score threshold for an internal match to be a strong match
my $MODERATE_SCORE_UPPER_TERMINAL = 23; #upper score threshold for a terminal match to be a moderate match
my $MODERATE_SCORE_UPPER_INTERNAL = 29; #upper score threshold for an internal match to be a moderate match
my $MODERATE_SCORE_LOWER_TERMINAL = 19; #lower score threshold for a terminal match to be a moderate match
my $MODERATE_SCORE_LOWER_INTERNAL = 25; #lower score threshold for an internal match to be a moderate match
my $WEAK_SCORE_LOWER_TERMINAL = 16; #lower score threshold for a terminal match to be a weak match
my $WEAK_SCORE_LOWER_INTERNAL = 23; #lower score threshold for an internal match to be a weak match

my $State_Naive = 0;
my $State_FoundBLASTN = 1;
my $State_FoundQuery = 2;
my $State_FoundAtLeastOne = 3;
my $State_FoundMatch = 4;

my $input_file; #input file of vecscreen output
my $output_file_internal; #output file for sequences from Genbank with internal matches
my $output_file_terminal; #output file for sequences from Genbank with only terminal matches
my $state; #state of the finite automaton
my $candidate_score; #test return value for whether a line has a score on it
my $final_score; #score of the match
my $match_adjective; #Strong, Moderate, or Weak
my $num_queries; #number of queries for which output has been processed
my $nextline; #one line of input file
my $query; #the identifier for the query
my $matchType; #type of one match
my $suspectOutputString; #string for suspect origins
my $needToPrintMatch; #0-1 valued variable that determines whether we print a match or not
my $startPos; #starting position of a match to vector
my $endPos; #ending position of a match to vector
my $vectorId; #string identifying a vector that matches a query
my $queryLength; #length of current query
my $testLength; #test value that will be either 0 or the length of the current query
my $internalMatch; #is this match an internal match
my $seenVector = 0; #have we seen the vector identifier for this match or not
my $subject_start; #first aligned position in one block of the subject (database) sequence
my $subject_end; #last aligned position in one block of the subject (database) sequence
my $query_start; #first aligned position in one block of the query sequence
my $query_end; #last aligned position in one block of the query  sequence
my $i; #loop index
my $numMatches; #number of matches for this query
my $processing_alignment = 0; #are we possibly in the middle of processing an alignment
my $overall_subject_start; #subject start for an alignment that can have multiple blocks
my $overall_subject_end; #subject end for an alignment that can have multiple blocks
my $overall_query_start; #query start for an alignment that can have multiple blocks
my $overall_query_end; #query end for an alignment that can have multiple blocks
my $has_strong_match; #does the query have a strong match
my $has_moderate_match; #does the query have a moderate match
my $has_weak_match; #does the query have a weak match
my $overall_adjective; #what is the strongest adjective for which this query has a match of that adjective
my $one_match_adjective; #what is the adjective for a single alignment
my $alignments_one_query; #how many alignments have we seen for one query
my $alignments_one_matched_pair; #how many alignments have we seen for one query-subject pair
my $identifier_filehandle; #filehandle for file of core identifiers
my $terminal_filehandle; #filehandle for terminal matches
my $internal_filehandle; #filehandle for terminal matches
my $test_query_debug = "NONE";
my $verbose_mode; #did user specify verbose mode, in which extra columns are printed

my %opt_HH = ();
my @opt_order_A = ();
my %opt_group_desc_H = ();

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
# The opt_Add() function is the way we add options to %opt_HH.
# It takes values of for each of the 2nd dim keys listed above.
$opt_group_desc_H{"1"} = "basic options";
#       option          type       default               group   requires incompat   preamble-outfile   help-outfile
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,     undef,            "display this help",                  \%opt_HH, \@opt_order_A);
opt_Add("--input",           "string", undef,                        1,    undef, undef,     "input fasta file",  "File name <s> with vecscreen output",     \%opt_HH, \@opt_order_A);
opt_Add("--verbose",           "boolean", 0,                        1,    undef, undef,      "be verbose",                                   "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--outfile_internal",           "string", undef,                        1,    undef, undef,     "output of sequences with internal matches",  "File name <s> to hold output sequences, with internal matches",     \%opt_HH,    \@opt_order_A);
opt_Add("--outfile_terminal",           "string", undef,                        1,    undef, undef,     "output of sequences with terminal matches",  "File name <s> to hold output sequences, with terminal matches",     \%opt_HH,    \@opt_order_A);


# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "parse_vecscreen.pl: convert vecscreen output file to one or more tab-delimited output files\n";
$usage    .= "parse_vecscreen.pl --input <input file> (optional) --verbose  --outfile_internal <output internal> --outfile_terminal <output terminal>\n";
$usage   .= "\nFor example:\nparse_vecscreen.pl --input vecscreen_output.txt  --outfile_internal output_internal.txt --outfile_terminal output_terminal.txt\n";

my $options_okay =
    &GetOptions('h'            => \$GetOptions_H{"-h"},
# basic options
                'input=s'            => \$GetOptions_H{"--input"},
                'verbose'            => \$GetOptions_H{"--verbose"},
                'outfile_internal=s'            => \$GetOptions_H{"--outfile_internal"},
                'outfile_terminal=s'            => \$GetOptions_H{"--outfile_terminal"});


# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) {
    opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
    if(! $options_okay) { die "ERROR, unrecognized option;"; }
    else                { exit 0; } # -h, exit with 0 status
}

# set options in %opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

$input_file = opt_Get("--input", \%opt_HH);
$verbose_mode = opt_Get("--verbose", \%opt_HH);
$output_file_internal = opt_Get("--outfile_internal", \%opt_HH);
$output_file_terminal = opt_Get("--outfile_terminal", \%opt_HH);
if (!defined($input_file)) {
    die("--input was not specified; exiting\n");
} 

if (!defined($output_file_internal)) {
    die("--outfile_internal was not specified; exiting\n");
} 
if (!defined($output_file_terminal)) {
    die("--outfile_terminal was not specified; exiting\n");
} 



open(INPUT, "<$input_file") or die "Cannot open 1 $input_file\n"; 
open($internal_filehandle, ">$output_file_internal") or die "Cannot open 5 $output_file_internal\n"; 
open($terminal_filehandle, ">$output_file_terminal") or die "Cannot open 6 $output_file_terminal\n"; 

$state = $State_Naive;
$num_queries = 0;

#The following line types need to be recognized and distinguished:
#lines that include the token BLASTN, indicating the start of output for a new query 
#lines that have the word Length indicating the length of the query
#lines that have the word Query on them to indicate a segment of the query in an alignment
#lines that have the word Subject on them to indicate a segment of the subject (here a vector sequence) in an alignment
#lines that distinguish whether a Query/Subject pair are the first part of an alignment or the continuation of an alignment
#positions that determine whether a query has a terminal match or has an internal match
while(defined($nextline = <INPUT>)) {
    chomp($nextline);
    $testLength = getLength($nextline);
    if ($testLength > 0 && (!$seenVector)) {
	$queryLength = $testLength;
    }
    if (isQueryAlignmentLine($nextline)) {
      ($query_start,$query_end) = getQueryAlignmentPositions($nextline); 
      if ($processing_alignment) {
        # EPN QUESTION: isn't it an error if the following line is not true?
        #               if we're processing an alignment shouldn't we always
        #               update $oveall_query_end here? If so, do you either 
        #               want to put a die() statement for the unlikely event
        #               that the following if is not true, OR remove the if()
        #               if you're certain that it is always true?
	  if (1 == abs($query_start - $overall_query_end)) { 
	      $overall_query_end = $query_end;
	  }
      }
      else {
	  $overall_query_start = $query_start;
	  $overall_query_end = $query_end;
      }
    } 
    if (isSubjectLine($nextline)) {
      ($subject_start,$subject_end) = getSubjectAlignmentPositions($nextline); 
      if ($processing_alignment) {
        # EPN QUESTION: isn't it an error if the following line is not true?
        #               if we're processing an alignment shouldn't we always
        #               update $oveall_query_end here? If so, do you either 
        #               want to put a die() statement for the unlikely event
        #               that the following if is not true, OR remove the if()
        #               if you're certain that it is always true?
	  if (1 == abs($subject_start - $overall_subject_end)) {
	      $overall_subject_end = $subject_end;
	  }
      }
      else {
	  $overall_subject_start = $subject_start;
	  $overall_subject_end = $subject_end;
	  $processing_alignment = 1;
	  $alignments_one_query++;
	  $alignments_one_matched_pair++;
      }
    } 
    #starting processing of next set of BLASTN results for the next query; resetting some state variables accordingly
    if (isBLASTNLine($nextline)) {
      $seenVector  = 0;	
      $processing_alignment = 0;
      $state = $State_FoundBLASTN;
      if ($needToPrintMatch) {
	  if (isInternalMatch($overall_query_start, $overall_query_end, $queryLength)) {
	      $internalMatch = 1;
	  }
	  else {
	      $internalMatch = 0;
	  }
	  $overall_adjective = getStrongestAdjective($has_strong_match, $has_moderate_match, $has_weak_match); 
	  $one_match_adjective = getMatchAdjective($final_score,$internalMatch);
	  selectivePrinting(1);
      }
      $needToPrintMatch = 0;
      if (defined($query) && ($query =~m/$test_query_debug/)) {
	  print "set needToPrintMatch to $needToPrintMatch\n";
      }
      $numMatches = 0;
      $suspectOutputString = "None";
      $internalMatch = 0;
      $has_strong_match = $has_moderate_match = $has_weak_match = 0;
      $alignments_one_query = 0;
      $alignments_one_matched_pair = 0;
      $num_queries++;
    }
    if (($State_FoundMatch == $state) || ($State_FoundAtLeastOne == $state)) {
	if (isPositionLine($nextline)) {
	    ($startPos, $endPos) = getPositions($nextline);
	    $state = $State_FoundAtLeastOne;
	    if (defined($query) && ($query =~m/$test_query_debug/)) {
		print "setting needToPrintMatch to 1 at location A; triggering line is $nextline\n"; 
	    }
	    $needToPrintMatch = 1;
	    if ("Suspect" eq $matchType) {
		if ("None" eq $suspectOutputString) {
		    $suspectOutputString = $matchType . '[' . $startPos . ',' . $endPos . ']' . ';' ;
		}
		else {
		    $suspectOutputString = $suspectOutputString . $matchType . '[' . $startPos . ',' . $endPos . ']' . ';' ;
		}
	    }
	}
	else {
	    if ($matchType = getMatch($nextline)) {
		($has_strong_match, $has_moderate_match, $has_weak_match) = updateMatchAdjectives($has_strong_match, $has_moderate_match, $has_weak_match, $matchType);
		$state = $State_FoundAtLeastOne;
		if (defined($query) && ($query =~m/$test_query_debug/)) {
		    print "setting needToPrintMatch to 1 at location B; triggering line is $nextline\n"; 
		}
		$needToPrintMatch = 1;
	    }
	    # print "$matchType\n";
	}
    }
    if (($State_FoundBLASTN == $state)) {
	if ($matchType = getMatch($nextline)) {
	    ($has_strong_match, $has_moderate_match, $has_weak_match) = updateMatchAdjectives($has_strong_match, $has_moderate_match, $has_weak_match, $matchType);
	    $state = $State_FoundMatch;
	    # print "$matchType\n";
	}
        # EPN QUESTION: looks like you're updating numMatches here even if getMatch() doesn't return a non-zero value
        # I think you want numMatches inside the above block, but you aren't using numMatches anywhere else so maybe
        # you just want to delete it.
	$numMatches++;
    }
    if (($State_FoundAtLeastOne == $state)) {
	if (isQueryLine($nextline)) {
	    $query = getQuery($nextline);
	    # print "query line: $query\n";
	    $state = $State_FoundQuery;
	}
    }
    if ($State_FoundQuery == $state) {
	$candidate_score = isScoreLine($nextline);
	if ($candidate_score > 0) {
	    if ($alignments_one_matched_pair > 0) {
		if (isInternalMatch($overall_query_start, $overall_query_end, $queryLength)) {
		    $internalMatch = 1;
		}
		else {
		    $internalMatch = 0;
		}
		$overall_adjective = getStrongestAdjective($has_strong_match, $has_moderate_match, $has_weak_match); 
		$one_match_adjective = getMatchAdjective($final_score,$internalMatch);
                selectivePrinting(2);
		$processing_alignment = 0;
	    }
	    $final_score = $candidate_score;
	}
    }
    if (($State_FoundQuery == $state) && ($needToPrintMatch)) {
	if (isVectorLine($nextline)) {
	    $seenVector = 1;
	    if ($alignments_one_query > 0) {
		if (isInternalMatch($overall_query_start, $overall_query_end, $queryLength)) {
		    $internalMatch = 1;
		}
		else {
		    $internalMatch = 0;
		}
		$overall_adjective = getStrongestAdjective($has_strong_match, $has_moderate_match, $has_weak_match); 
		$one_match_adjective = getMatchAdjective($final_score,$internalMatch);
                selectivePrinting(3);
		$alignments_one_matched_pair = 0;
		$alignments_one_query = 0;
	    }
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

# Subroutine: isQueryLine() 
# Synopsis: Determines whether a line starts with the string Query= or not  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: 1 if yes, 0 if no
sub isQueryLine {
    my $sub_name = "isQueryLine()";
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

# Subroutine: isSubjectLine() 
# Synopsis: Determines whether a line starts with the string Sbjct or not  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: 1 if yes, 0 if no
sub isSubjectLine {
    my $sub_name = "isSubjectLine()";
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

# Subroutine: isVectorLine() 
# Synopsis: Determines whether a line starts with the string > gnl| or not  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: 1 if yes, 0 if no
sub isVectorLine {
    my $sub_name = "isVectorLine()";
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

# Subroutine: isScoreLine() 
# Synopsis: Determines whether a line includes the word Score  or not  
#
# Args: $local_one_line: the line of data       
#             
#             
#
# Returns: score if yes, 0 if no
sub isScoreLine {
    my $sub_name = "isScoreLine()";
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
    if ($query =~m/$test_query_debug/) {
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
