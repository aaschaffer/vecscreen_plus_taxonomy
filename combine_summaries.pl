#!/usr/bin/perl -w
# the first line of perl code has to be above
#
# Author: Alejandro Schaffer
# Code to combine summaries of matches from terminal and internal specific summaries; if strenght available, then Strong goes bove Moderate goes above Weak
# Usage: combine_summaries.pl --input_terminal <input terminal file> --input_internal <input internal file> (optional) --verbose  --outfile <combined output summary file>

use strict;
use warnings;
use Getopt::Long;

require "epn-options.pm";

my $input_terminal_file; #input fasta file 
my $input_internal_file; #input fasta file 
my $output_file; #output fasta file for sequences that do not have the forbidden terms
my $nextline; #one line of the file
my @terminal_lines; #array of forbidden terms
my @internal_lines; #number of forbidden terms
my @terminal_accessions; #accessions for each terminal line
my @internal_accessions; #accessions for each internal line
my @terminal_strengths; #strengths for each terminal line
my @internal_strengths; #strengths for each internal line
my @terminal_overall_strengths; #strengths for each match
my @internal_overall_strengths; #strengths for each match
my @terminal_printed; #has each terminal line been printed
my @internal_printed; #has each internal line been printed
my $num_terminal_lines;
my $num_internal_lines;
my $line_index; #counter for lines;
my $t; #loop index
my %internal_accession_hash; #first lines of each accession in the internal set
my %overall_strength_hash; #overall strength of an accession
my @fields; #entries in one row;
my $accession; #accession part of an identifier  
my $prev_accession; #previously seen accession
my $ACC_COLUMN =0; #the column with the identifier
my $STR_COLUMN = 6; #the column with the strength of the match
my $OVERALL_STR_COLUMN = 7; #the column with the overal highest strength for this match
my $next_terminal; #index of next terminal line to be considered for printing
my $next_internal; #index of next internal line to be considered for printing
my $special_next_internal; #index of next internal line to be considered for printing when accession matches terminal
my $num_printed; #number of lines printed so far
my $repeated_accession; #an accession that occurs in both terminal and interminal files
my @strength_seeking; #what strength of match are we seeking in this round
my $s; #loop index over strengths  
my $verbose_mode; #did user specify verbose mode, in which case match strengths are known and extra columns are printed

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
opt_Add("--input_terminal",           "string", undef,                        1,    undef, undef,     "input terminal summary file",  "File name <s> of terminal summary",     \%opt_HH, \@opt_order_A);
opt_Add("--input_internal",           "string", undef,                        1,    undef, undef,     "input internal summary file",  "File name <s> of internal summary",     \%opt_HH, \@opt_order_A);
opt_Add("--verbose",           "boolean", 0,                        1,    undef, undef,      "be verbose",                                   "be verbose in output",  \%opt_HH, \@opt_order_A);
opt_Add("--outfile",           "string", undef,                        1,    undef, undef,     "output file",  "Name <s> of output file",     \%opt_HH, \@opt_order_A);


# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "combine_summaries.pl: combine summaries of terminal and internal matches\n";
$usage   .= "Usage: --input_terminal <input terminal file> --input_internal <input intenral  file> (optional) --verbose --outfile <output file>\n";
$usage   .= "\nFor example:\ncombine_summaries.pl --input_terminal summary_terminal.txt --input_internal summary_internal.txt (optional --verbose  --outfile combined_summary.txt\n";

my $options_okay =
    &GetOptions('h'            => \$GetOptions_H{"-h"},
# basic options
                'input_terminal=s'            => \$GetOptions_H{"--input_terminal"},
                'input_internal=s'            => \$GetOptions_H{"--input_internal"},
                'verbose'            => \$GetOptions_H{"--verbose"},
                'outfile=s'            => \$GetOptions_H{"--outfile"});


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

$input_terminal_file = opt_Get("--input_terminal", \%opt_HH);
$input_internal_file = opt_Get("--input_internal", \%opt_HH);
$verbose_mode = opt_Get("--verbose", \%opt_HH);
$output_file = opt_Get("--outfile", \%opt_HH);

open(TERMINAL, "<$input_terminal_file") or die "Cannot open 1 $input_terminal_file\n"; 
open(INTERNAL, "<$input_internal_file") or die "Cannot open 2 $input_internal_file\n"; 
open(OUTPUT, ">$output_file") or die "Cannot open 3 $output_file\n"; 

$num_terminal_lines = 0;
while(defined($nextline = <TERMINAL>)) {
    chomp($nextline);
    $terminal_lines[$num_terminal_lines] = $nextline;
    @fields = split /\t/, $nextline;
    $accession = $fields[$ACC_COLUMN];
    $terminal_accessions[$num_terminal_lines] = $accession;
    if ($verbose_mode) {
	$terminal_strengths[$num_terminal_lines] = $fields[$STR_COLUMN];
    }
    else {
	$terminal_strengths[$num_terminal_lines] = "Strong";
    }
    $terminal_overall_strengths[$num_terminal_lines] = $fields[$OVERALL_STR_COLUMN];
    if (!defined($overall_strength_hash{$accession})) {
	if ($verbose_mode) {
	    $overall_strength_hash{$accession} = $fields[$OVERALL_STR_COLUMN];
	}
	else {
	    $overall_strength_hash{$accession} = "Strong";
	}
    }
    $num_terminal_lines++;
}
close(TERMINAL);

for ($line_index = 0; $line_index < $num_terminal_lines; $line_index++) {
    $terminal_printed[$line_index] = 0;
}

$num_internal_lines = 0;
$prev_accession = "NONE";
while(defined($nextline = <INTERNAL>)) {
    chomp($nextline);
    $internal_lines[$num_internal_lines] = $nextline;
    @fields = split /\t/, $nextline;
    $accession = $fields[$ACC_COLUMN];
    $internal_accessions[$num_internal_lines] = $accession;
    if ($verbose_mode) {
	$internal_strengths[$num_internal_lines] = $fields[$STR_COLUMN];
    }
    else {
	$internal_strengths[$num_internal_lines] = "Strong";
    }
    $internal_overall_strengths[$num_internal_lines] = $fields[$OVERALL_STR_COLUMN];
    if (!defined($overall_strength_hash{$accession})) {
	if ($verbose_mode) {
	    $overall_strength_hash{$accession} = $fields[$OVERALL_STR_COLUMN];
	}
	else {
	    $overall_strength_hash{$accession} = "Strong";
	}
    }
    else {
	if ($verbose_mode) {
	    if(compareStrengths($fields[$OVERALL_STR_COLUMN], $overall_strength_hash{$accession})) {
		$overall_strength_hash{$accession} = $fields[$OVERALL_STR_COLUMN];
	    }
	}
    }
    if (! ($prev_accession eq $accession)) {
	$internal_accession_hash{$accession} = $num_internal_lines; 
    } 
    $num_internal_lines++;
    $prev_accession = $accession;
}
close(INTERNAL);

for ($line_index = 0; $line_index < $num_internal_lines; $line_index++) {
    $internal_printed[$line_index] = 0;
}

$num_printed = 0;

$strength_seeking[0] = "Strong";
$strength_seeking[1] = "Moderate";
$strength_seeking[2] = "Weak";

for ($s = 0; $s < 3; $s++){
    $line_index = 0;
    $next_terminal = 0;
    $next_internal = 0;
    while (($next_terminal < $num_terminal_lines) && ($next_internal < $num_internal_lines)) {
	#print "terminal $next_terminal  $terminal_accessions[$next_terminal]\n";
	#print "strength is $overall_strength_hash{$terminal_accessions[$next_terminal]}\n";
	if ($next_terminal < $num_terminal_lines) {
	    if ($overall_strength_hash{$terminal_accessions[$next_terminal]} eq $strength_seeking[$s]) {
		if (defined($internal_accession_hash{$terminal_accessions[$next_terminal]})) {
		    $special_next_internal = $internal_accession_hash{$terminal_accessions[$next_terminal]};
		    $repeated_accession = $terminal_accessions[$next_terminal];
		    #print "Repeated accession is $repeated_accession  $next_terminal  $num_terminal_lines $special_next_internal $num_internal_lines\n";
		    while (($next_terminal < $num_terminal_lines) && ($terminal_accessions[$next_terminal] eq $repeated_accession) &&
			   ($special_next_internal < $num_internal_lines) && ($internal_accessions[$special_next_internal] eq $repeated_accession)){
			if (compareStrengths($terminal_strengths[$next_terminal], $internal_strengths[$special_next_internal])) {
			    print OUTPUT "$terminal_lines[$next_terminal]";
			    print OUTPUT "\n";
			    $terminal_printed[$next_terminal] = 1;
			    $num_printed++; 
			    $next_terminal++;
			} 
			else {
			    print OUTPUT "$internal_lines[$special_next_internal]";
			    print OUTPUT "\n";
			    $internal_printed[$special_next_internal] = 1;
			    $num_printed++; 
			    $special_next_internal++;
			}
		    }
		    while (($next_terminal < $num_terminal_lines) && ($terminal_accessions[$next_terminal] eq $repeated_accession)) {
			print OUTPUT "$terminal_lines[$next_terminal]";
			print OUTPUT "\n";
			$terminal_printed[$next_terminal] = 1;
			$num_printed++; 
			$next_terminal++;
		    }
		    while (($special_next_internal < $num_internal_lines) && ($internal_accessions[$special_next_internal] eq $repeated_accession)){
			print OUTPUT "$internal_lines[$special_next_internal]";
			print OUTPUT "\n";
			$internal_printed[$special_next_internal] = 1;
			$num_printed++; 
			$special_next_internal++;
		    }
		}
		else {
		    if ((! ($terminal_printed[$next_terminal])) && ($overall_strength_hash{$terminal_accessions[$next_terminal]} eq $strength_seeking[$s])) {
			print OUTPUT "$terminal_lines[$next_terminal]";
			print OUTPUT "\n";
			$terminal_printed[$next_terminal] = 1;
			$num_printed++; 
		    }
		    $next_terminal++;
		}
	    }
	    else {
		$next_terminal++;
	    }
	}
	while (($next_terminal == $num_terminal_lines) && ($next_internal < $num_internal_lines)) {
	    #print "Internal $next_internal\n";
	    #print "$internal_accessions[$next_internal]\n";
	    #print "Internal strength is $overall_strength_hash{$internal_accessions[$next_internal]}\n";
	    if ((! ($internal_printed[$next_internal])) && ($overall_strength_hash{$internal_accessions[$next_internal]} eq $strength_seeking[$s]) ) {
		print OUTPUT "$internal_lines[$next_internal]";
		print OUTPUT "\n";
		$internal_printed[$next_internal] = 1;
		$num_printed++; 
	    }
	    $next_internal++;
	}
    }
}
close (OUTPUT);


# Subroutine: compare two strengths()
# Synopsis: return 1 if first stregth is greater than or equal to secon strength
#
# Args: $local_strength1
#       $local_strength2
#
#
# Returns: 1 if the first strength is greater than or equal to the second strength
sub compareStrengths {
    my $sub_name = "compareStrengths()";
    my $nargs_exp = 2;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_strength1, $local_strength2) = @_;
    my $local_s1_number = 0;
    my $local_s2_number = 0; #convert strength words to numbers

    if ("Strong" eq $local_strength1) {
	$local_s1_number = 3;
    }
    if ("Moderate" eq $local_strength1) {
	$local_s1_number = 2;
    }
    if ("Weak" eq $local_strength1) {
	$local_s1_number = 1;
    }
    if ("Strong" eq $local_strength2) {
	$local_s2_number = 3;
    }
    if ("Moderate" eq $local_strength2) {
	$local_s2_number = 2;
    }
    if ("Weak" eq $local_strength2) {
	$local_s2_number = 1;
    }
    return($local_s1_number >= $local_s2_number);
}
