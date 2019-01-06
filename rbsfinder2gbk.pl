#!/usr/bin/perl

#===============================================================================
#   Author: Alejandro Manzano Marin
#
#   File: rbsfinder2gbk.pl
#   Date: 22-06-2014
#   Version: 1.0
#
#   Usage:
#      perl rbsfinder2gbk.pl --infile|-i <tRNAscanSEout structure file> --genetab <tabular gene file> [options]
#
#      Check out 'perl rbsfinder2gbk.pl -h' for short usage manual and info on the software.
#
#    Description: This program will take as input the structure file output of tRNAscanSE and will print
#                 the annotation of tRNA genes in genbank format.
#                 
#
#    Contact: Contact the author at alejandro.manzano@uv.es using 'rbsfinder2gbk: ' as
#             as begining for subject for bug reporting, feature request or whatever
#             (of course realted to the software).
#
#    COPYRIGHT: Copyright (C) 2014  Alejandro Manzano-Marin.
#
#    LICENCE: This program is free software: you can redistribute it and/or modify it under the terms
#             of the GNU General Public License as published by the Free Software Foundation, either
#             version 3 of the License, or (at your option) any later version.
#             This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#             without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#             See the GNU General Public License for more details.
#             You should have received a copy of the GNU General Public License along with this program.
#             If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


# Load modules
use Getopt::Long;
use Pod::Usage;


# Define subroutines
sub printVersion {
	my $software= $_[0];
	my $version= $_[1];
	my $fileHandle= $_[2];
	print $fileHandle "$software v$version\n";
	exit (0);
}

sub max {
	my @array=@_;
	my $max=$array[0];
	for (my $i=1; $i<scalar(@array); $i++){
		if ($array[$i]>$max){
			$max=$array[$i];
		}
	}
	return ($max);
}

sub min {
	my @array=@_;
	my $min=$array[0];
	for (my $i=1; $i<scalar(@array); $i++){
		if ($array[$i]<$min){
			$min=$array[$i];
		}
	}
	return ($min);
}

sub switchVal {
	my $aRef= $_[0];
	my $bRef= $_[1];
	my $temp= 0;
	$temp= $$aRef;
	$$aRef= $$bRef;
	$$bRef= $temp;
	return (0);
}

sub revComp {
	my $seqRef= $_[0];
	my $newSeq= '';
	$newSeq= scalar(reverse($$seqRef));
	$newSeq=~ tr/atgcATGC/tacgTACG/;
	return ($newSeq);
}

sub printRBSfeat{
	my $fileHandle= $_[0];
	my $locus_tag= $_[1];
	my $RBS_start= $_[2];
	my $RBS_pattern= $_[3];
	my $strand= $_[4];
	my $gene_name= $_[5];
	my $RBS_motif= $_[6];
	
	my $RBSpos_string= '';
	
	if ($strand == -1){
		$RBSpos_string= 'complement(' . ($RBS_start-length($RBS_pattern)+1) . '..' . $RBS_start . ')';
	}
	else {
		$RBSpos_string= $RBS_start . '..' . ($RBS_start+length($RBS_pattern)-1);
	}
	print $fileHandle '     regulatory      ' . $RBSpos_string . "\n";
	print $fileHandle "                     /regulatory_class=\"ribosome_binding_site\"\n";
	print $fileHandle '                     /locus_tag="' . $locus_tag . "\"\n";
	if ($gene_name ne 'NULL'){
		print $fileHandle '                     /gene="' . $gene_name . "\"\n";
	}
	print $fileHandle '                     /inference="COORDINATES:nucleotide motif:RBSfinder:' . $RBS_motif . '"' . "\n";
#	print $fileHandle '                     /note="putative RBS"' . "\n";
}


# Variable definition

## Define other variables
my $file= '';
my $line= '';
my $id= '';
my $flag= 0;
my %reads= ();
my $inHead= '>';
my %aaNames=(
	'Ala' => ['trnA', 'Alanine'],
	'Arg' => ['trnR', 'Arginine'],
	'Asn' => ['trnN', 'Asparagine'],
	'Asp' => ['trnD', 'Aspartic acid'],
	'Cys' => ['trnC', 'Cysteine'],
	'Glu' => ['trnE', 'Glutamic acid'],
	'Gln' => ['trnQ', 'Glutamine'],
	'Gly' => ['trnG', 'Glycine'],
	'His' => ['trnH', 'Histidine'],
	'Ile' => ['trnI', 'Isoleucine'],
	'Leu' => ['trnL', 'Leucine'],
	'Lys' => ['trnK', 'Lysine'],
	'Met' => ['trnM', 'Methionine'],
	'Phe' => ['trnF', 'Phenylalanine'],
	'Pro' => ['trnP', 'Proline'],
	'Ser' => ['trnS', 'Serine'],
	'Thr' => ['trnT', 'Threonine'],
	'Trp' => ['trnW', 'Tryptophan'],
	'Tyr' => ['trnY', 'Tyrosine'],
	'Val' => ['trnV', 'Valine'],
	'Sec' => ['trnU', 'Selenocysteine'],
	'Pyl' => ['trnO', 'Pyrrolysine'],	
);

#my %geneticCode=('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11');

#my %{$geneticCode{'1'}}=(
#	'1' => ,
#	'2' => ,
#	'3' => ,
#	'4' => ,
#	'5' => ,
#	'6' => ,
#	'7' => ,
#	'8' => ,
#	'9' => ,
#	'10' => ,
#	'11' => ,
#	
#);

## General variables
my $PROGRAMNAME= 'rbsfinder2gbk';
my $VERSION= '1.0';

## Define options default values
my $opt_inFile= '';

my $opt_geneTab= '';

my $opt_verbose= 0;
my $opt_man= 0;
my $opt_help= 0;
my $opt_printVersion= 0;

## Define options hash
GetOptions(\%opts, 
	'infile|i=s' => \$opt_inFile, 
	'geneTab=s' => \$opt_geneTab, 
	'verbose|v!' => \$opt_verbose, 
	'help|h!' => \$opt_help, 
	'man!'  => \$opt_man, 
	'version!' => \$opt_printVersion) || pod2usage(-exitval => 1,  -verbose => 2);

if ($opt_help){
	pod2usage(-exitval => 1,  -verbose => 1);
}
if ($opt_man){
	pod2usage(-exitval => 0, -verbose => 2);
}
if ($opt_printVersion){
	&printVersion($PROGRAMNAME, $VERSION, \*STDERR);
}

# Script documetation

=pod

=head1 NAME

rbsfinder2gbk

=head1 VERSION

rbsfinder2gbk v1.0

=head1 SYNOPSIS

perl rbsfinder2gbk.pl --infile|-i <tRNAscanSEout structure file> --genetab <tabular gene file> [--verbose|-v] [--help|-h] [--man] [--version]

=head1 DESCRIPTION

This program will take as input the structure file output of tRNAscanSE and will print the annotation of tRNA genes in genbank format

=head1 OPTIONS

=head2 INPUT

=over 8

=item B<-i> | B<--infile> <string> (mandatory)

File containing sequences to be filtered in FAST(A/Q) format.

=item B<--genetab> <string> (mandatory)

Tabular gene file in the format '<locus_tag>\t<start>\t<end>\t<strand>\t<gene name>' as the ne used as input file for rbsfinder.pl, but including in the fourth column the gene name.

=back

=head2 INFO AND HELP

=over 8

=item B<-v> | B<--verbose> <boolean> (default: 0)

Prints status and info messages while processing.

=item B<-h> | B<--help> <boolean>

Print useful help on using this script.

=item B<--man> <boolean>

Print the full documentation.

=item B<--version> <boolean>

Print program version.

=back

=head1 AUTHOR

Alejandro Manzano-Marin, C<< <alejandro_dot_manzano_at_uv_dot_es> >>

=head1 BUGS

If you find a bug please email me at C<< <alejandro_dot_manzano_at_uv_dot_es> >> so that I can keep making rbsfinder2gbk better.

=head1 COPYRIGHT

Copyright (C) 2013  Alejandro Manzano-Marin.

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut


# Check options
if (!$opt_inFile){
	print STDERR "ERROR:\n";
	if (!$opt_inFile){
		print STDERR "rbsfinder output file missing\n";
	}
	print STDERR "Please check usage manual\n";
	pod2usage(-exitval => 1,  -verbose => 1);
}


# Assign values to variables dependant on options


# Main program

## Read list of accepted reads
if ($opt_verbose){
	print STDERR "Reading tRNAscanSE structure outfile " . $opt_inFile . " and printing tRNA annotations..";
}

open (GENETAB, '<', $opt_geneTab) || die "Unable to open file $opt_geneTab for reading\n$!\n";
while ($line= <GENETAB>){
	if ($line=~ m/^(\S+)\t(\d+)\t(\d+)\t(\-?1)\t(\S+)/){
		$geneFeat{$1}= [$2, $3, $4, $5];
	}
}
close (GENETAB);

if ($opt_verbose){
	print STDERR 'Captured ' . scalar(keys %geneFeat) . ' features from file: ' . $opt_geneTab . "\n";
}

open (RBSOUT, '<', $opt_inFile) || die "Unable to open file $opt_inFile for reading\n$!\n";
while ($line= <RBSOUT>){
	if ($line =~ m/^\s*(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)/){
		$locus_tag= $1;
		$new_start= $2;
		$stop= $3;
		$RBS_pattern= $4;
		$RBS_start= $5;
		$new_start_codon= $6;
		$start_shift= $7;
		$old_start_codon= $8;
		$old_start= $9;
		if ($old_start > $stop){
			$strand= -1;
		}
		else {
			$strand= 1;
		}
		if ($geneFeat{$locus_tag}[3]=~ m/\w+_\w+/){
			$geneFeat{$locus_tag}[3]= 'NULL';
		}
		if ($RBS_pattern!~ m/^\-+$/){
			if ($old_start ne $new_start){
				print STDERR '# ' . $locus_tag . ' at ' . $old_start . ' with predicted new start codon at ' . $new_start . ' (' . $old_start_codon . '->' . $new_start_codon . ')' ."\n";
				printRBSfeat(\*STDERR, $locus_tag, $RBS_start, $RBS_pattern, $strand, $geneFeat{$locus_tag}[3], 'AGGAG');
			}
			else {
				printRBSfeat(\*STDOUT, $locus_tag, $RBS_start, $RBS_pattern, $strand, $geneFeat{$locus_tag}[3], 'AGGAG');
			}
		}
	}
}
close (RBSOUT);

if ($opt_verbose){
	print STDERR "DONE\n";
}


exit (0);

