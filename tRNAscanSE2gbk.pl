#!/usr/bin/perl

#===============================================================================
#   Author: Alejandro Manzano Marin
#
#   File: tRNAscanSE2gbk.pl
#   Date: 27-05-2013
#   Version: 1.0
#
#   Usage:
#      perl tRNAscanSE2gbk.pl --infile|-i <tRNAscanSEout structure file> [options]
#
#      Check out 'perl tRNAscanSE2gbk.pl -h' for short usage manual and info on the software.
#
#    Description: This program will take as input the structure file output of tRNAscanSE and will print
#                 the annotation of tRNA genes in genbank format.
#                 
#
#    Contact: Contact the author at alejandro.manzano@uv.es using 'tRNAscanSE2gbk: ' as
#             as begining for subject for bug reporting, feature request or whatever
#             (of course realted to the software).
#
#    COPYRIGHT: Copyright (C) 2013  Alejandro Manzano-Marin.
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


# Variable definition

## Define other variables
my $file= "";
my $line= "";
my $id= "";
my $flag= 0;
my %reads= ();
my $inHead= ">";
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
my $PROGRAMNAME= "tRNAscanSE2gbk";
my $VERSION= "1.0";

## Define options default values
my $opt_inFile= '';

my $opt_fasta= 0;

my $opt_verbose= 0;
my $opt_man= 0;
my $opt_help= 0;
my $opt_printVersion= 0;

## Define options hash
GetOptions(\%opts, 
	'infile|i=s' => \$opt_inFile, 
	'fasta!' => \$opt_fasta, 
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

tRNAscanSE2gbk

=head1 VERSION

tRNAscanSE2gbk v1.0

=head1 SYNOPSIS

perl tRNAscanSE2gbk.pl --infile|-i <tRNAscanSEout structure file> [-verbose|v] [-help|h] [-man] [-version]

=head1 DESCRIPTION

This program will take as input the structure file output of tRNAscanSE and will print the annotation of tRNA genes in genbank format

=head1 OPTIONS

=head2 INPUT

=over 8

=item B<-i> | B<--infile> <string> (mandatory)

File containing sequences to be filtered in FAST(A/Q) format.

=back

=head2 OTHER OPTIONS

=over 8

=item B<--fasta> <boolean> (default: 0)

If use, make a fasta file of the tRNAs.

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

If you find a bug please email me at C<< <alejandro_dot_manzano_at_uv_dot_es> >> so that I can keep making tRNAscanSE2gbk better.

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
		print STDERR "tRNAscanSE structure output file missing\n";
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

$c= '0001';
if ($opt_fasta){
	open (TRNAFASTA, ">$opt_inFile.fasta") || die "Unable to open file $opt_inFile.fasta for reading\n$!\n";
}

open (TRNASTRUCT, "<$opt_inFile") || die "Unable to open file $opt_inFile for reading\n$!\n";
while ($line= <TRNASTRUCT>){
	if ($line =~ m/^\S+\.trna\d+\s+\((\d+)\-(\d+)\)\s+Length:\s+\d+\s+bp/){
		$tRNAstart= $1;
		$tRNAEnd= $2;
		$flag{'pseudo'}= 0;
		$line= <TRNASTRUCT>;
		$line=~ m/^Type:\s+(\w+)\s+Anticodon:\s+(\w+)\s+at\s+\d+\-\d+\s+\((\d+)\-(\d+)\)\s+Score:\s+\d+\.\d+/;
		$anticodonThreeAA= $1;
		$anticodonString= $2;
		$anticodonStart= $3;
		$anticodonEnd= $4;
		$anticodonStringRNA= uc($anticodonString);
		$anticodonStringRNA=~ s/T/U/g;
		if ($tRNAstart < $tRNAEnd){
			$tRNAposString= $tRNAstart . '..' . $tRNAEnd;
			$anticodonAnnotString= '(pos:' . $anticodonStart . '..' . $anticodonEnd . '),aa:' . $anticodonThreeAA . ',seq:' . lc($anticodonString) . ')';
		}
		else {
			&switchVal(\$tRNAstart, \$tRNAEnd);
			&switchVal(\$anticodonStart, \$anticodonEnd);
			$tRNAposString= 'complement(' . $tRNAstart . '..' . $tRNAEnd . ')';
			$anticodonAnnotString= '(pos:complement(' . $anticodonStart . '..' . $anticodonEnd . '),aa:' . $anticodonThreeAA . ',seq:' . lc($anticodonString) . ')';
		}
		$codonString= &revComp(\$anticodonString);
		print STDOUT '     tRNA            ' . $tRNAposString . "\n";
		print STDOUT '                     /gene="' . $aaNames{$anticodonThreeAA}[0] . "\"\n";
		print STDOUT '                     /product="transfer RNA ' . $aaNames{$anticodonThreeAA}[1] . "\"\n";
		$line= <TRNASTRUCT>;
		print STDOUT '                     /anticodon="' . $anticodonAnnotString . "\"\n";
		print STDOUT '                     /inference="profile:tRNAscan-SE:1.3.1"' . "\n";
		if ($line=~ m/^Possible pseudogene/){
			print STDOUT '                     /pseudogene="unknown"' . "\n";
			print STDOUT '                     /note="putative tRNA-' . $anticodonThreeAA . " pseudogene\"\n";
		}
		else {
			print STDOUT '                     /inference="profile:RFAM:RF00005"' . "\n";
			print STDOUT '                     /note="tRNA-' . $anticodonThreeAA . '-' . $anticodonStringRNA . "\"\n";
		}
		
		
	}
	if ($opt_fasta && $line=~ m/^Seq:\s+(\w+)/){
		print TRNAFASTA ">tRNA-" . $anticodonThreeAA . "-" . $c . "\n" . $1 . "\n";
		$c++;
	}
	else {
		next;
	}
}
close (TRNASTRUCT);

if ($opt_fasta){
	close (TRNAFASTA);
}

if ($opt_verbose){
	print STDERR "DONE\n";
}


exit (0);

