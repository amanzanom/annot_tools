#!/usr/bin/perl

#===============================================================================
#   Author: Alejandro Manzano Marin
#
#   File: MITOS2gbk.pl
#   Date: 05-11-2021
#   Version: 0.3
#
#   Usage:
#      perl MITOS2gbk.pl --infile|-i <MITOS_result_file> [options]
#
#      Check out 'perl MITOS2gbk.pl -h' for short usage manual and info on the software.
#
#    Description: This program uses the output from the MITOS v1.x pipeline
#                 (http://mitos.bioinf.uni-leipzig.de) and produces a complete annotation file (so far,
#                 only GenBank supported).
#                 
#
#    Contact: Contact the author at alejandro.manzano.marin@gmail.com using 'MITOS2gbk: ' as
#             as begining for subject for bug reporting, feature request or whatever
#             (of course realted to the software).
#
#    COPYRIGHT: Copyright (C) 2021  Alejandro Manzano-Marin.
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
use Data::Dumper;
use File::Basename;
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

sub loadGenCodeNCBI {
	my $genCodeFile= $_[0];

#	my $tempString= ();
	my $genCodeFile_fH= '';
	my $genCodeName= '';
	my $genCodeNum= 0;
	my @genCode_aa= ();
	my @genCode_startTERM= ();
#	my @genCode_TERM= ();
	my @genCode_1st= ();
	my @genCode_2nd= ();
	my @genCode_3rd= ();
	my $codonSeq= '';
	my $anticodonSeq= '';
	my %genCode= ();
	my $i= 0;
	
	# Define full amino acid info
	my %aaNames=(
		'A' => ['Ala', 'Alanine'],
		'R' => ['Arg', 'Arginine'],
		'N' => ['Asn', 'Asparagine'],
		'D' => ['Asp', 'Aspartic acid'],
		'C' => ['Cys', 'Cysteine'],
		'E' => ['Glu', 'Glutamic acid'],
		'Q' => ['Gln', 'Glutamine'],
		'G' => ['Gly', 'Glycine'],
		'H' => ['His', 'Histidine'],
		'I' => ['Ile', 'Isoleucine'],
		'L' => ['Leu', 'Leucine'],
		'K' => ['Lys', 'Lysine'],
		'M' => ['Met', 'Methionine'],
		'F' => ['Phe', 'Phenylalanine'],
		'P' => ['Pro', 'Proline'],
		'S' => ['Ser', 'Serine'],
		'T' => ['Thr', 'Threonine'],
		'W' => ['Trp', 'Tryptophan'],
		'Y' => ['Tyr', 'Tyrosine'],
		'V' => ['Val', 'Valine'],
		'U' => ['Sec', 'Selenocysteine'],
		'O' => ['Pyl', 'Pyrrolysine'],
		'*' => ['Stop', 'Stop']
	);

	open ($genCodeFile_fH, '<', $genCodeFile) || die "Unable to open file $genCodeFile for reading\n$!\n";
	while ($line= <$genCodeFile_fH>){
		if ($line=~ m/^\-\-/){
			next;
		}
		# Get the genetic code
		if ($line=~ m/^ \{$/){			
			
			# Get the genetic code name
			$line= <$genCodeFile_fH>;
			$line=~ m/^\s+name \"([^\"]+).+$/;
			$genCodeName= $1;

			while ($line!~/\" ,$/){
				$line= <$genCodeFile_fH>;
				$line=~ m/^ ([^\"]+)/;
				$genCodeName= $genCodeName . ' ' . $1;
			}
		
			# Get the genetic code number
			$line= <$genCodeFile_fH>;
			if ($line=~ m/^  name \"SGC\d+\" ,$/){ # Skip an extra line if double name line is encountered
				$line= <$genCodeFile_fH>;
			}
			$line=~ m/^  id (\d+) ,$/;
			$genCodeNum= $1;
			
			# Get the Genetic code amino acids
			$line= <$genCodeFile_fH>;
			$line=~ m/^  ncbieaa  \"(\S+)\",$/;
			@genCode_aa= split('', $1);
			
			# Get the Genetic code start codons
			$line= <$genCodeFile_fH>;
			$line=~ m/^  sncbieaa \"(\S+)\"$/;
#			$tempString= $1;
			@genCode_startTERM= split('', $1);
#			@genCode_TERM= split('', $1);
			
			# Get the Genetic code first codon nucleotide
			$line= <$genCodeFile_fH>;
			$line=~ m/^  \-\- Base1  (\S+)$/;
			@genCode_1st= split('', $1);
			
			# Get the Genetic code second codon nucleotide
			$line= <$genCodeFile_fH>;
			$line=~ m/^  \-\- Base2  (\S+)$/;
			@genCode_2nd= split('', $1);
			
			# Get the Genetic code third codon nucleotide
			$line= <$genCodeFile_fH>;
			$line=~ m/^  \-\- Base3  (\S+)$/;
			@genCode_3rd= split('', $1);
			
			# Fill the %genCode hash with the available genetic codes
			%{$genCode{$genCodeNum}}=();
#			print "Genetic code: $genCodeNum\n";
			for ($i=0; $i<scalar(@genCode_aa); $i++){
				$codonSeq= $genCode_1st[$i].$genCode_2nd[$i].$genCode_3rd[$i];
				$codonSeq=~ tr/ATGC/atgc/;
				$anticodonSeq= &revComp(\$codonSeq);
				$genCode{$genCodeNum}{$anticodonSeq}{'aa_oneLetter'}= $genCode_aa[$i];
				$genCode{$genCodeNum}{$anticodonSeq}{'aa_threeLetter'}= $aaNames{$genCode_aa[$i]}[0];
				$genCode{$genCodeNum}{$anticodonSeq}{'aa_fullName'}= $aaNames{$genCode_aa[$i]}[1];
				$genCode{$genCodeNum}{$anticodonSeq}{'as_startTERM'}= $genCode_startTERM[$i];
#				print "\t$anticodonSeq\n\t\t$genCode{$genCodeNum}{$anticodonSeq}{'as_startTERM'}\n";$pito=<STDIN>;
#				$genCode{$genCodeNum}{$anticodonSeq}{'as_TERM'}= $genCode_TERM[$i];
			}
		}
#		print Data::Dumper %genCode; $pito=<STDIN>;
	}
	close ($genCodeFile_fH);

	return(\%genCode);
}

sub printSeqGBK {
	my $sequence= $_[0];
	my $fileHandle= $_[1];
	
	my $seqLength= length($sequence);
	my $i= 0;
	my $j= 0;
	my $substring= '';
	my $flag_last= 0;
	
	print $fileHandle "ORIGIN\n";
	
	if ($charPerLine < 0){
		print STDERR "printFasta ERROR: Illegal number of characters per line, number must be >=0\n";
		return (1);
	}

	for ($i=0; $i<$seqLength; $i+=60){
		printf  $fileHandle sprintf("%9d", ($i+1));
		for ($j=$i; $j<($i+60); $j+=10){
			if (length(substr($sequence, $j, $seqLength-$j)) < 10){
				$substring= substr($sequence, $j, $seqLength-$j);
				$flag_last= 1;
			}
			else {
				$substring= substr($sequence, $j, 10);
			}
			print $fileHandle " " . $substring;
			if ($flag_last){
				last;
			}
		}
		print $fileHandle "\n";
	}
	
	print $fileHandle "\/\/";
	
	return (0);
}

# Variable definition

## Define other variables
my $file= '';
my $line= '';
my $id= '';
my $flag= 0;
my %qualifiers= ();
my $inHead= '>';
my $mitSeq= '';
my $dateString= '';
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
	'Pyl' => ['trnO', 'Pyrrolysine']	
);

my %names_gene=(
	'trnA'  => 'transfer RNA alanine',
	'trnR'  => 'transfer RNA arginine',
	'trnN'  => 'transfer RNA asparagine',
	'trnD'  => 'transfer RNA aspartic acid',
	'trnC'  => 'transfer RNA cysteine',
	'trnE'  => 'transfer RNA glutamic acid',
	'trnQ'  => 'transfer RNA glutamine',
	'trnG'  => 'transfer RNA glycine',
	'trnH'  => 'transfer RNA histidine',
	'trnI'  => 'transfer RNA isoleucine',
	'trnL1' => 'transfer RNA leucine',
	'trnL2' => 'transfer RNA leucine',
	'trnK'  => 'transfer RNA lysine',
	'trnM'  => 'transfer RNA methionine',
	'trnF'  => 'transfer RNA phenylalanine',
	'trnP'  => 'transfer RNA proline',
	'trnS1' => 'transfer RNA serine',
	'trnS2' => 'transfer RNA serine',
	'trnT'  => 'transfer RNA threonine',
	'trnW'  => 'transfer RNA tryptophan',
	'trnY'  => 'transfer RNA tyrosine',
	'trnV'  => 'transfer RNA valine',
	'trnU'  => 'transfer RNA selenocysteine',
	'trnO'  => 'transfer RNA pyrrolysine',
	'rrnS'  => 'small subunit ribosomal RNA',
	'rrnL'  => 'large subunit ribosomal RNA',
	'atp6'  => 'ATP synthase F0 subunit a',
	'atp8'  => 'ATP synthase F0 subunit 8',
	'cob'   => 'apocytochrome b',
	'cox1'  => 'cytochrome c oxidase subunit 1',
	'cox2'  => 'cytochrome c oxidase subunit 2',
	'cox3'  => 'cytochrome c oxidase subunit 3',
	'nad1'  => 'NADH dehydrogenase subunit 1',
	'nad2'  => 'NADH dehydrogenase subunit 2',
	'nad3'  => 'NADH dehydrogenase subunit 3',
	'nad4'  => 'NADH dehydrogenase subunit 4',
	'nad4L' => 'NADH dehydrogenase subunit 4L',
	'nad5'  => 'NADH dehydrogenase subunit 5',
	'nad6'  => 'NADH dehydrogenase subunit 6'
);

#my %geneticCode=('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11');

#my %{$geneticCode{'5'}}=(
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
my $PROGRAMNAME= "MITOS2gbk";
my $VERSION= "0.3";

#my $NCBIGENCODEFILE= "/media/manzanoa/DATA/DBs/NCBI/NCBI_genCode_v18-11-2016.prt";
my $NCBIGENCODEFILE= dirname(__FILE__) . "/NCBI_genCode_v4.6.txt";
my %genCode= %{&loadGenCodeNCBI($NCBIGENCODEFILE)}; # Load NCBI genetic code table
#print STDERR Data::Dumper->Dumper(%genCode);

## Define options default values
my $opt_inFile= '';
my $opt_inFASTA= '';
#my $opt_nctRNAfile= '';
#my $opt_ncrRNAfile= '';

my $opt_genCode= 5;
my $opt_noCDS=0;
my $opt_outFormat= 0;
my $opt_outPrefix= '';

my $opt_MITOSversion= 'rev6b33f95';

my $opt_verbose= 0;
my $opt_man= 0;
my $opt_help= 0;
my $opt_printVersion= 0;

## Define options hash
GetOptions(\%opts, 
	'infile|i=s' => \$opt_inFile, 
	'infasta|if=s' => \$opt_inFASTA, 
#	'nctRNAfile|nt:s' => \$opt_nctRNAfile, 
#	'ncrRNAfile|nr:s' => \$opt_ncrRNAfile, 
	'gencode|gc:i' => \$opt_genCode, 
	'noCDS!' => \$opt_noCDS, 
	'outformat|of!' => \$opt_outFormat, 
	'outprefix|o=s' => \$opt_outPrefix, 
	'MITOSversion:s' => \$opt_MITOSversion, 
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

MITOS2gbk

=head1 VERSION

MITOS2gbk v0.3

=head1 SYNOPSIS

perl MITOS2gbk.pl --infile|-i <result> [--infasta|--if] <sequence.fas> [--verbose|-v] [--help|-h] [--man] [--version]

=head1 DESCRIPTION

This program uses the output from the MITOS v1.x and v2.x pipeline (http://mitos.bioinf.uni-leipzig.de) and produces a complete annotation file (so far, only GenBank supported).

=head1 OPTIONS

=head2 INPUT

=over 8

=item B<-i> | B<--infile> <string> (mandatory)

MITOS "result" file containing the annotation, found inside the "mitfi-global" directory.

=item B<--if> | B<--infasta> <string> (optional)

MITOS "sequence.fas" file containing the FASTA-formatted sequence, found inside the "mitfi-global" directory.

=back

=head2 OUTPUT

=over 8

=item B<--gencode> | B<--gc> <integer> (optional)

Genetic code to use for annotations. If unsure which one to use see the NCBI website for genetic codes (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

=item B<--noCDS> <boolean> (default: 0)

If used, do not print out annotations for CDS features.

=item B<--of> | B<--outformat> <boolean> (default: gbk)

Choose (so far, only Genbank format supported).

=back

=head2 OTHER OPTIONS

=over 8

=item B<--MITOSversion> <string> (default: 917)

If use, make a fasta file of the tRNAs (default version up to date to 09-02-2017).

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

Alejandro Manzano Marin, C<< <alejandro_dot_manzano_dot_marin_at_gmail_dot_com> >>

=head1 BUGS

If you find a bug please email me at C<< <alejandro_dot_manzano_dot_marin_at_gmail_dot_com> >> so that I can keep making MITOS2gbk better.

=head1 COPYRIGHT

Copyright (C) 2017  Alejandro Manzano Marin.

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut


# Check options
#if (!$opt_inFile || ($opt_inFile && !$opt_noCDS && !$opt_nctRNAfile && !$opt_ncrRNAfile) ){
if (!$opt_inFile || !$opt_inFASTA){
	print STDERR "ERROR:\n";
	if (!$opt_inFile){
		print STDERR "MITOS result file missing\n";
	}
#	if ($opt_inFile && !$opt_noCDS && !$opt_nctRNAfile && !$opt_ncrRNAfile){
#		print STDERR "MITOS result file given, but all writing are off!\n";
#	}
	if (!$opt_inFASTA){
		print STDERR "MITOS FASTA file missing\n";
	}
	print STDERR "Please check usage manual\n";
	pod2usage(-exitval => 1,  -verbose => 1);
}


# Assign values to variables dependant on options


# Main program

## Read FASTA-fomratted input file
if ($opt_verbose){
	print STDERR 'Reading MITOS input sequence file ' . $opt_inFile . ' and start print of genbank-formatted file..';
}

### Get date for header
$dateString= localtime();

if ($dateString=~ m/^\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)/){
	$dateString= $2 . '-' . $1 . '-' . $3;
}
else {
	print STDERR "Shit! Seems like something in localtime() function in perl has changed, please check the call to this function in the script\n";
	pod2usage(-exitval => 1,  -verbose => 0);
}

### Read the FASTA-formatted file
open (MITOSFASTA, '<', $opt_inFASTA) || die "Unable to open file $opt_inFASTA for reading\n$!\n";
while ($line= <MITOSFASTA>){
	if ($line =~ m/^>(\S+)/){
		print STDOUT 'LOCUS       ' . $1;
		#"	15647 bp               linear 13-FEB-2018\n";
	}
	elsif ($line =~ m/^([ATGCN]+)/i){
		$mitSeq= $mitSeq . $1;
	}
}
close (MITOSFASTA);

print STDOUT '	' . length($mitSeq) . ' bp               linear ' . $dateString . "\n";

if ($opt_verbose){
	print STDERR "DONE\n";
}

## Read list of accepted reads
if ($opt_verbose){
	print STDERR 'Reading MITOS result file ' . $opt_inFile . ' and printing annotations to genbank-formatted file..';
}

#if ($opt_outFormat eq 'gbk'){
#	open (TRNAFASTA, '>', "$opt_outPrefix.$outFormat") || die "Unable to open file $opt_outPrefix.$outFormat for writing\n$!\n";
#}

print STDOUT "FEATURES             Location/Qualifiers\n";
open (MITOSRESULT, '<', $opt_inFile) || die "Unable to open file $opt_inFile for reading\n$!\n";
while ($line= <MITOSRESULT>){
	if ($line =~ m/^\S+\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+(\t\S+)?\t(\S+)\t(\S+\t)?\S+\t\S+\t\S+\t(\S+)$/){
		$feat_type= $1;
		$qualifiers{'gene'}= $2;
		$qualifiers{'inference'}= $3; # This one assumes MiTFI v0.1, MITOS revision is controlled by the "MITOSversion" option
		$feat_start= $4;
		$feat_end= $5;
		$feat_strand= $6;
		$unkField1= $7;
		$trna_codon= $8; # add check
		$unkField2= $9;
		$codon_start= $10; # add check
		
		if ($feat_type eq 'gene' && $opt_noCDS){
			next;
		}
		if ($feat_type eq 'gene'){
			$feat_type= 'CDS';
		}
		if ($qualifiers{'gene'} eq 'nad4l'){
			$qualifiers{'gene'}= 'nad4L';
		}
		if ($qualifiers{'inference'} eq 'mitfi'){
			$qualifiers{'inference'}= 'COORDINATES:profile:MiTFi:0.1';
		}
		elsif ($qualifiers{'inference'} eq 'mitos'){
			$qualifiers{'inference'}= "EXISTENCE:similar to AA sequence:MITOS:$opt_MITOSversion";
		}

		if ($feat_strand==-1){
			$feat_location= 'complement(' . ($feat_start+1) . '..' . ($feat_end+1) . ')';
		}
		else {
			$feat_location= ($feat_start+1) . '..' . ($feat_end+1);
		}
		
		# Check if fields 7 and 9 do not exist, if the do, then the result file comes form MITOS v2.x
		if (!$unkField1 && !$unkField2){
			$tRNA_anticodon= &revComp(\$trna_codon);
		}
		else {
			$tRNA_anticodon= $trna_codon;
		}
		$tRNA_anticodon=~ tr/ATGC/atgc/;
		if ($feat_type eq 'tRNA'){
			if ($feat_strand==-1){
				$anticodon_location= 'complement(' . ($feat_end+1-$codon_start-2) . '..' . ($feat_end+1-$codon_start) . ')';
			}
			else {
				$anticodon_location= ($feat_start+1+$codon_start) . '..' . ($feat_start+1+$codon_start+2);
			}
			$qualifiers{'anticodon'}= '(pos:' . $anticodon_location . ',aa:' . $genCode{$opt_genCode}{$tRNA_anticodon}{'aa_threeLetter'} . ',seq:' . $tRNA_anticodon . ')';
		}
		if ($codon_start eq '.'){
			$codon_start= 0;
		}
	
		### Print feature and qualifiers
		print STDOUT sprintf("     %-16s%s\n", $feat_type, $feat_location);
		print STDOUT sprintf("                     %s\n", "/gene=\"$qualifiers{'gene'}\"");
		print STDOUT sprintf("                     %s\n", "/product=\"$names_gene{$qualifiers{'gene'}}\"");
		print STDOUT sprintf("                     %s\n", "/inference=\"$qualifiers{'inference'}\"");
		if ($feat_type eq "CDS"){
			print STDOUT sprintf("                     %s\n", "/transl_table=$opt_genCode");
			print STDOUT sprintf("                     %s\n", '/codon_start=' . ($codon_start+1));
		}
		if ($feat_type eq "tRNA"){
			print STDOUT sprintf("                     %s\n", "/anticodon=$qualifiers{'anticodon'}");
		
			# Check for consistency between anticodon and supposed amino acid assignment
			$qualifiers{'gene'}=~ m/^trn(\w)/;
			$MiTFI_aa= $1;
			if ($genCode{$opt_genCode}{$tRNA_anticodon}{'aa_oneLetter'} ne $MiTFI_aa){
				print STDERR "Gencode: $opt_genCode\nAnticodon: $tRNA_anticodon\nOneLetter_AA: $genCode{$opt_genCode}{$tRNA_anticodon}{'aa_oneLetter'}\nMiTfI_AA: $MiTFI_aa\n";
				print STDOUT sprintf("                     %s\n", "/note=\"check, anticodon and MiTFI amino acid assignment differ\"");
			}
		}
	}
}
close (MITOSRESULT);

if ($opt_verbose){
	print STDERR "DONE\n";
}

### Print the sequence
printSeqGBK($mitSeq, \*STDOUT);


exit (0);

