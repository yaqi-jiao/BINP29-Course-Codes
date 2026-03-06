#!/usr/bin/perl -w

# Author: Bjorn Canback
# License: No license.
# Usage: see -h

use strict;
use Getopt::Long;
use File::Path qw/make_path/;
use File::Basename;
use File::Spec;
use Cwd;

my $programName = File::Basename::fileparse($0,('.pl'));
my $program = $programName . '.pl';
my $version = '1.1';
my @ARGVcopy = @ARGV;

my ($fastaFile,$gffFile,$basename,$directory,$feature,$attribute,
    $protein,$report,$changePhase,$length,$force,$quiet,$versionNumber);

Getopt::Long::Configure('no_ignore_case');
my $optionSuccess =
  GetOptions('i=s' => \$fastaFile,
	     'g=s' => \$gffFile,
	     'b=s' => \$basename,
	     'd=s' => \$directory,
	     'f=s' => \$feature,
	     'a=s' => \$attribute,
	     'c' => \$changePhase,
	     'p' => \$protein,
	     'r' => \$report,
	     'l' => \$length,
	     'F' => \$force,
	     'Q' => \$quiet,
	     'v' => \$versionNumber,
	     'h' => \&usage 
	   );
exit(__LINE__) if ! $optionSuccess;

# Print version
if ($versionNumber) {
    print STDERR "$program version $version\n";
    exit(__LINE__);
}

##### CHECK OPTIONS AND SET DEFAULTS ######

if (! $length) {
    if ( ! $fastaFile ) {
	print STDERR "ERROR: An input fasta file has to be specified by -i.\n";
	exit(__LINE__);
    }
    if (! -e $fastaFile) {
	print STDERR "ERROR: The input file $fastaFile doesn't exist.\n";
	exit(__LINE__);
    }
    if (-z $fastaFile) {
	print STDERR "ERROR: The input file $fastaFile is empty.\n";
	exit(__LINE__);
    }
}
if ( ! $gffFile ) {
    print STDERR "ERROR: A gff/gtf file has to be specified by -g.\n";
    exit(__LINE__);
}
if (! -e $gffFile) {
    print STDERR "ERROR: The input file $gffFile doesn't exist.\n";
    exit(__LINE__);
}
if (-z $gffFile) {
    print STDERR "ERROR: The input file $gffFile is empty.\n";
    exit(__LINE__);
}
if ($changePhase && ! $protein) {
    print STDERR "ERROR: If -c is set, also -p has to be set.\n";
    exit(__LINE__);
}

$feature ||= 'CDS';
$attribute ||= 'gene_id';
$basename ||= $programName;
$directory ||= Cwd::getcwd() if ! defined($directory);
if (-d $directory && $directory ne Cwd::getcwd()) {
    print STDERR "ERROR: The output directory already exists. See $program -h.\n";
    exit(__LINE__);
}
File::Path::make_path($directory) or die $! if $directory ne Cwd::getcwd();

if (! $force) {
    if (-e "$basename.fna") {
	print STDERR "ERROR: The output file $basename.fna already exists. Use -F if you want to overwrite.\n";
	exit(__LINE__);
    }
    if (-e "$basename.log") {
	print STDERR "ERROR: The output file $basename.log already exist. Use -F if you want to overwrite.\n"; 
	exit(__LINE__);
    }
    if (-e "$basename.faa") {
	print STDERR "ERROR: The output file $basename.faa already exists. Use -F if you want to overwrite.\n";
	exit(__LINE__);
    }
    if (-e "$basename.report") {
	print STDERR "ERROR: The output file $basename.report already exists. Use -F if you want to overwrite.\n";
	exit(__LINE__);
    }
}

# Create output files
my $ntFile  = File::Spec->catfile($directory,"$basename.fna");
my $logFile = File::Spec->catfile($directory,"$basename.log");
my $aaFile = File::Spec->catfile($directory,"$basename.faa") if $protein;
my $reportFile = File::Spec->catfile($directory,"$basename.report") if $report;

##### PARSE GFF/GTF FILE AND KEEP IN MEMORY ######

# gff definition from: www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml

# 0 seqname
# 1 source
# 2 feature
# 3 start
# 4 end
# 5 score
# 6 strand
# 7 phase
# 8 [attributes]
# 9 [comments]

my $gene;                    # Current gene name
my $countGenes = 0;          # Number of genes
my $countFeature = 0;        # Number of found features
my $countFoundScaffolds = 0; # Number of scaffolds with genes
my $warning;                 # Current warning message.
my $countAdjustedGenes = 0;  # GLOBAL. Number of phase adjusted genes
my @warning1 = ();           # Warning messages.
my @warning2 = ();           # Warning messages.
my @warning3 = ();           # GLOBAL. Warning messages.
my @line = ();               # Split on tab of current line
my @gene = ();               # GLOBAL. Gene names
my @scaffold = ();           # Scaffold names
my @strand = ();             # GLOBAL. Gene direction
my @firstPos = ();           # GLOBAL. Gene index as first element, first positions as values. 2D.
my @lastPos = ();            # GLOBAL. Gene index as first element, first positions as values. 2D.
my @lastPhase = ();          # GLOBAL. Phase for the first feature in reverse strand. 1D.
my @foundScaffold = ();      # GLOBAL. Scaffold index as first element, gene indexes as values. 2D.
my @matrix;
my %geneIndex = ();          # Maps gene name to array element number
my %scaffoldIndex = ();      # GLOBAL. Maps scaffold name to array element number


open(GFF,$gffFile) or die $!;

calculate_length() if $length; # Ends with exit

while (my $line = <GFF>) {
    chomp $line;
    @line = split("\t",$line);
    if ($line[2] eq $feature) {

	$line[8] =~ m/$attribute[ =]\"*(.*?)\"*;/;
	$gene = $1;

	# Check gene id
	if (! $gene) {
	    $warning = "Line $. contained the feature $feature but not the attribute $attribute set by -a.";
	    push(@warning1,$warning);
	    next;
	}

	# Check that positions are increasing
	if ($line[3] > $line[4]) {
	    print STDERR "ERROR: Positional error in line $. for gene $gene.\n";
	    exit(__LINE__);
	}

	# Make a scaffold index
	if (! defined($scaffoldIndex{$line[0]})) {
	    push(@scaffold,$line[0]);
	    $scaffoldIndex{$line[0]} = $countFoundScaffolds;
	    $countFoundScaffolds++;
	}

	# Make a gene index
	if (! defined($geneIndex{$gene})) {
	    push(@gene,$gene);
	    push(@strand,$line[6]);
	    $geneIndex{$gene} = $#gene;
	    push(@{$foundScaffold[$scaffoldIndex{$line[0]}]},$geneIndex{$gene});
	    $countGenes++;
	}

	# Check that strands are not conflicting within the same gene
	if ($strand[$geneIndex{$gene}] ne $line[6]) {
	    print STDERR "ERROR: Gene $gene have conflicting directions.\n";
	    exit(__LINE__);
	}

	$countFeature++;

	if ($line[6] eq '+') {
	    if (! defined($firstPos[$geneIndex{$gene}])) {
		$line[7] = 0 if $line[7] eq '.';
		if ($line[7] == 0) {
		    push(@{$firstPos[$geneIndex{$gene}]},$line[3] - 1);
		} elsif ($line[7] == 1) {
		    push(@{$firstPos[$geneIndex{$gene}]},$line[3] - 1 + 2);
		} elsif ($line[7] == 2) {
		    push(@{$firstPos[$geneIndex{$gene}]},$line[3] - 1 + 1);
		} else {
		    print STDERR "ERROR: The $feature feature has a not supported phase ($line[7] in line $..\n";
		    exit(__LINE__);
		}
		push(@{$lastPos[$geneIndex{$gene}]},$line[4] - 1);
	    } else {
		push(@{$firstPos[$geneIndex{$gene}]},$line[3] - 1);
		push(@{$lastPos[$geneIndex{$gene}]},$line[4] - 1);
	    }
	} elsif ($line[6] eq '-') {
	    # For later use (overwrites):
	    push(@{$firstPos[$geneIndex{$gene}]},$line[3] - 1);
	    push(@{$lastPos[$geneIndex{$gene}]},$line[4] - 1);
	    if ($line[7] eq '.') {
		$lastPhase[$geneIndex{$gene}] = 0;
	    } else {
		$lastPhase[$geneIndex{$gene}] = $line[7];
	    }
	} else {
	    print STDERR "ERROR: The $feature feature has a not supported strand ($line[6] in line $..\n";
	    exit(__LINE__);
	}
    }
}
close(GFF) or die $!;

##### PARSE SCAFFOLD, ONLY KEEP CURRENT SCAFFOLD IN MEMORY ######

my $foundScaffold = 0;      # Boolean. Only read in sequence if true
my $countScaffolds = 0;     # Number of scaffolds in fasta file
my $scaffoldSequence = '';  # Current scaffod sequence
my $scaffold;               # Current scaffold name
my $oldScaffold;            # GLOBAL. The most recent scaffold name
my @scaffoldSequence = ();  # GLOBAL. Contain segments (lines) of current scaffold sequence
my @softMask = ();          # GLOBAL. An array of the fractions of soft masked nucleotides.
my %softWarning = ();       # GLOBAL. Gene index as key, 1 as value if internal stop codons.

open(GENOME,$fastaFile) or die $!;
open(NTOUT,">$ntFile") or die $!;
open(AAOUT,">$aaFile") or die $! if $aaFile;

while (my $line = <GENOME>) {
    chomp $line;
    if ($line =~ m/^>/) {

	$countScaffolds++;
	($scaffold = $line) =~ s/\s.*//;
	$scaffold =~ tr/>//d;
	$foundScaffold = 1; # Initially true

	# Check if scaffold contains genes
	if (! defined($scaffoldIndex{$scaffold})) {
	    $foundScaffold = 0;
	    next;
	}

	# Check if scaffold is found in the gff/gtf file
	if (! defined($foundScaffold[$scaffoldIndex{$scaffold}])) {
	    print STDERR "ERROR: Scaffold $scaffold not found in the gff/gtf file.\n";
	    exit(__LINE__);
	}

	# Output
	parse_and_print() if $countScaffolds > 1 && $foundScaffold == 1;
	$oldScaffold = $scaffold;
	@scaffoldSequence = ();
    } else {
	if ($foundScaffold == 1) {
	    next if $line !~ m/./;
	    $line =~ s/\s//g;
	    if ($line !~ m/^[a-zA-Z]+$/) {
		$warning = "Line $. contains illegal characters: $line";
		push(@warning2,$warning);
	    }
	    push(@scaffoldSequence,$line);
	}
    }
}

# Print last scaffold and close files
parse_and_print();
close(GENOME) or die $!;
close(NTOUT) or die $!;
close(AAOUT) or die $! if $aaFile;

##### PRINT SOFTMASK REPORT #####
if ($report) {
    my @out;
    open(REPORT,">$reportFile") or die $!;
    print REPORT "#Fraction of soft masked nucleotides.\n";
    for (my $i = 0; $i < scalar(@gene); $i++) {
	@out = ();
	push(@out,$gene[$i],$softMask[$i]);
	if (defined($softWarning{$i})) {
	    push(@out,'inframe-stop-codons');
	} else {
	    push(@out,'');
	}
	print REPORT join("\t",@out), "\n";
    }
    close(REPORT) or die $!;
}

##### PRINT LOG #####

open(LOG,">$logFile") or die $!;
print LOG "INFO: Program was run as $program ", join(' ',@ARGVcopy), "\n";
my $log = <<HERE;

INFO: The gff or gtf file $gffFile has successfully been parsed.
      There were $countFoundScaffolds scaffolds containing genes.
      The scaffolds contained $countGenes genes.
      The genes contained the feature $feature $countFeature times.

INFO: The scaffold/genome file $fastaFile was successfully parsed.
      There were $countScaffolds scaffolds. This number may be higher than the one above.

INFO: Following files were output in the directory
      $directory:
      $basename.fna (fasta file containing the genes)
      $basename.log (log file)
HERE

$log .= "      $basename.faa (fasta file containing translated genes)\n" if $aaFile;
$log .= "      $basename.report (report on fraction of soft masked nucleotides per gene)\n" if $report;
print STDERR $log, "\n" if ! $quiet;
print LOG $log;
my $countWarnings = scalar(@warning1) + scalar(@warning2) + scalar(@warning3);
if ($countWarnings > 0) {
    print STDERR "NOTICE: There are $countWarnings warnings in the log file.\n";
    print STDERR "        Reading frames have been adjusted for $countAdjustedGenes genes (with the use of -c option).\n\n" if $changePhase;
    print LOG "\nWARNINGS:\n\n";
    print LOG join("\n",@warning1), "\n" if scalar(@warning1) > 0;
    print LOG join("\n",@warning2), "\n" if scalar(@warning2) > 0;
    print LOG join("\n",@warning3), "\n" if scalar(@warning3) > 0;
}

##### SUBROUTINES #####

sub parse_and_print {

    my $geneSequence = '';      # Current gene sequence
    my $frame2;                 # Gene sequence counted from frame 2
    my $frame3;                 # Gene sequence counted from frame 2
    my $geneLength = 0;         # Current gene length
    my $lowQuality = '';        # Current gene sequence with soft masked nucleotides removed
    my $scaffoldLength = 0;     # Length of scaffold
    my $aa;                     # Amino acid sequence
    my $aa2;                    # Amino acid sequence for phase 2 if -c is used
    my $aa3;                    # Amino acid sequence for phase 3 if -c is used
    my $segment = '';           # Current segment sequence
    my $warning;                # Current warning message.
    my $adjusted = 0;           # If phase is adjusted
    my @geneSequence = ();      # Contain segments of current gene sequence

    $scaffoldSequence = join('',@scaffoldSequence);
    $scaffoldLength = length($scaffoldSequence);

    foreach my $geneIndex (@{$foundScaffold[$scaffoldIndex{$oldScaffold}]}) {
	$geneSequence = '';
	@geneSequence = ();
	for (my $i = 0; $i < scalar(@{$firstPos[$geneIndex]}); $i++) {

	    # Test if positions are found in scaffolds
	    if ($firstPos[$geneIndex][$i] < 0 ||
		$firstPos[$geneIndex][$i] > $scaffoldLength) {
		print STDERR "ERROR: First position of $firstPos[$geneIndex][$i] in $feature of gene $gene[$geneIndex] is outside the boundaries of scaffold $oldScaffold.\n";
		exit(__LINE__);
	    }
	    if ($lastPos[$geneIndex][$i] < 0 ||
		$lastPos[$geneIndex][$i] > $scaffoldLength) {
		print STDERR "ERROR: Last position of $lastPos[$geneIndex][$i] in $feature of gene $gene[$geneIndex] is outside the boundaries of scaffold $oldScaffold.\n";
		exit(__LINE__);
	    }
	    $segment = substr($scaffoldSequence,$firstPos[$geneIndex][$i],$lastPos[$geneIndex][$i]-$firstPos[$geneIndex][$i]+1);
	    push(@geneSequence,$segment);
	}

	$geneSequence = join('',@geneSequence);

	if ($strand[$geneIndex] eq "-") { # Complement strand
	    $geneSequence = reverse($geneSequence);
	    if ($lastPhase[$geneIndex] == 1) {
		$geneSequence = substr($geneSequence,2);
	    } elsif ($lastPhase[$geneIndex] == 2) {
		$geneSequence = substr($geneSequence,1);
	    }
	    $geneSequence =~ tr/acgtuACGTU/tgcaaTGCAA/;
	}

	$geneLength = length($geneSequence);
	if ($report) {
	    ($lowQuality = $geneSequence) =~ tr/ACGTU//d;
	    $softMask[$geneIndex] = sprintf("%.3f",length($lowQuality) / $geneLength);
	}

	if ($aaFile) {
	    $adjusted = 0;
	    $aa = translate_sequence($geneSequence);
	    if ( $aa =~ m/\*./ ) {
		$warning = "Gene $gene[$geneIndex] ($strand[$geneIndex] strand) contains internal stop codons.";
		$softWarning{$geneIndex} = 1 if $report;
		if ($changePhase) {
		    $frame2 = substr($geneSequence,1);
		    $aa2 = translate_sequence($frame2);
		    if ( $aa2 !~ m/\*./ ) {
			$warning .= " Phase 2 of sequence is used instead.";
			$geneSequence = $frame2;
			$geneLength--;
			$aa = $aa2;
			$countAdjustedGenes++;
			$adjusted = 2;
		    } else {
			$frame3 = substr($frame2,1);
			$aa3 = translate_sequence($frame3);
			if ( $aa3 !~ m/\*./ ) {
			    $warning .= " Phase 3 of sequence is used instead.";
			    $geneSequence = $frame3;
			    $geneLength -= 2;
			    $aa = $aa3;
			    $countAdjustedGenes++;
			    $adjusted = 3;
			}
		    }
		}
		push(@warning3,$warning);
	    }
	    if ($adjusted) {
		print AAOUT ">", join("\t",$gene[$geneIndex],'length=' . length($aa),'scaffold=' . $oldScaffold,'strand=' . $strand[$geneIndex],'adjustedToFrame=' . $adjusted), "\n", $aa, "\n";
	    } else {
		print AAOUT ">", join("\t",$gene[$geneIndex],'length=' . length($aa),'scaffold=' . $oldScaffold,'strand=' . $strand[$geneIndex]), "\n", $aa, "\n";
	    }
	}
	if ($adjusted) {
	    print NTOUT ">", join("\t",$gene[$geneIndex],'length=' . $geneLength,'scaffold=' . $oldScaffold,'strand=' . $strand[$geneIndex],'adjustedToFrame=' . $adjusted), "\n", $geneSequence, "\n";
	} else {
	    print NTOUT ">", join("\t",$gene[$geneIndex],'length=' . $geneLength,'scaffold=' . $oldScaffold,'strand=' . $strand[$geneIndex]), "\n", $geneSequence, "\n";
	}
    }
}

sub translate_sequence {

    my $sequence = $_[0];

    $sequence = uc($sequence);
    my($codon,@translated,$aa);

    for ( my $i = 0; $i < length($sequence) - 2; $i += 3) {

	$codon = substr($sequence,$i,3);

	# Sort according to the most used codons in E. coli (http://www.sci.sdsu.edu/~smaloy/MicrobialGenetics/topics/in-vitro-genetics/codon-usage.html)

	if ($codon eq 'CTG') {      # 5.2%
	    push(@translated,'L');
	} elsif ($codon eq 'GAA') { # 4.4%
	    push(@translated,'E');
	} elsif ($codon eq 'AAA') { # 3.8%
	    push(@translated,'K');
	} elsif ($codon eq 'GAT') { # 3.3%
	    push(@translated,'D');
	} elsif ($codon eq 'GCG') { # 3.2%
	    push(@translated,'A');
	} elsif ($codon eq 'GGC') { # 3.0%
	    push(@translated,'G');
	} elsif ($codon eq 'CAG') { # 2.9%
	    push(@translated,'Q');
	} elsif ($codon eq 'GGT') { # 2.8%
	    push(@translated,'G');
	} elsif ($codon eq 'ATT') { # 2.7%
	    push(@translated,'I');
	} elsif ($codon eq 'ATC') { # 2.7%
	    push(@translated,'I');
	} elsif ($codon eq 'ATG') { # 2.6%
	    push(@translated,'M');
	} elsif ($codon eq 'AAC') { # 2.6%
	    push(@translated,'N');
	} elsif ($codon eq 'CGT') { # 2.4%
	    push(@translated,'R');
	} elsif ($codon eq 'CCG') { # 2.4%
	    push(@translated,'P');
	} elsif ($codon eq 'GTG') { # 2.4%
	    push(@translated,'V');
	} elsif ($codon eq 'ACC') { # 2.4%
	    push(@translated,'T');
	} elsif ($codon eq 'GCC') { # 2.3%
	    push(@translated,'A');
	} elsif ($codon eq 'GAC') { # 2.3%
	    push(@translated,'D');
	} elsif ($codon eq 'CGC') { # 2.2%
	    push(@translated,'R');
	} elsif ($codon eq 'GCA') { # 2.1%
	    push(@translated,'A');
	} elsif ($codon eq 'GTT') { # 2.0%
	    push(@translated,'V');
	} elsif ($codon eq 'TTT') { # 1.9%
	    push(@translated,'F');
	} elsif ($codon eq 'GAG') { # 1.9%
	    push(@translated,'E');
	} elsif ($codon eq 'TTC') { # 1.8%
	    push(@translated,'F');
	} elsif ($codon eq 'GCT') { # 1.8%
	    push(@translated,'A');
	} elsif ($codon eq 'TAT') { # 1.6%
	    push(@translated,'Y');
	} elsif ($codon eq 'AAT') { # 1.6%
	    push(@translated,'N');
	} elsif ($codon eq 'AGC') { # 1.5%
	    push(@translated,'S');
	} elsif ($codon eq 'TAC') { # 1.4%
	    push(@translated,'Y');
	} elsif ($codon eq 'GTC') { # 1.4%
	    push(@translated,'V');
	} elsif ($codon eq 'TGG') { # 1.4%
	    push(@translated,'W');
	} elsif ($codon eq 'ACG') { # 1.3%
	    push(@translated,'T');
	} elsif ($codon eq 'CAA') { # 1.3%
	    push(@translated,'Q');
	} elsif ($codon eq 'CAT') { # 1.2%
	    push(@translated,'H');
	} elsif ($codon eq 'GTA') { # 1.2%
	    push(@translated,'V');
	} elsif ($codon eq 'ACT') { # 1.2%
	    push(@translated,'T');
	} elsif ($codon eq 'AAG') { # 1.2%
	    push(@translated,'K');
	} elsif ($codon eq 'TCT') { # 1.1%
	    push(@translated,'S');
	} elsif ($codon eq 'TTG') { # 1.1%
	    push(@translated,'L');
	} elsif ($codon eq 'CAC') { # 1.1%
	    push(@translated,'H');
	} elsif ($codon eq 'TTA') { # 1.0%
	    push(@translated,'L');
	} elsif ($codon eq 'CTT') { # 1.0%
	    push(@translated,'L');
	} elsif ($codon eq 'TCC') { # 1.0%
	    push(@translated,'S');
	} elsif ($codon eq 'GGG') { # 0.9%
	    push(@translated,'G');
	} elsif ($codon eq 'CTC') { # 0.9%
	    push(@translated,'L');
	} elsif ($codon eq 'CCA') { # 0.8%
	    push(@translated,'P');
	} elsif ($codon eq 'TCG') { # 0.8%
	    push(@translated,'S');
	} elsif ($codon eq 'GGA') { # 0.7%
	    push(@translated,'G');
	} elsif ($codon eq 'AGT') { # 0.7%
	    push(@translated,'S');
	} elsif ($codon eq 'TCA') { # 0.7%
	    push(@translated,'S');
	} elsif ($codon eq 'CCT') { # 0.7%
	    push(@translated,'P');
	} elsif ($codon eq 'TGC') { # 0.6%
	    push(@translated,'C');
	} elsif ($codon eq 'CGG') { # 0.5%
	    push(@translated,'R');
	} elsif ($codon eq 'TGT') { # 0.4%
	    push(@translated,'C');
	} elsif ($codon eq 'ATA') { # 0.4%
	    push(@translated,'I');
	} elsif ($codon eq 'CCC') { # 0.4%
	    push(@translated,'P');
	} elsif ($codon eq 'CGA') { # 0.3%
	    push(@translated,'R');
	} elsif ($codon eq 'CTA') { # 0.3%
	    push(@translated,'L');
	} elsif ($codon eq 'AGA') { # 0.2%
	    push(@translated,'R');
	} elsif ($codon eq 'AGG') { # 0.2%
	    push(@translated,'R');
	} elsif ($codon eq 'TAA') { # 0.2%
	    push(@translated,'*');
	} elsif ($codon eq 'TGA') { # 0.1%
	    push(@translated,'*');
	} elsif ($codon eq 'ACA') { # 0.1%
	    push(@translated,'T');
	} elsif ($codon eq 'TAG') { # 0.03
	    push(@translated,'*');
	} else {
	    push(@translated,"X");
	}
    }
    $aa = join("",@translated);
    return($aa);
}

sub usage {

    print STDERR <<HERE;

Description:   The program retrieves the nucleotide sequences from the 
               fasta input file based on the positions in the gff or gtf 
               file for a particular feature. If wanted these are 
               translated (only with the standard code). Attributes in the 
               gff-file are seperated by semicolons. The script looks for 
               the specified tag set by -a and hopefully retrieves the
               correct identifier.

Usage:         $program -i fasta_file -g gff_file [-b basename]
               [-d directory] [-f feature] [-a attribute] 
               [-p] [-s] [-l] [-F] [-Q] [-v] [-h]

Output files:  basename.fna, basename.log and optionally basename.faa,
               basename.report and basename.length.

Short options:

-i (string)    Name of fasta sequence input file.
-g (string)    Name of gff or gtf file.
-b (string)    Basename of output and log files. Default: $program.
-d (string)    Output directory. Default: Current working directory.
-f (string)    Sequence feature of interest, like exon or CDS. Default: CDS.
-a (string)    Tag used for identifying gene id. Default: gene_id.
-p (boolean)   Also output an amino acid fasta file.
-r (boolean)   Print a report of the proportion of soft masked (lower case)
               nucleotdies for each gene.
-c (boolean)   If internal stop codons are found in a protein (when -p is
               used), try the two other frames to see if the internal 
               stop codons disappears. If so, both the basename.fna and
               basename.aa will contain the adjusted gene and protein
               sequences respectively. The described problem is not rare 
               if the GFF file is produced by a gene prediction program
               like Genemark.
-l (boolean)   Just calculate length for the feature. Output is a two
               column list with gene id and length. No consistency checks
               are made. Phase information is not included.
-F (boolean)   Force potentially overwrite of existing file.
-Q (boolean)   Surpress output to screen.
-v (boolean)   Print version number and exit.
-h (boolean)   This help.

Example usage: $program -i laccaria.allmasked -g Lbicolor_BestModelsv1.0.gtf

HERE
exit(1);
}

sub calculate_length {

    my $gene;
    my @line;
    my @gene;
    my %gene;
    my %length;

    while (my $line = <GFF>) {
	chomp $line;
	@line = split("\t",$line);
	if ($line[2] eq $feature) {
	    $line[8] =~ m/$attribute[ =]\"*(.*?)\"*;/;
	    $gene = $1;
	    if (! defined($gene{$gene})) {
		push (@gene,$gene);
		$gene{$gene} = 1;
		$length{$gene} = 0;
	    }
	    $length{$gene} += $line[4] - $line[3] + 1;
	}
    }
    close(GFF) or die $!;

    $directory ||= Cwd::getcwd() if ! defined($directory);
    if (-d $directory && $directory ne Cwd::getcwd()) {
	print STDERR "ERROR: The output directory already exists. See $program -h.\n";
	exit(__LINE__);
    }
    File::Path::make_path($directory) or die $! if $directory ne Cwd::getcwd();
    my $lengthFile  = File::Spec->catfile($directory,"$basename.length");
    open(LENGTHOUT,">$lengthFile") or die $!;
    foreach my $gene (@gene) {
	print LENGTHOUT join("\t",$gene,$length{$gene}), "\n";
    }
    close(LENGTHOUT) or die $!;

    exit(0);
}

__END__

KEYWORDS gff, gtf, gff3, parse, assembly, genome
DEPENDENCIES 
