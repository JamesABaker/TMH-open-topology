#---------------------------------------------------------------------------
# Singapore, 21st October 2011

# ACADEMIC SOFTWARE LICENSE for TMSOC
# ***********************************

# Copyright:
# Wing-Cheong Wong, Sebastian Maurer-Stroh, Georg Schneider, Frank Eisenhaber
# Bioinformatics Institute (BII) A*STAR Singapore

# This software is subjected to license restrictions. The software is
# provided as is. In its present form and without written consent of
# the copyright holder, the software is provided in accordance with GPL
# (www.gnu.org/licenses/gpl.html). Among other issues, this excludes any
# warranty and indemnity when using this software. For a commercial license,
# please approach the authors (via email to wongwc@bii.a-star.edu.sg).

# When publishing results with this software, please refer to:
# Wing-Cheong Wong, Sebastian Maurer-Stroh, Frank Eisenhaber, 2011,
# "Not all transmembrane helices are born equal: Towards the extension of
# the sequence homology concept to the membrane proteins", Biology Direct

# When observing bugs or strange behavior of the software, the user is
# encouraged to contact the authors.

#---------------------------------------------------------------------------
#!/usr/bin/perl
#use warnings;
#use strict;	# Forces variables to be declared

# Avoids awkward @INC path errors on some systems. Fixes Issue #5.
use FindBin;
use lib $FindBin::Bin;
use generateTMclassification;

my $i;
my $j;
my $k;
my @tmp1 = ();
my @tmp2 = ();
my @seqname = ();
my @FASTAseq = ();
my @TMsegments = ();
my $sequence = "";
my $segments = "";
my $resultsRef;
my @results;
my $newFASTAseq;
if (scalar(@ARGV)>=2 && scalar(@ARGV)<=3) {
	# Read FASTA sequences
	open (MYFILE1, $ARGV[0]);
	while (<MYFILE1>) {
		$line = $_;
		$line =~ s/\r//g;	# remove linefeed
		$line =~ s/\n//g;	# remove linefeed

		if ($line =~ /^>/) {
			if ($sequence ne "") {
				push(@FASTAseq, $sequence);
				$sequence = "";
			}
			@tmp1 = split(/\s+/, $line);
			push(@seqname, $tmp1[0]);
		}
		else { $sequence = $sequence.$line; }
	}
	push(@FASTAseq, $sequence);
	close(MYFILE1);

	# Read TM segments associated to FASTA sequences
	open (MYFILE2, $ARGV[1]);
	while (<MYFILE2>) {
		$line = $_;
		$line =~ s/\r//g;	# remove linefeed
		$line =~ s/\n//g;	# remove linefeed

		@tmp1 = split(/\s+/, $line);
		if (scalar(@tmp1) > 0) {
			push(@TMsegments, $line);
		}
	}
	close(MYFILE2);


	if (scalar(@FASTAseq)==scalar(@TMsegments) && scalar(@FASTAseq)>0) {

		for ($i=0; $i<scalar(@FASTAseq); $i++) {

			# Needs to offset the position values by -1; assumes first position starts at 1
			@tmp1 = split(/\s+/, $TMsegments[$i]);
			for ($j=0,$segment=""; $j<scalar(@tmp1); $j++) {
				@tmp2 = split(/\,/, $tmp1[$j]);
				if ($segment eq "") { $segment = ($tmp2[0]-1).",".($tmp2[1]-1); }
				else { $segment = $segment." ".($tmp2[0]-1).",".($tmp2[1]-1); }
			}

			# Classify the TM region(s)
			($resultsRef, $newFASTAseq) = generateTMclassification($segment, $FASTAseq[$i]);
			@results = @$resultsRef;
			if (uc($FASTAseq[$i]) eq uc($newFASTAseq)) {
				$newFASTAseq = "none";
			}

			# Output the results
			if ($ARGV[2] eq '-m') {
				print $seqname[$i]."\n".$newFASTAseq."\n\n";
			}
			else {
				print "1. TM segment(s) summary:\n";
				for ($j=0; $j<scalar(@results); $j++) {
					@tmp1 = split(/\;/,$results[$j]);
					print $results[$j]."\n";
				}
				print "2. Masked FASTA sequence:\n".$seqname[$i]."\n".$newFASTAseq."\n\n";
			}
		}
	}
	else {
		print "There are ".scalar(@FASTAseq)." sequences but ".scalar(@TMsegments)." associated TM segments. Please check input files(s).\n";
	}
}
else {
	print "Usage: perl TMSOC.pl [sequence.fasta] [TMsegments.txt] [options]\n";
	print "[options]: -m denotes output masked sequences only"
}
