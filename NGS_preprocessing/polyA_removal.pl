#!/usr/bin/perl

###########
# MODULES #
###########
use Data::Dumper;
use Getopt::Std;
use Data::Dumper;
use XML::Simple;
use Cwd 'abs_path';

##########################
# COMMAND LINE ARGUMENTS #
##########################
my %opts;
# Options
# f: input FASTQ
# o: output FASTQ
getopts('f:o:', \%opts);  # option are in %opts
if (! defined($opts{'f'}) || $opts{'f'} eq '') {
        die("No input file specified");
}
my $inputfile = $opts{'f'};
if (! -e $inputfile) {
        die("The specified inputfile '$inputfile' does not exist");
}

if (! defined($opts{'o'}) || $opts{'o'} eq '') {
        die("No input file specified");
}
my $outputfile = $opts{'o'};


#########
# START #
#########
# 1. Trimm polyA if regex match
# 2. Print original read if no polyA trimmed
# 3. Filter out read if length polyA trimmed read < 20

my $unzipped_inputfile = `mktemp`;
chomp($unzipped_inputfile);
my $unzipped_outputfile = `mktemp`;
chomp($unzipped_outputfile);

open OUT, ">$unzipped_outputfile";

system("zcat $inputfile > $unzipped_inputfile");
open IN, $unzipped_inputfile or die "Cannot open $unzipped_inputfile";


# make sure length file is divisable by 4
my $lengthfile = `wc -l $unzipped_inputfile`;
chomp($lengthfile);
if ($lengthfile % 4 != 0) {
	print "Error in number of lines FASTQ file; not divisable by 4\n";
}

my $minlength = 20;
my $trim_count = 0;
my $skip_count = 0;
while (<IN>) {
	my $headerline = $_; 
	my $readline = <IN>;
	my $plusline = <IN>;
	my $qcline = <IN>;
	chomp($readline);
	chomp($qcline);
	
	my ($trimmed_readline, $trimmed_qcline);

	## Search for polyA tail (minimum 7 A's); non-greedy group to maximize tail
        if ($readline =~ /(.*?)A{7,}$/) {
		#print "polyA hit: $readline\n";
		$trimmed_readline = $1;
	}

	## Read was trimmed, print trimmed line
	if ($trimmed_readline) {
		$trim_count++;
		my $length_trimmed = length($trimmed_readline);		
		# Shorter than minimum length, skip
		if ($length_trimmed < $minlength) {
			$skip_count++;
			next;
		}

		$trimmed_qcline = substr($qcline, 0, $length_trimmed);

		print OUT $headerline;
		print OUT $trimmed_readline."\n";
		print OUT $plusline;
		print OUT $trimmed_qcline."\n";
	}

	## Not trimmed, print original line
	else {
		print OUT $headerline;
		print OUT $readline."\n";
		print OUT $plusline;
		print OUT $qcline."\n";
	}

}
 
close IN;
close OUT;

print "Number of reads trimmed: $trim_count; $skip_count discarded, smaller than $minlength bp\n\n";

my $zipped_outputfile = $unzipped_outputfile.".gz";
my $command = "gzip $unzipped_outputfile && mv $zipped_outputfile $outputfile";
print "$command\n";
system($command);

## clean up
$command = "rm $unzipped_inputfile";
print "$command\n";
system($command);
