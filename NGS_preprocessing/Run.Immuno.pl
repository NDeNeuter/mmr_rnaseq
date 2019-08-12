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
# d: runfolder
getopts('d:', \%opts);  # option are in %opts

if (! defined($opts{'d'}) || $opts{'d'} eq '') {
	die("Specify run directory");
}

my $projectdir = $opts{'d'};
if (! -d $projectdir) {
	die("The specified project dir '$projectdir' does not exist");
}


#####################
# LOCATION BINARIES # 
#####################
my $javabin = "/opt/software/sun-jre-8/bin/java";
my $rbin = "/opt/software/R/3.2.1/bin/R";
my $pythonbin = "/opt/software/python2.7.5/bin/python";
my $trimmomaticjar = "/opt/software/Trimmomatic-0.36/trimmomatic-0.36.jar";
my $hisatbin = "/opt/software/hisat2-2.0.4/hisat2";

my $GRCh38_index_dir = "/opt/NGS/References/GRCh38/HISAT2_index/";
my $hisat2_base_indexname = "GRCh38_noalt_h2";
my $annotationfile = "/opt/NGS/References/GRCh38/GCF_000001405.31_GRCh38.p5_genomic.gff";
my $countfile = "readcounts.txt";
my $statsfile = "Statistics.csv";


#############################
# MAKE OUTPUT DIR STRUCTURE #
#############################
my $jodir = "$projectdir/Job_output";
my $tmpfilesdir = "$projectdir/tmp_files";
my $tmpbindir = "$projectdir/tmp_binaries";
my $resultsdir = "$projectdir/Results";

foreach my $dir ($jodir, $tmpfilesdir, $tmpbindir, $resultsdir) {
	if (! -d $dir) {
		mkdir $dir;
	}

}


###################
# OTHER VARIABLES #
###################
my $polyAscript = "/home/shared_data_immuno/bin/polyA_removal.pl";
my $countscript = "/home/shared_data_immuno/bin/CountTable.py";
my $statsscript = "/home/shared_data_immuno/bin/DeSeq2_mixed.r";

my $countjob = "$tmpbindir/Counts.sh";
my $statsjob = "$tmpbindir/Statistics.sh";

my $email = 'pieter.meysman@uantwerpen.be';
my $emailfile = "$tmpfilesdir/finishmail.txt";



############
# GET DATA #
############
my @files = split(/\n/, `cd $projectdir && ls *_R1_*fastq.gz`);
my @samples;
my @done;


##################
# START ANALYSIS #
##################
# Per Sample/input file
#     1. Trim with Trimmomatic
#     2. Trim polyA
#     3. Map with HISAT2
# Wait for all to finish
# On all files
#     4. Do counts
#     5. Do statistics


foreach my $file (@files) {
	$sample = `echo $file | cut -d_ -f1`;
	chomp($sample);
	#push(@samples, $sample);
	#print "Process sample $sample\n";
	
	$file =~ /(.*).fastq.gz/;
	my $stem = $1;
	#print "stem file: $stem\n";
	push(@samples, $stem);
        print "Process sample $stem\n";

	my $polyA_fastq ="$stem"."_polyA.fastq.gz";
	my $trimmomatic_fastq ="$stem"."_Trimmomatic.fastq.gz";
	my $samfile = "$stem".".sam";
	my $stemdone = $tmpfilesdir."/$stem.done";
	push(@done, $stemdone);

	my $trimmomaticjob = "$projectdir/tmp_binaries/$stem.Trimmomatic.sh";
	my $polyAtrimjob = "$projectdir/tmp_binaries/$stem.polyA-trim.sh";
	my $hisatjob = "$projectdir/tmp_binaries/$stem.HISAT2.sh";


	# 1. Make Trimmomatic trim
	open OUT, ">$trimmomaticjob";
	#&WriteHeader("OUT", $email, 6, 4, $projectdir, "Immuno.$stem.Trimmomatic", "Development", "");
	&WriteHeader("OUT", $email, 6, 6, $projectdir, "Immuno.$stem.Trimmomatic", "batch", "Immuno");
	# Write body
	&WriteCommand("OUT", "$javabin -jar $trimmomaticjar  SE -threads 6 -phred33 $projectdir/$file /tmp/$trimmomatic_fastq HEADCROP:20 SLIDINGWINDOW:4:15 MINLEN:30");
	&WriteCommand("OUT", "cp /tmp/$trimmomatic_fastq $tmpfilesdir");
	&WriteCommand("OUT", "rm /tmp/$trimmomatic_fastq");
	&WriteNextJobSubmit("OUT", $polyAtrimjob);
	close OUT;


	# 2. Make PolyA trim job file
	open OUT, ">$polyAtrimjob";
	#&WriteHeader("OUT", $email, 1, 2, $projectdir, "Immuno.$stem.polyAtrim", "Development", "");
	&WriteHeader("OUT", $email, 1, 2, $projectdir, "Immuno.$stem.polyAtrim", "batch", "Immuno");
	# Write body
	&WriteCommand("OUT", "perl $polyAscript -f $tmpfilesdir/$trimmomatic_fastq -o /tmp/$polyA_fastq");
	&WriteCommand("OUT", "cp /tmp/$polyA_fastq $tmpfilesdir/");
	&WriteCommand("OUT", "rm /tmp/$polyA_fastq");
	&WriteNextJobSubmit("OUT", $hisatjob);
	close OUT;


	# 3. Make HISAT2 mapping job file
	open OUT, ">$hisatjob";
	#&WriteHeader("OUT", $email, 1, 2, $projectdir, "Immuno.$stem.HISAT2", "Development", "");
	&WriteHeader("OUT", $email, 6, 8, $projectdir, "Immuno.$stem.HISAT2", "batch", "Immuno");
	# Write body
	&WriteCommand("OUT", "export HISAT2_INDEXES=$GRCh38_index_dir");
	&WriteCommand("OUT", "$hisatbin -p 6 -x $hisat2_base_indexname -U $tmpfilesdir/$polyA_fastq -S /tmp/$samfile --time");
	&WriteCommand("OUT", "cp /tmp/$samfile $tmpfilesdir/");
	&WriteCommand("OUT", "rm /tmp/$samfile");
	&WriteCommand("OUT", "echo 1 > $stemdone");
	close OUT;


	# submit first trim step
	my $submitcommand = "qsub $trimmomaticjob";
	print "$submitcommand\n";
	system($submitcommand);



}


#############################
## WAIT FOR JOBS TO FINISH ##
#############################
## torque/pbs does not have the -sync options, hence this workaround:
my $finished = 0;
while ($finished == 0) {
	$finished = 1;
	@queue = @done;
	foreach my $checkifdone (@queue) {
		if (! -e $checkifdone) {
			$finished = 0;
			sleep 60;
			last;
		}
	}
	#print "Wait...\n";
}

print "Proceed to count\n";

# 4. Make Counting job
open OUT, ">$countjob";
#&WriteHeader("OUT", $email, 1, 2, $projectdir, "Immuno.Count", "Development", "");
&WriteHeader("OUT", $email, 1, 6, $projectdir, "Immuno.Count", "batch", "Immuno");
&WriteCommand("OUT", "$pythonbin $countscript -i $tmpfilesdir -a $annotationfile -o /tmp/$countfile");
&WriteCommand("OUT", "cp /tmp/$countfile $resultsdir");
&WriteCommand("OUT", "rm /tmp/$countfile");
&WriteCommand("OUT", "echo 1 > $resultsdir/analysis.done");
#&WriteNextJobSubmit("OUT", $statsjob);
close OUT;


## 5. Make Statistics job
#open OUT, ">$statsjob";
#&WriteHeader("OUT", $email, 1, 2, $projectdir, "Immuno.Statistics", "Development", "");
##&WriteHeader("OUT", $email, 1, 2, $projectdir, "Immuno.Statistics", "batch", "Immuno");
#&WriteCommand("OUT", "$rbin $statsscript -i $resultsdir/$countfile -o /tmp/$statsfile");
#&WriteCommand("OUT", "cp /tmp/$statsfile $resultsdir");
#&WriteCommand("OUT", "rm /tmp/$statsfile");
#&WriteCommand("OUT", "echo 1 > $resultsdir/analysis.done");
#close OUT;


## All files made, submit count job and wait for analysis to finish
my $submitcommand = "qsub $countjob";
print "$submitcommand\n";
system($submitcommand);
while (! -e "$resultsdir/analysis.done") {
	print "Analysis done\n";
	sleep 60;
}


## Analysis done, send email to user
open OUT, ">$emailfile";
print OUT "To: $email\r\n";
print OUT "subject: Immuno Analysis Done\r\n";
print OUT "from: Immuno.Pipeline\r\n";
print OUT "The analysis is done, you know where to look for data:\r\n";
#print OUT "\r\n";
close OUT;

## send mail
system("sendmail -t < $emailfile");

#DONE



###############
# SUBROUTINES #
###############

## Function to write command to PBS script
sub WriteCommand {
	my ($filehandle, $print_command) = @_;
	if (tell($filehandle) == -1) {
		print "Outputfile is not open for writing! Exit script\n";
		exit;
	}
	print $filehandle "echo '$print_command'\n";
        print $filehandle "$print_command\n";
}


## Print out header to PBS script
sub WriteHeader {
	my ($filehandle, $tomail, $cpu, $mem, $dir, $jobname, $queue, $hpc_account) = @_;
	if (tell($filehandle) == -1) {
		print "Outputfile is not open for writing! Exit script\n";
		exit;
	}
        print $filehandle "#!/usr/bin/env bash\n";
        print $filehandle "#PBS -m a\n";
        print $filehandle "#PBS -M $tomail\n";
        print $filehandle "#PBS -d $dir\n";
        print $filehandle "#PBS -l nodes=1:ppn=$cpu,mem=$mem"."g\n";
        print $filehandle "#PBS -N $jobname\n";
        print $filehandle "#PBS -o $jodir/$jobname.o.txt.\$PBS_JOBID\n";
        print $filehandle "#PBS -e $jodir/$jobname.e.txt.\$PBS_JOBID\n";
        print $filehandle "#PBS -V\n";
        print $filehandle "#PBS -q $queue\n";
	if ($hpc_account) {
	        print $filehandle "#PBS -A $hpc_account\n";
	}
        print $filehandle "echo 'Running on : ' `hostname`\n";
        print $filehandle "echo 'Start Time : ' `date`\n";
        print $filehandle "echo 'Command:'\n";
        print $filehandle "echo '========'\n";
}

# Print out bash code for next PBS job submit (PBS script already in place)
sub WriteNextJobSubmit {
	my ($filehandle, $scriptsub) = @_;	
	if (tell($filehandle) == -1) {
		print "Outputfile is not open for writing! Exit script\n";
		exit;
	}

	print $filehandle "\n\n## submit next step\n";
	print $filehandle "JID=`qsub $scriptsub`\n";
	print $filehandle 'while  [[ ! $JID =~ ^[[:digit:]] ]] ; do'."\n";
	print $filehandle "  sleep 5\n";
	print $filehandle '  tmpFile=$(mktemp)'."\n";
	print $filehandle '  TRASH=`qstat -x 2>$tmpFile | grep '."$scriptsub`\n";
	print $filehandle '  # can be non-zero on grep-not-found and qstat-failure. only resumbit on valid qstat check.'."\n";
	print $filehandle '  if  [[ $? -ne 0 ]] ; then'."\n";
	print $filehandle '     if [ ! -s $tmpFile ]; then'."\n";
	print $filehandle "         JID=`qsub $scriptsub`\n";
	print $filehandle "     fi\n";
	print $filehandle "  fi\n";
	print $filehandle '  rm $tmpFile'."\n";
	print $filehandle "done\n";
}

