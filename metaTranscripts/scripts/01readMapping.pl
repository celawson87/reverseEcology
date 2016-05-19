###############################################################################
# uncompReadMapping.pl
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This is a wrapper funtion to map metagenomics reads to reference genomes
# using BWA. It was developed to map reads from OMD-TOIL to our SAGs and MAGs.
################################################################################
#!/usr/bin/perl
use strict;
use warnings;

# Import packages
use File::Path qw(make_path);
use Parallel::ForkManager;

################################################################################
### User-defined files and folder structure
################################################################################

my $mtFolder = '../rawData';
my $genomeFolder = '../../refGenomes/concat';
my $mapFolder = '../bamFiles';
my $numThreads = 24;
my $numChildren = 1;

# Check if the output directory exists and create it if necessary
if (-d $mapFolder) {
  print "Output directory exists \n";
}
# If not, create it
else {
  print "Creating output directory \n";
  make_path($mapFolder);
}

my @genomeList = glob($genomeFolder.'/*.fna');
my @sampleList = glob($mtFolder.'/*.fastq');

################################################################################
### Step 2: Generate a list of commands to run
################################################################################

my @commandArray;

foreach my $samplePath (@sampleList) {
  foreach my $genomePath (@genomeList) {
    if ($genomePath =~ /.+\/(.+).fna/) {
      my $genome = $1;
      if ($samplePath =~ /.+\/(.+).fastq/) {
	my $sample = $1;	
       	push @commandArray, "bbmap.sh t=".$numThreads." in=".$mtFolder."/".$sample.".fastq outm=".$mapFolder."/".$sample."-".$genome.".sam ref=".$genomeFolder."/".$genome.".fna ambig=all nodisk sam=1.3"
      }
    }
  }
}

################################################################################
### Step 3: Simultaneously Index and Map Metatranscriptomes to Reference Genomes
################################################################################

my $pm = new Parallel::ForkManager($numChildren);

foreach my $command (@commandArray) {
  # Forks and returns the pid for the child:
  my $pid = $pm->start and next;

  system($command);
  
  $pm->finish; # Terminates the child process
}

$pm->wait_all_children;
print "Mapping complete!\n";
