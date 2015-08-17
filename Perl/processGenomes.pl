###############################################################################
# processGenomes.pl
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Converts a KBase model in .readable format to .TSV format for use with
# RE model-processing scripts.
################################################################################

#!/usr/bin/perl
use strict;
use warnings;
use Archive::Tar;
use File::Copy;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Path;

my $usage  = "Command sequence: perl readableToTSV.pl genomeDir rawModelDir \n";

my $inputDirectory  = shift or die $usage;
my $outputDirectory = shift or die $usage;

# Create an array of all *.fna files in the genome directory
my @modelFiles;

opendir (DIR, $inputDirectory) or die $!;
while (my $file = readdir(DIR)) {
  # Only include files ending in .fna
  next unless ($file =~ m/\.fna$/);
  # Remove file extension
  push (@modelFiles, $file);
}
closedir(DIR);

# Extract the tarball from KBase, move it, and clean up
# Extract
my $tarPath    = $outputDirectory.'/modelFiles.tar.gz';
my $tar        = Archive::Tar->new($tarPath);
$tar->extract_archive($tarPath, $outputDirectory);
# Move
my $sourceDir = 'tmp/narrative/modelDir';
dirmove($sourceDir, $outputDirectory);
# Clean
unlink($tarPath);
rmtree('tmp');

# Loop over the array of files and process each one
my $iter = 1;
my $numFiles = @modelFiles;

foreach my $file (@modelFiles) {

  # Hash to store visited metabs
  my %seenMetabs = ();

  # Trigger for processing metabs
  my $trigger = 0;
  
  # Grab just the model name
  $file =~ s/\.fna//;

  print "Processing ", $file, ", file ", $iter, " of ", $numFiles, "\n";
  
# Open input file. Create output file and directory.
  open INFILE, "<", $outputDirectory."/".$file.".readable" or die "Cannot open input file\n";
  mkdir "$outputDirectory/$file";
  open OUTFILE, ">", $outputDirectory."/".$file."/".$file."Compounds.tsv" or die "Cannot open output file\n";

  select STDOUT;

  while (<INFILE>) {
    # Split each line by tabs
    my @fields = split /\t/, $_;

    # If the first entry of the line matches '' then we have reached the metabolites
    if ($fields[0] =~ /modelcompounds \(formula/) {
      select STDOUT;
      print "Reached metab info\n";
      $trigger = 1;
      # Print header info
      #id	name	formula	charge	aliases
      select OUTFILE;
      print "id\tname\tformula\tcharge\n";
      next;
    }

    # If the first entry of the line matches '' then we're done
    if ($fields[0] =~ /modelreactions \(id/) {
      select STDOUT;
      print "Done processing metabs\n";
      $trigger = 0;
      next;
    }

    
    if (($trigger == 1) and exists($fields[1])) {
      # Iterate over string elements and replace _c0 or _e0
      s/_c0// for @fields;
      s/_e0// for @fields;	

      # Check if a metabolite has been seen before
      if (exists $seenMetabs{$fields[2]}) {
	next;
      } else {
	$seenMetabs{$fields[2]} = $fields[2];
      }

      # Print metabolite info in the proper order
      # id name modelcompounds (formula charge
      # If the charge is neutral, explicitly make it zero
      if ($fields[4]) {} else {$fields[4] = "0";}
      select OUTFILE;
      print $fields[1], "\t", $fields[2], "\t", $fields[0], "\t", $fields[4], "\n";
    }
    
  }

  $iter = $iter + 1;
  close INFILE;
  close OUTFILE;

  # Move .xml file and delete .readable file
  move("$outputDirectory/$file.xml", "$outputDirectory/$file/$file.xml");
  unlink("$outputDirectory/$file.readable");
  
}
