#!/usr/bin/perl

use strict;
use List::Util qw(max);

sub usage();

my @array = ();

# no command line arguments?
if ( $#ARGV < 0  ) { 
  usage();
  exit;
}

print "Conormalizing @ARGV\n";
foreach my $argnum (0 .. $#ARGV) {

  open(LIST, '<', $ARGV[$argnum]) or die "Can't open results file ($ARGV[$argnum]): $!\n";
  {
    # put perl into slurp mode: this one line reads the entire contents of the file pointed to by LIST
    local $/;
    my $content=<LIST>;
    
    #print "####################################################################################################\n";
    #print "$ARGV[$argnum]\n";
    #print "$content\n";
    push @array, split /\n/, $content;
    #local $"="::";
    #print "array: @array\n";
  }
  close LIST;
}

my $max = max(@array);
#print "Max: $max\n";

foreach my $argnum (0 .. $#ARGV) {
  my $inputFile = $ARGV[$argnum];
  $ARGV[$argnum] =~ s/wavered\.1D/normalized\.1D/;
  
  open(my $waveredFile, '<', $inputFile) or die "Can't open results file $inputFile): $!\n";
  open(my $normalizedFile,  '>', $ARGV[$argnum]) or die "Can't open results file ($ARGV[$argnum]): $!\n";

  while (<$waveredFile>) {
    chomp;
    if ( $max != 0 ) { 
      print $normalizedFile $_ / $max . "\n";
    } else {
      print $normalizedFile $_ . "\n";
    }
  }
  
  close $waveredFile;
  close $normalizedFile;
}

sub usage() {

  print <<END
NAME

	$0 - normalises, to a maximum amplitude of 1, the wavered
	regressors in each of the files specified on the command line.

DESCRIPTION 

	$0  normalises, to a maximum amplitude of 1, the wavered
	regressors in each of the files specified on the command line.

	The program determines the maximum peak amplitued across all
	the files provided on the command line and uses this when
	scaling so that the maximum amplitude of all the output
	files is normalized to 1.

	The program expects the wavered files to have the text
	"wavered.1D" at the end of the filename. The normliased
	version of each filename will be the same as the input file
	*EXCEPT* that "wavered" is replaced with "normalized".

USAGE
	$0 regressor1

	Normalise the contents of regressor1 to have maximum amplitude of 1

	$0 regressor1 regressor2

	Conormalise regressor1 and regressor2 so that the maximum
	amplitude is 1 across both files.

	For example, to noramalize all wavered regressors for subject
	1221200502 the following command is appropriate:

	$0 1221200502*wavered.1D
	
END
}
