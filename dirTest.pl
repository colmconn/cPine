#!/Volumes/opt/local/bin/perl

use Data::Dumper;

## these variables are used to control the experimental phase being
## processed
my $ACQUISITION=0x8ABEABFA;
my $RECALL     =0x8ABEABAA;

my $MDD_ROOT='/Volumes/PROMISEPEGASUS/yangdata';

my $subjectId="396_A";
##my $subjectFolder=;

my $phase=$ACQUISITION;

my $subjectFolder = "$MDD_ROOT/$subjectId/$subjectId" . "_beh";

opendir(DIR, $subjectFolder) or die $!;

my @list = readdir(DIR);

foreach my $f (@list) {
  print "\$f = \"$f\"\n";
}

my @candidateFiles=();
#if ($phase == $ACQUISITION) {

print "Now trying to grep the list\n";
@candidateFiles= grep { /fac2sd-.*-\d*\.txt/		  # Begins with a period
			  && -f "$subjectFolder/$_" } @list;



#}
# elsif ($phase == $RECALL) {
#   @candidateFiles= grep { 
#     /rcl2sd-.*-[0-9]*\.txt/		  # Begins with a period
#       && -f "$subjectFolder/$_"	  # and is a file
#     } readdir(DIR);
# } else {
#   die "Unable to handlle that phase. Please check your command line and ensure you have selected either acuqisition or recall. Exiting\n";
# }
print "The length of the candidateFiles array is " . ( @candidateFiles ) . "\n";

print "Candidate Files are (in reverse order): ";
#@candidateFiles = sort {$b cmp $a} @candidateFiles;

print "@candidateFiles\n";

if ( scalar(@candidateFiles) > 1 ) { 
  print "Found more than one file returning only the first one\n";
}

close(DIR);
