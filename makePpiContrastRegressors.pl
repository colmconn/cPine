#!/Volumes/opt/local/bin/perl

use strict;
use warnings;
use English;
use POSIX;
use Getopt::Long;

# use Scalar::Util 'reftype';
# use File::Basename;
# use List::MoreUtils qw/uniq/;
# use Data::Dumper;

## the locaiton where the data archive is stored. If you change this
## you may need to alter findTrainingSubjectDataFile and
## checkSubjectBehavioralDataDirectoryExists
my $MDD_ROOT='/Volumes/PROMISEPEGASUS/yangdata';

## the root of the direcotry hierarchy where the Pine regressors are to be stored
my $REGRESSOR_ROOT="$MDD_ROOT/cPine/data/regressors/";

## used for the construction of the 1deval commands
my @letters=qw(a b c d e f g h i j k l m n o p q r s t u v w x y z);

my %stimulusClassToColumnNumberMapping= (
					 ## These are the column indices of the regressors 
					 "instrA"          =>  1,
					 "instrB"          =>  2,
					 "instrC"          =>  3,
					 "instrD"          =>  4,
					 "instrR"          =>  5,
					 "fearful"         =>  6,
					 "fearfulRem"      =>  7,
					 "fearfulNotRem"   =>  8,
					 "fearfulOmission" =>  9,
					 "fixation"        => 10,
					 "happy"           => 11,
					 "happyRem"        => 12,
					 "happyNotRem"     => 13,
					 "happyOmission"   => 14,
					 "instructions"    => 15,
					 "neutral"         => 16,
					 "neutralRem"      => 17,
					 "neutralNotRem"   => 18,
					 "neutralOmission" => 19,
					 "sad"             => 20,
					 "sadRem"          => 21,
					 "sadNotRem"       => 22,
					 "sadOmission"     => 23
				);


my %contrastToStimuliMapping=(
    "fearfulVsHappy"            => ["fearful",    "happy"],
    "fearfulVsNeutral"          => ["fearful",    "neutral"],
    "fearfulVsSad"              => ["fearful",    "sad"],
    "happyVsNeutral"            => ["happy",      "neutral"],
    "happyVsSad"                => ["happy",      "sad"],
    "neutralVsSad"              => ["neutral",    "sad"],
    "allEmotiveVsNeutral"       => ["fearful",    "happy", "sad", "neutral"],
    "happyRemVsHappyNotrem"     => ["happyRem",   "happyNotrem"],
    "fearfulRemVsFearfulNotrem" => ["fearfulRem", "fearfulNotrem"],
    "neutralRemVsNeutralNotrem" => ["neutralRem", "neutralNotrem"],
    "sadRemVsSadNotrem"         => ["sadRem",     "sadNotRem"],
    "allRemVsAllNotrem"         => ["happyRem",   "happyNotrem", "fearfulRem", "fearfulNotrem", "sadRem", "sadNotrem", "neutralRem", "neutralNotrem"]
    );

##my @contrasts=("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad", "allEmotiveVsNeutral", "happyRemVsHappyNotrem",
##	       "fearfulRemVsFearfulNotrem", "neutralRemVsNeutralNotrem", "sadRemVsSadNotrem", "allRemVsAllNotrem");
## my @contrasts=("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad", "allEmotiveVsNeutral");
my @contrasts=("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad");

## number of fMRI volumes in the Pine task
my $runLength=860;
my $tr=2;
my $min_frac=0.6;

##print Dumper \%contrastToStimuliMapping;

my $subject;
my $help=0;
my $man=0;
GetOptions ('help|?' => \$help, 
	    'man' => \$man,
	    'subject=s' => \$subject) or pod2usage (-message => "Problem with the options specified" , -verbose => 2, -exitval => 1);

#	    'phase=s' => \$phase) or pod2usage (-message => "Problem with the options specified" , -verbose => 2, -exitval => 1);

pod2usage(1) if $help;
pod2usage(-verbose => 4) if $man;

if ( ! defined($subject) ) { 
  die ("You did not provide the ID of the subject on the command line.\n");
}

foreach my $contrast ( @contrasts ) {

  my @fileStack=();
  foreach my $emotion ( @{ $contrastToStimuliMapping{$contrast} } ) { 
    #print "$contrast: $emotion\n";

    my $binarizedTimingFile="${subject}.$emotion.ppi.binary.1D";
    push @fileStack, $binarizedTimingFile;

    my $awkCommand=
      "( cd $REGRESSOR_ROOT ; cat ${subject}.acquisition.onsetDuration.regressors.tab | sed 1d | awk \'{print \$$stimulusClassToColumnNumberMapping{$emotion}}\' |" .
      " grep -v 9999:9999 | tr \'\\n\' \' \' > ${subject}.$emotion.ppi.tab )";
    print "$awkCommand\n";
    system $awkCommand;

    my $timingCommand = "( cd $REGRESSOR_ROOT ; timing_tool.py -timing ${subject}.$emotion.ppi.tab -timing_to_1D $binarizedTimingFile -tr $tr -min_frac $min_frac -run_len $runLength )";
    print "$timingCommand\n";
    system $timingCommand;

    ## if the block of code below is uncommented the following block
    ## should be commented out, i.e, they are mutually exclusive
    ## this command converts the 1/0 binary files created by the timingCommand to be 1/-1
    
    ## start of block to be commented or uncommented 

    # my $conversionCommand = "( cd $REGRESSOR_ROOT ; 1deval -a ${binarizedTimingFile} -expr \"a - not(a)\" > ${binarizedTimingFile}.new ; mv -f ${binarizedTimingFile}.new $binarizedTimingFile )";
    # print "$conversionCommand\n";
    # system $conversionCommand;

    ## end of block to be commented or uncommented 

    # #rmdir "${subject}.$emotion.tab";
  }

  ## uncomment the following code to generate contrast regressors that are for example 1 for fearful and -1 for sad
  
  ## start of block to be commented or uncommented 
  
  my $command="1deval";
  my $expr="";
  my $i;
  my $contrastRegressorFilename="${subject}.$contrast.contrast.ppi.binary.1D";
  
  for ($i=0; $i <= $#fileStack; $i=$i+1) {
    $command=join(" ", ( $command, "-" . $letters[$i], $fileStack[$i] ) );
  }
  if ($contrast eq "allEmotiveVsNeutral" ) {
    $expr="a+b+c+(-d)";
  } else {
    $expr="a+(-b)";
  }
  
  my $contrastGenerationCommand = "( cd $REGRESSOR_ROOT ; $command -expr \"$expr\" > $contrastRegressorFilename )";
  print "$contrastGenerationCommand\n";
  system $contrastGenerationCommand;
  @fileStack=();
  
  ## end of block to be commented or uncommented 
}
