#!/Volumes/opt/local/bin/perl

use strict;
use warnings;
use English;
use POSIX;
use File::Basename;
use File::Grep qw( fgrep );
use Getopt::Long;
use List::MoreUtils qw( uniq );

####################################################################################################
########## BEGINNING OF PROTOTYPES #################################################################
####################################################################################################

sub findSubjectBehavioralDataFile($$$);
sub checkSubjectBehavioralDataDirectoryExists($);

sub processHeader($);
sub processLogEntries($$$$);
sub processLogEntry($$$$$$);
sub processInstructionsLogEntry($$$$$);
sub processEpochLogEntry($$$$$);

sub printSummaryLine($$$);
sub printRegressorLines($$);
sub printInstructionsInEffectLine($);

sub makeInstructionsInEffectRegressorComponent($);

sub makeStimulusRegressorLine($$$$$$$$);

sub extractLogFrames($);

sub trim($);
sub ltrim($);
sub rtrim($);

####################################################################################################
########## END OF PROTOTYPES #######################################################################
####################################################################################################


####################################################################################################
########## BEGINNING OF CONSTNTS ##################################################################
####################################################################################################

my $ILLEGAL_OFFSET_VALUE=-10;
my $ILLEGAL_VALUE=-1;

## these variables are used to control which type of stimulus
## regressor line is created by the makeStimulusRegressorLine function
my $FEARFUL     =0xAFBAEBA8;
my $FIXATION    =0xAFBAEBB8;
my $HAPPY       =0xAFBAEBC8;
my $INSTRUCTIONS=0xAFBAEBD8;
my $NEUTRAL     =0xAFBAEBE8;
my $SAD         =0xAFBAEBF8;

## these variables are used to control the experimental phase being
## processed
my $ACQUISITION=0x8ABEABFA;
my $RECALL     =0x8ABEABAA;

## these are strings used as keys into the recall summary hash that
## holds the number of hits, misses correct rejections, and false
## alarms.
my $FEARFUL_RECALL_HITS="RecallFearfulHits";
my $HAPPY_RECALL_HITS="RecallHappyHits";
my $NEUTRAL_RECALL_HITS="RecallNeutralHits";
my $SAD_RECALL_HITS="RecallSadHits";

my $TOTAL_RECALL_HITS="TotalRecallHits";
my $TOTAL_RECALL_MISSES="TotalRecallMisses";
my $TOTAL_RECALL_CORRECT_REJECTIONS="TotalRecallCorrectRejections";
my $TOTAL_RECALL_FALSE_ALARMS="TotalRecallFalseAlarms";
my $TOTAL_RECALL_OMISSIONS="TotalRecallOmissions";

# my $TOTAL_ACQ_FEARFUL_RECALL_HITS="AcqTotalFearfulHits";
# my $TOTAL_ACQ_HAPPY_RECALL_HITS="AcqTotalHappyHits";
# my $TOTAL_ACQ_NEUTRAL_RECALL_HITS="AcqTotalNeutralHits";
# my $TOTAL_ACQ_SAD_RECALL_HITS="AcqTotalSadHits";

# my $TOTAL_ACQ_FEARFUL_RECALL_MISSES="AcqTotalFearfulMisses";
# my $TOTAL_ACQ_HAPPY_RECALL_MISSES="AcqTotalHappyMisses";
# my $TOTAL_ACQ_NEUTRAL_RECALL_MISSES="AcqTotalNeutralMisses";
# my $TOTAL_ACQ_SAD_RECALL_MISSES="AcqTotalSadMisses";

## these two constants are used to indicate whether a face presented
## at the recall phase was correctly remembered or not
my $REMEMBERED="Remembered";
my $NOT_REMEMBERED="NotRemembered";
my $OMISSION="Omission";
my $FALSE_ALARM="FalseAlarm";

## these constants hold the duration (in milliseconds) of the
## instruction and epoch blocks respectively. These durations are
## taken from Roberson-Nay R, McClure EB, Monk CS, et al. Increased
## amygdala activity during successful memory encoding in adolescent
## major depressive disorder: An FMRI study. Biol
## Psychiatry. 2006;60(9):966–973.

my $INSTRUCTIONS_DURATION=3000;
my $TRIAL_DURATION=4000;

# these two constants code the type of regressor that should be
# written. The first indicates that time pairs written in the
# regressor lines should be in the format onset:duration. The second
# indicates that the time pairs should be in the format onset:offset.
my $ONSET_DURATION_REGRESSOR=0;
my $ONSET_OFFSET_REGRESSOR=1;
my $ONSET_ONLY_REGRESSOR=2;


## the locaiton where the data archive is stored. If you change this
## you may need to alter findSubjectBehavioralDataFile and
## checkSubjectBehavioralDataDirectoryExists
my $MDD_ROOT='/Volumes/PROMISEPEGASUS/yangdata';

## the root of the direcotry hierarchy where the Pine regressors are to be stored
my $REGRESSOR_ROOT="$MDD_ROOT/cPine/data/regressors/";

# a flag to indicate whether debuging info should be written to stdout
# (a value of 1 tends to make it very verbose)
my $DEBUG=0;


####################################################################################################
########## END OF CONSTANTS ########################################################################
####################################################################################################


####################################################################################################
########## BEGINNING OF VARIABLES ##################################################################
####################################################################################################
# indicates the type of regressor that should be written to the
# regressor output files
my $regressorType = $ONSET_DURATION_REGRESSOR;

# indicates the experimental phase that should be processed: acquisition or recall
#my $phase;

## the ID of the subject
my $subjectId;

# the name of this file is stored in $subjectDataFile
my $subjectDataFile;

my $programName = $ARGV[0];
my $help=0;
my $man=0;
my $forceDuration;

my $firstStimulusPresentedAt=$ILLEGAL_VALUE;

## a hash to keep track of the onset times of each instruction block
## the recall phase of the expriment
my %globalRecallInstructionsOnsetTimes;

# a hash to keep track of whihc actor presented which emotion for the
# subjects
my %globalEmotionToActorMap;

## an array of arrays: each of the subarrays consists of a triplet as
## follows: (trialNumber, targetId, targetEmotion) target is the actor
## depicting the emotion
my @globalRecallTrialsArray=();
####################################################################################################
########## END OF VARIABLES ########################################################################
####################################################################################################

# my $MDD_ROOT=$ENV{'MDD_ROOT'};
# if ( ! defined($MDD_ROOT) ) { 
#   die "Could not find the MDD_ROOT environment variable. Cannot continue. Exiting.\n";
# } else {
#   print "MDD_ROOT is set to $MDD_ROOT\n";
# }

if ( $#ARGV + 1 < 1 ) { 
  pod2usage(1);
}

$forceDuration=1;
my $regressor="onsetDuration";

GetOptions ('debug' => \$DEBUG, 
	    'help|?' => \$help, 
	    'man' => \$man,
	    'regressor=s' => \$regressor,
	    'subject=s' => \$subjectId) or pod2usage (-message => "Problem with the options specified" , -verbose => 2, -exitval => 1);

#	    'phase=s' => \$phase) or pod2usage (-message => "Problem with the options specified" , -verbose => 2, -exitval => 1);

pod2usage(1) if $help;
pod2usage(-verbose => 4) if $man;

if ( ! defined($subjectId) ) { 
  print ("You did not provide the ID of the subject on the command line.\nWill process all subjects available in $MDD_ROOT\n.");
}

if ($regressor =~ m/onsetOf\w*/ ) {
  #$regressorType = $ONSET_OFFSET_REGRESSOR;
  #print "Creating onset:offset regressors\n";
  die "Sorry can't create offset:offset regressors\n";
} elsif ($regressor =~ m/onsetOn\w*/ ) {
  $regressorType = $ONSET_ONLY_REGRESSOR;
  print "Creating onset only regressors\n";
} elsif ($regressor =~ m/onsetD\w*/ ) {
  $regressorType = $ONSET_DURATION_REGRESSOR;
  print "Creating onset:duration regressors\n";
} else {
  pod2usage(-message => "Sorry I don't know how to create regressors of type $regressorType\n",
	    -verbose => 3);
}

if ($DEBUG) {
  print "DEBUG is $DEBUG\n";
  print "regressorType is $regressorType\n";
  print "forceDuration is " . ( defined($forceDuration) ? "$forceDuration" : "undefined" ) . "\n";
}


if ( ! -d $REGRESSOR_ROOT ) { 
  print "$REGRESSOR_ROOT does not exist. Creating it\n";
  mkdir $REGRESSOR_ROOT;
}

my $subjectFolder=checkSubjectBehavioralDataDirectoryExists($subjectId);
print "Got subjectFolder: $subjectFolder\n";
if ( $subjectFolder ne "0" ) { 

  #my @phases = ( $RECALL );
  my @phases = ( $RECALL, $ACQUISITION );

  ## a hash to keep track of phaces that were recalled or not during
  ## the recall phase of the expriment
  my %recallFacesHash = ();

  my %recallSummaryHash = (
			   $FEARFUL_RECALL_HITS => 0,
			   $HAPPY_RECALL_HITS => 0,
			   $NEUTRAL_RECALL_HITS => 0,
			   $SAD_RECALL_HITS => 0,
			   $TOTAL_RECALL_HITS => 0,
			   $TOTAL_RECALL_MISSES => 0,
			   $TOTAL_RECALL_CORRECT_REJECTIONS => 0,
			   $TOTAL_RECALL_FALSE_ALARMS => 0,
			   $TOTAL_RECALL_OMISSIONS => 0
			  );
  
  foreach my $phase ( @phases ) {
    %globalRecallInstructionsOnsetTimes = ();
    print "Processing the " . ( $phase == $ACQUISITION ? "acquisition" : "recall" ) . " experimental phase.\n";
    
    print "Attempting to locate the subject's behavioral data file.\n";
    print "Looking in: $subjectFolder\n";
    my $behavioralFile=findSubjectBehavioralDataFile($subjectId, $subjectFolder, $phase);
    if ( defined ($behavioralFile) ) {
      print "Selected subject behavioral file is: $behavioralFile\n";

      #my $subjectFolder = "$MDD_ROOT/$subjectId/$subjectId" . "_beh" . "/$behavioralFile";
      $behavioralFile = "$subjectFolder" . "/$behavioralFile";

      open(my $fh, "<:crlf", $behavioralFile) or die "Cannot open $behavioralFile: $!\n";
      # read file into an array
      ## my @pineData = <$fh>;

      #undef my $firstStimulusPresentedAt;
    
      processHeader($fh);

      ## this is where all the work of processing the log entries occurs

      my @logFrameEntries = @{ extractLogFrames($fh) };
      my $regressorLines=processLogEntries(\@logFrameEntries, \%recallSummaryHash, \%recallFacesHash, $phase);
    
      my $regressorFilename="$REGRESSOR_ROOT/$subjectId" . "." . ( $phase == $ACQUISITION ? "acquisition" : "recall" );
      if ($regressorType == $ONSET_OFFSET_REGRESSOR) {
	$regressorFilename = $regressorFilename . ".onsetOffset";
      } elsif ($regressorType == $ONSET_ONLY_REGRESSOR) {
	$regressorFilename = $regressorFilename . ".onsetOnly";
      } elsif ($regressorType == $ONSET_DURATION_REGRESSOR) {
	$regressorFilename = $regressorFilename . ".onsetDuration";
      }
      $regressorFilename = $regressorFilename . ".regressors.tab";
      my $instructionsInEffectRegressorFilenamePrefix="$REGRESSOR_ROOT/$subjectId" . "." . ( $phase == $ACQUISITION ? "acquisition" : "recall" );# . ".instructions.in.effect.regressors.tab";
    
      ## now print the results

      if ( $phase == $ACQUISITION ) {

	print "Emotion to actor mapping for this subject\n";
	for my $key ( sort keys %globalEmotionToActorMap ) {
	  my @value = uniq sort @{ $globalEmotionToActorMap{$key} };
	  $globalEmotionToActorMap{$key} = \@value;
	  print "$key => ";
	  print join(" ", @{ $globalEmotionToActorMap{$key}});
	  print " (# " . ( 0 +  @value ) . " )\n";
	}

	print "Writing regressors to $regressorFilename\n";
	print "Note that there will be 160 regressor lines for the emotion trials and 16 for the intstructions presentations, making a total of 176 lines\n";
	printRegressorLines($regressorLines, $regressorFilename);
	# print "Onset of instruction blocks\n";
	# for my $key ( sort keys %globalRecallInstructionsOnsetTimes ) {
	#   my $value = $globalRecallInstructionsOnsetTimes{$key};
	#   print "$key => @{ $value }\n";
	# }
	# print "Writing instructions in effect regressors to $instructionsInEffectRegressorFilenamePrefix\n";
	printInstructionsInEffectLine($instructionsInEffectRegressorFilenamePrefix);
      } elsif ( $phase == $RECALL ) {
    
    	print "Memory of faces summary (# " . ( 0 + keys %recallFacesHash  ) . " )\n";
    	for my $key ( sort keys %recallFacesHash ) {
    	  my $value = $recallFacesHash{$key};
    	  print "$key => $value\n";
    	}

      } ## end of if ( $phase == $RECALL )

      ## close the log file
      close ($fh);

    } else {
      print "Could not find a behavioral file for this subject\n";
    }## end of if ( defined ($behavioralFile) )
  } ## end of foreach $phase in ( $RECALL, $ACQUISITION )

  ## now we need to reverse the mapping from emotion to actor in the
  ## globalEmotionToActorMap
  my %actorToEmotionMap=();
  for my $emotion ( sort keys %globalEmotionToActorMap ) {
    my @actors = @{ $globalEmotionToActorMap{$emotion} };
    for my $actor ( @actors ) {
      $actorToEmotionMap{$actor} = $emotion;
    }
  }
  print "Actor to emotion mapping\n";
  for my $key ( sort keys %actorToEmotionMap ) {
    my $value = $actorToEmotionMap{$key};
    print "$key => $value\n";
  }
  
  ##print "globalRecallTrialArray (# ". ( 0 + @globalRecallTrialsArray ) . " ) :\n";
  for my $i (0 .. @globalRecallTrialsArray - 1 ) {
    my $targetId=$globalRecallTrialsArray[$i][0];
    my $targetEmotion=$globalRecallTrialsArray[$i][1];
    my $recalled = $recallFacesHash{$targetId} || "UNKN";
    ##print "$i: $targetEmotion $targetId $recalled\n";
    
    if ($recalled eq $REMEMBERED) {
      if (exists($actorToEmotionMap{$targetId}) ) {
	if ($actorToEmotionMap{$targetId}  eq "f") {
	  $recallSummaryHash{ $FEARFUL_RECALL_HITS } = $recallSummaryHash{ $FEARFUL_RECALL_HITS } + 1;
	  if ($DEBUG) {
	    print "### FEARFUL face CORRECTLY detected\n";
	  }
	} elsif ($actorToEmotionMap{$targetId}  eq "h") {
	  $recallSummaryHash{ $HAPPY_RECALL_HITS } = $recallSummaryHash{ $HAPPY_RECALL_HITS } + 1;
	  if ($DEBUG) {
	    print "### HAPPY face CORRECTLY detected\n";
	  }
	} elsif ($actorToEmotionMap{$targetId}  eq "n") {
	  $recallSummaryHash{ $NEUTRAL_RECALL_HITS } = $recallSummaryHash{ $NEUTRAL_RECALL_HITS } + 1;
	  if ($DEBUG) {
	    print "### NEUTRAL face CORRECTLY detected\n";
	  }
	} elsif ($actorToEmotionMap{$targetId}  eq "s") {
	  $recallSummaryHash{ $SAD_RECALL_HITS } = $recallSummaryHash{ $SAD_RECALL_HITS } + 1;
	  if ($DEBUG) {
	    print "### SAD face CORRECTLY detected\n";
	  }
	}
      } ## end of if (exists($actorToEmotionMap{$targetId}) ) {
    } ## end of if ($recalled eq $REMEMBERED) {
  } ## end of for my $i (0 .. @globalRecallTrialsArray - 1 ) {
  
  my $summaryFilename="$REGRESSOR_ROOT/summary.tab";
  
  if ( fgrep { /$subjectId/ } $summaryFilename ) {
    print "Found subject $subjectId in $summaryFilename. Not adding it again.\n"
  } else {
    printSummaryLine($subjectId, \%recallSummaryHash, $summaryFilename);
  }
  
  print "Recall summary\n";
  while ( my ($key, $value) = each(%recallSummaryHash) ) {
    print "$key => $value\n";
  }
} else { 
  print "Could not find subjects behavioral directory. Skipping.\n";
} ## end of if ( $subjectFolder ne "0" ) { 
 
####################################################################################################
########## Procedures are after this point #########################################################
####################################################################################################

sub processHeader($) {

  my $pineDataFile=shift;
  my $experimentDate="";
  my $experimentTime="";

  while ( <$pineDataFile> ) {
    chomp $_;
    if ( /\*\*\* Header Start \*\*\*/ .. /\*\*\* Header End \*\*\*/) {
      ## print "@@@@@@@@@@ Line: $_\n";
      if ( /SessionDate:\s([-0-9]*)/) {
	$experimentDate=$1;
      }
      if ( /SessionTime:\s([:0-9]*)/) {
	$experimentTime=$1;
      }
    } else {
      last;
    }
  }
  
  print "Experiment was conducted on $experimentDate at $experimentTime.\n";
}

sub processLogEntries($$$$) {

  my $logEntries=shift;
  my $recallSummaryHash_ref=shift;
  my $facesRecalledHash_ref=shift;
  my $phase=shift;

  my $regressorLine;
  my @regressorLinesArray=();
  
  print "Retrieved " . @$logEntries . " log entries.\n";
  print "Currently the final log entry is skipped. The number of lines in the regressor file will be one less than the number of indicated as retrieved above.\n";
  
  ## $instructionsInEffect can take on the values ABCDR and U. The
  ## latter is used as an initializing constant and should never
  ## survive past processing the behavioral file.
  my $instructionsInEffect="U";
  
  ## the -1 skips the last entry until i know what to do with it
  for my $i (0 .. $#$logEntries - 1 ) {
    if ($DEBUG) {
      print "####################################################################################################\n";
      print "### processLogEntries - Processing LogFrame $i\n";
    }
    ($instructionsInEffect, $regressorLine)=processLogEntry($recallSummaryHash_ref, $facesRecalledHash_ref, $instructionsInEffect, @{$logEntries}[$i], $i, $phase);
    push @regressorLinesArray, $regressorLine;
  }
  
  if ($DEBUG) {
    for my $i ( 0 .. $#regressorLinesArray ) {
      print "### processLogEntries - regressor line $i: " . $regressorLinesArray[$i] . "\n";
    }
  }

  return (\@regressorLinesArray);
}


## this function only process a single log entry at a time. Stop forgetting this, dumb ass!
sub processLogEntry($$$$$$) {

  my $recallSummaryHash_ref=shift;
  my $facesRecalledHash_ref=shift;
  my $instructionsInEffect=shift;
  my $logEntry=shift;
  my $logEntryNumber=shift;
  my $phase=shift;

  my $regressorLine = "**** ERROR default value assigned in processLogEntry. If you can see this in the text file some thing really bad has gone wrong. ****";
  ##
 

  ## $i indexes the line in the log entry
  for my $i ( 0 .. $#$logEntry ) {
    if ($DEBUG) {
      print "### processLogEntry - Processing LogFrame line: $i\n";
      print @{$logEntry}[$i] . "\n";
    }
    if ( @{$logEntry}[$i] =~ /\s*Procedure:\sInstr([ABCDR])/ ) {
      $instructionsInEffect=$1;
      if ($instructionsInEffect eq "R" && $phase != $RECALL) {
	die "Trying to process the recall phase of the experiment, but the instructions set in effect is not the recall set. Something is wrong with your command line arguments, the subject's behavioral data file, or both. This should never happen. Exiting\n";
      }
      if ($DEBUG) {
	print "### processLogEntry - INSTRUCTIONS detected\n";
	print "### processLogEntry - Instructions set in effect is: $instructionsInEffect\n";	
      }
      $regressorLine=processInstructionsLogEntry($instructionsInEffect, $logEntry, $phase, $facesRecalledHash_ref, $recallSummaryHash_ref);
      if ($DEBUG) {
	print "### processLogEntry - Found an instructions block\n";
	print "### processLogEntry - Regressor line: $regressorLine\n";
      }
    } elsif ( @{$logEntry}[$i] =~ /\s*Procedure:\sEpochR?/ ) {
      if ($DEBUG) {
	print "### processLogEntry - EPOCH detected\n";
      }
      
      $regressorLine=processEpochLogEntry($recallSummaryHash_ref, $facesRecalledHash_ref, $instructionsInEffect, $logEntry, $phase);
      if ($DEBUG) {
	print "### processLogEntry - Regressor line: $regressorLine\n";
      }
    }
    # else {
    #   print "### processLogEntry -  I have no idea what to do with this log entry. Oops :-)\n";
    # }
  }
  
  if ($DEBUG) {
    print "### processLogEntry - firstStimulusPresentedAt is $firstStimulusPresentedAt\n";
  }

  return ($instructionsInEffect, $regressorLine);
}

sub processInstructionsLogEntry($$$$$) {
  if ($DEBUG) {
    print "### ENTRY processInstructionsLogEntry\n";
  }
  
  my $instructionsInEffect=shift;
  my $logEntry=shift;
  my $phase=shift;
  my $facesRecalledHash_ref=shift;
  my $recallSummaryHash_ref=shift;
  
  my $onset=$ILLEGAL_VALUE;

  if ($DEBUG) {
    print "### processInstructionsLogEntry - Instructions set in effect is: $instructionsInEffect\n";	
  }
  for my $i ( 0 .. @{$logEntry} - 1 ) {
    if ($DEBUG) {
      print "### processInstructionsLogEntry - Log line: " . @{$logEntry}[$i] . "\n";
    }
    
    if ( @{$logEntry}[$i] =~ /\s*(Instructions([ABCD]).OnsetTime):\s*(\d*)/ ) {
      my $instructionsBlock=$2;
      $onset=$3;
      ## print "### instructionsBlock=$instructionsBlock, onset=$onset\n";

      push @{ $globalRecallInstructionsOnsetTimes{$instructionsBlock} }, $onset;
      
      if ($DEBUG) {
	print "### processInstructionsLogEntry - Found the onset of the first instruction block\n";
	print "### processInstructionsLogEntry - The first stimulus was presented at: $onset\n"; 
	print "### processInstructionsLogEntry - firstStimulusPresentedAt is $firstStimulusPresentedAt\n";
      }
      if ($firstStimulusPresentedAt == $ILLEGAL_VALUE ) {
       	$firstStimulusPresentedAt=$3;
	if ($DEBUG) {
	  print "### processInstructionsLogEntry - Set firstStimulusPresentedAt to $firstStimulusPresentedAt\n";
	}
      }
    }
  }

  my $offset=$ILLEGAL_OFFSET_VALUE;
  my $line=makeStimulusRegressorLine($instructionsInEffect, "", $INSTRUCTIONS, $onset, $offset, $phase, $facesRecalledHash_ref, $recallSummaryHash_ref);

  if ($DEBUG) {
    print "### processInstructionsLogEntry - Regressor line is: $line\n";
    print "### EXIT processInstructionsLogEntry\n";
  }

  return ($line);
}


sub processEpochLogEntry($$$$$) {
  if ($DEBUG) {
    print "### ENTRY processEpochLogEntry\n";
  }

  my $recallSummaryHash_ref=shift;
  my $facesRecalledHash_ref=shift;
  my $instructionsInEffect=shift;
  my $logEntry=shift;
  my $phase=shift;

  my $targetEmotion="";
  my $targetId=$ILLEGAL_VALUE;
  my $run=$ILLEGAL_VALUE;
  my $itiDur=$ILLEGAL_VALUE;
  my $onset=$ILLEGAL_VALUE;
  my $correctAnswer=$ILLEGAL_VALUE;
  my $response=$ILLEGAL_VALUE;

  if ($DEBUG) {
    print "### processEpochLogEntry - Instruction set in effect is $instructionsInEffect\n";
  }

  #if ($DEBUG) {
  #  print "### processEpochLogEntry - LogEntry is: @{$logEntry}\n";
  #}

  ## $i iterates over the lines in a log frame
  for my $i ( 0 .. @{$logEntry} - 1 ) {
    if ($DEBUG) {
      print "### processEpochLogEntry - Log line: " . @{$logEntry}[$i] . "\n";
    }

    ########## Target ##########
    if ( @{$logEntry}[$i] =~ /\s*(Target:\s(\w)(\d+))|(fix\.bmp)/ ) {
      if (! defined ($4) ) {
	$targetEmotion=lc("$2");
	$targetId=$3;
	if ($phase == $ACQUISITION) {
	  push @{ $globalEmotionToActorMap{$targetEmotion} }, $targetId;
	}
	if ($DEBUG) {
	  print "### processEpochLogEntry - Target for this epoch is: $2$3\n";
	}
      } else { 
	## short for fixation
	$targetEmotion="fix";
	if ($DEBUG) {
	  print "### processEpochLogEntry - Target for this epoch is: $4\n";
	}
      }
    }

    ########## Run number  ##########
    if ( @{$logEntry}[$i] =~ /\s*Running:\s(Run(\d))/ ) {
      $run=$2;
      if ($DEBUG) {
	print "### processEpochLogEntry - This is run: " . $2 . "\n";
      }
    }

    ########## ITIDur  ##########
    if ( @{$logEntry}[$i] =~ /\s*ITIDur:\s(\d*)/ ) {
      $itiDur=$1;
      if ($DEBUG) {
	print "### processEpochLogEntry - ITIDur is: " . $1 . "\n";
      }
    }

    ########## ShowFace.OnsetTime  ##########
    if ( @{$logEntry}[$i] =~ /\s*ShowFace\.OnsetTime:\s(\d*)/ ) {
      $onset=$1;
      if ($DEBUG) {
	print "### processEpochLogEntry - OnsetTime is: " . $1 . "\n";
      }
    }
    
    if ($phase == $RECALL) {
      if ( @{$logEntry}[$i] =~ /\s*CorrectAnswer:\s(\d*)/ ) {
	$correctAnswer=$1;
	if ($DEBUG) {
	  print "### processEpochLogEntry - Correct Answer is: " . $1 . "\n";
	}
      }
      
      if ( @{$logEntry}[$i] =~ /\s*ShowFace\.RESP:\s(\d*)/ ) {
	$response=$1;
	if ($DEBUG) {
	  print "### processEpochLogEntry - Response is: " . $1 . "\n";
	}
      }
    }

  } ## end of for my $i ( 0 .. @{$logEntry} - 1 ) {

  if ($phase == $RECALL) {

    if ($correctAnswer eq "1" && $response eq "1" ) {
      ## Hit
      if ($DEBUG) {
	print "### processEpochLogEntry - Image was correctly recalled\n";
      }

      $recallSummaryHash_ref->{ $TOTAL_RECALL_HITS } = $recallSummaryHash_ref->{ $TOTAL_RECALL_HITS } + 1;
      $facesRecalledHash_ref->{ $targetId } = $REMEMBERED;

    } elsif ($correctAnswer eq "1" && $response eq "2" ) {
      ## Miss
      if ($DEBUG) {
	print "### processEpochLogEntry - Image was missed\n";
      }
      
      $recallSummaryHash_ref->{ $TOTAL_RECALL_MISSES } = $recallSummaryHash_ref->{ $TOTAL_RECALL_MISSES } + 1;
      $facesRecalledHash_ref->{ $targetId } = $NOT_REMEMBERED;
    } elsif ($correctAnswer eq "2" && $response eq "2" ) {
      ## correct rejection
      if ($DEBUG) {
	print "### processEpochLogEntry - Image was correctly rejected\n";
      }
      
      $recallSummaryHash_ref->{ $TOTAL_RECALL_CORRECT_REJECTIONS } = $recallSummaryHash_ref->{ $TOTAL_RECALL_CORRECT_REJECTIONS } + 1;
      $facesRecalledHash_ref->{ $targetId } = $NOT_REMEMBERED;

    } elsif($correctAnswer eq "2" && $response eq "1" ) {
      ## false alarm
      if ($DEBUG) {
	print "### processEpochLogEntry - Image was incorrectly recalled\n";
      }
      
      $recallSummaryHash_ref->{ $TOTAL_RECALL_FALSE_ALARMS } = $recallSummaryHash_ref->{ $TOTAL_RECALL_FALSE_ALARMS } + 1;
      $facesRecalledHash_ref->{ $targetId } = $FALSE_ALARM;
    } else {
      if ($DEBUG) {
	print "### processEpochLogEntry - No reponse\n";
      }

      # print "### processEpochLogEntry - I don't know what to do in this situation, assuming it's an ommission: correctAnswer: $correctAnswer response: $response\n";

      $recallSummaryHash_ref->{ $TOTAL_RECALL_OMISSIONS } = $recallSummaryHash_ref->{ $TOTAL_RECALL_OMISSIONS } + 1;
      $facesRecalledHash_ref->{ $targetId } = $OMISSION;
    }

    # for my $key ( sort keys %{ $facesRecalledHash_ref} ) {
    #   my $value = $facesRecalledHash_ref->{$key};
    #   print "$key => $value\n";
    # }
    
    ## save the target Id and emotion (which will always be neutral)
    ## from the recall phase to the global trial array
    my @ar=( $targetId, $targetEmotion );
    push @globalRecallTrialsArray, [ @ar ];

  } ## end of if ($phase == $RECALL)

  my $offset=$ILLEGAL_OFFSET_VALUE;

  my $line="### ERROR if you can see this something went wrong. Do not use these regressors.###";
  if ( $targetEmotion eq "f" ) {
    $line= makeStimulusRegressorLine($instructionsInEffect, $targetId, $FEARFUL, $onset, $offset, $phase, $facesRecalledHash_ref, $recallSummaryHash_ref);
    if ($DEBUG) {
      print "### processEpochLogEntry - FEARFUL face detected\n";
    }
  } elsif ( $targetEmotion eq "h" ) {
    $line = makeStimulusRegressorLine($instructionsInEffect, $targetId, $HAPPY, $onset, $offset, $phase, $facesRecalledHash_ref, $recallSummaryHash_ref);
    if ($DEBUG) {
      print "### processEpochLogEntry - HAPPY face detected\n";
    }
  } elsif ( $targetEmotion eq "s" ) {
    $line = makeStimulusRegressorLine($instructionsInEffect, $targetId, $SAD, $onset, $offset, $phase, $facesRecalledHash_ref, $recallSummaryHash_ref);
    if ($DEBUG) {
      print "### processEpochLogEntry - SAD face detected\n";
    }
  } elsif ( $targetEmotion eq "n" ) {
    $line = makeStimulusRegressorLine($instructionsInEffect, $targetId, $NEUTRAL, $onset, $offset, $phase, $facesRecalledHash_ref, $recallSummaryHash_ref);
    if ($DEBUG) {
      print "### processEpochLogEntry - NEUTRAL face detected\n";
    }
  } elsif (  $targetEmotion eq "fix" ) {
    $line = makeStimulusRegressorLine($instructionsInEffect, $targetId, $FIXATION, $onset, $offset, $phase, $facesRecalledHash_ref, $recallSummaryHash_ref);
    if ($DEBUG) {
      print "### processEpochLogEntry - FIXATION cross detected\n";
    }
  }
  
  if ($DEBUG) {
    print "### processEpochLogEntry - line is: $line\n";
  }

  if ($DEBUG) {
    print "### EXIT processEpochLogEntry\n";
  }

  return ($line);
}

sub printInstructionsInEffectLine($) {
  my $regressorFilenamePrefix=shift;
  if ($DEBUG) {
    print "### ENTRY printInstructionsInEffectLine\n";
  }
  
  #open(my $fh, ">", $regressorFilename) or die "Cannot open $regressorFilename: $!\n";

  my @onsetOffsetMat=();
  my @onsetDurationMat=();
  my ( $startA, $endA, $startB, $endB, $startC, $endC, $startD, $endD, $durationA, $durationB, $durationC, $durationD ); 
  # print $fh "#InstrA\tInstrB\tInstrC\tInstrD\n";
  ## $i counts the epochs of blocks in the experiment
  for ( my $i=0; $i<4; $i=$i+1 ) {
    #normalizeRelativeToFirstStimulusPresentationTime()
    ## convert to seconds by dividing by 1000
    $startA=normalizeRelativeToFirstStimulusPresentationTime      ( @{ $globalRecallInstructionsOnsetTimes{'A'} }[$i] )/1000;
    $endA=$startB=normalizeRelativeToFirstStimulusPresentationTime( @{ $globalRecallInstructionsOnsetTimes{'B'} }[$i] )/1000;
    $endB=$startC=normalizeRelativeToFirstStimulusPresentationTime( @{ $globalRecallInstructionsOnsetTimes{'C'} }[$i] )/1000;
    $endC=$startD=normalizeRelativeToFirstStimulusPresentationTime( @{ $globalRecallInstructionsOnsetTimes{'D'} }[$i] )/1000;
    if ( $i == 3 ) {
      $endD="859";
    } else {
      $endD=normalizeRelativeToFirstStimulusPresentationTime      ( @{ $globalRecallInstructionsOnsetTimes{'A'} }[$i+1] )/1000
    }
    
    $durationA=$endA-$startA;
    $durationB=$endB-$startB;
    $durationC=$endC-$startC;
    #if ( $i == 3 ) {
    #  $durationD=9999;
    #} else { 
    $durationD=$endD-$startD;
    #}
    
    #code to trasnpose matrix
    $onsetOffsetMat[0][$i] = "$startA:$endA";
    $onsetOffsetMat[1][$i] = "$startB:$endB";
    $onsetOffsetMat[2][$i] = "$startC:$endC";
    $onsetOffsetMat[3][$i] = "$startD:$endD";

    #code to trasnpose matrix
    $onsetDurationMat[0][$i] = "$startA:$durationA";
    $onsetDurationMat[1][$i] = "$startB:$durationB";
    $onsetDurationMat[2][$i] = "$startC:$durationC";
    $onsetDurationMat[3][$i] = "$startD:$durationD";

    # if ( $regressorType == $ONSET_OFFSET_REGRESSOR ) {
    # print $fh "$startA:$endA\t$startB:$endB\t$startC:$endC\t$startD:$endD\n";
    # } elsif ( $regressorType == $ONSET_DURATION_REGRESSOR )  {
    #   print $fh "$startA:$durationA\t$startB:$durationB\t$startC:$durationC\t$startD:$durationD\n";
    # } else {
    #   die "printInstructionsInEffectLine - Cant make a regressor of that type\n";
    # }
  }
  my $regressorPrefix;
  if ( $regressorType      == $ONSET_OFFSET_REGRESSOR ) {
    $regressorPrefix="onsetOffset";
  } elsif ( $regressorType == $ONSET_ONLY_REGRESSOR )  {
    $regressorPrefix="onsetOnly";
  } elsif ( $regressorType == $ONSET_DURATION_REGRESSOR )  {
    $regressorPrefix="onsetDuration";
  }

  ## code to print trasnposed matrix  
  ## $i iterates over the instruction blocks
  my $fh;
  for ( my $i=0; $i<4; $i++ ) {
    if ($i == 0 ) {
      open($fh, ">", "$regressorFilenamePrefix.$regressorPrefix.instrA.tab") or die "Cannot open $regressorFilenamePrefix.$regressorPrefix.instrA.tab: $!\n";
      print "Writing instructions in effect regressors to $regressorFilenamePrefix.$regressorPrefix.instrA.tab\n";
    } elsif ($i == 1 ) {
      open($fh, ">", "$regressorFilenamePrefix.$regressorPrefix.instrB.tab") or die "Cannot open $regressorFilenamePrefix.$regressorPrefix.instrB.tab: $!\n";
      print "Writing instructions in effect regressors to $regressorFilenamePrefix.$regressorPrefix.instrB.tab\n";
    } elsif ($i == 2 ) {
      open($fh, ">", "$regressorFilenamePrefix.$regressorPrefix.instrC.tab") or die "Cannot open $regressorFilenamePrefix.$regressorPrefix.instrC.tab: $!\n";
      print "Writing instructions in effect regressors to $regressorFilenamePrefix.$regressorPrefix.instrC.tab\n";
    } else {
      open($fh, ">", "$regressorFilenamePrefix.$regressorPrefix.instrD.tab") or die "Cannot open $regressorFilenamePrefix.$regressorPrefix.instrD.tab: $!\n";
      print "Writing instructions in effect regressors to $regressorFilenamePrefix.$regressorPrefix.instrD.tab\n";
    }

    ## $j iterates over the repetitions of the instruction blocks
    for (my $j=0; $j<4; $j++ ) { 
      if ( $regressorType == $ONSET_OFFSET_REGRESSOR ) {
	print $fh "$onsetOffsetMat[$i][$j] ";
      } elsif ( $regressorType == $ONSET_DURATION_REGRESSOR )  {
	print $fh "$onsetDurationMat[$i][$j] ";
      } elsif ( $regressorType == $ONSET_ONLY_REGRESSOR )  {
	print $fh " ";
      } else {
	die "printInstructionsInEffectLine - Can't make a regressor of that type\n";
      }
    }
    print $fh "\n";
    close $fh;
  }

  ## close $fh;

  if ($DEBUG) {
    print "### EXIT printInstructionsInEffectLine\n";
  }
}


sub printRegressorLines($$) {
  if ($DEBUG) {
    print "### ENTRY printRegressorLines\n";
  }

  my $regressorLinesRef=shift;
  my $regressorFilename=shift;

  my @regressorLines=@{$regressorLinesRef};
  
  my $header="#InstrA\tInstrB\tInstrC\tInstrD\tInstR\tFearful\tFearfulRemembered\tFearfulNotRemembered\tFearfulOmissions\tFixation\tHappy\tHappyRemembered\tHappyNotRemembered\tHappyOmissions\tInstructions\tNeutral\tNeutralRemembered\tNeutralNotRemembered\tNeutralOmissions\tSad\tSadRemembered\tSadNotRemembered\tSadOmissions";

  open(my $fh, ">", $regressorFilename) or die "Cannot open $regressorFilename: $!\n";

  if ($DEBUG) {
    print "$header\n";
  }
  print $fh "$header\n";

  for my $i ( 0 .. $#regressorLines ) {
    if ($DEBUG) {
      print "Regressor line $i -> " . $regressorLines[$i] . "\n";
    }
    #print $regressorLines[$i] . "\n";
    print $fh $regressorLines[$i] . "\n";

  }
  
  close ($fh);
  print "Wrote " . ( $#regressorLines + 1 ) . " lines to the regressor file.\n";

  if ($DEBUG) {
    print "### EXIT printRegressorLines\n";
  }
}


sub printSummaryLine($$$) {
  if ($DEBUG) {
    print "### ENTRY printSummaryLine\n";
  }
  
  my $subjectId=shift;
  my $recallSummaryHash_ref=shift;
  my $summaryFilename=shift;

  my $header=join ( "\t", ("Subject", 
			   ${TOTAL_RECALL_HITS}, 
			   ${FEARFUL_RECALL_HITS}, 
			   ${HAPPY_RECALL_HITS}, 
			   ${NEUTRAL_RECALL_HITS}, 
			   ${SAD_RECALL_HITS}, 
			   # ${TOTAL_ACQ_FEARFUL_RECALL_HITS}, 
			   # ${TOTAL_ACQ_HAPPY_RECALL_HITS}, 
			   # ${TOTAL_ACQ_NEUTRAL_RECALL_HITS}, 
			   # ${TOTAL_ACQ_SAD_RECALL_HITS}, 
			   # ${TOTAL_ACQ_FEARFUL_RECALL_MISSES}, 
			   # ${TOTAL_ACQ_HAPPY_RECALL_MISSES}, 
			   # ${TOTAL_ACQ_NEUTRAL_RECALL_MISSES}, 
			   # ${TOTAL_ACQ_SAD_RECALL_MISSES}, 
			   ${TOTAL_RECALL_MISSES}, 
			   ${TOTAL_RECALL_CORRECT_REJECTIONS}, 
			   ${TOTAL_RECALL_FALSE_ALARMS}, 
			   ${TOTAL_RECALL_OMISSIONS},
			   "Total Trials",
			   ${TOTAL_RECALL_HITS} . "Proportion", 
			   ${FEARFUL_RECALL_HITS} . "Proportion", 
			   ${HAPPY_RECALL_HITS} . "Proportion", 
			   ${NEUTRAL_RECALL_HITS} . "Proportion", 
			   ${SAD_RECALL_HITS} . "Proportion", 
			   # ${TOTAL_ACQ_FEARFUL_RECALL_HITS}, 
			   # ${TOTAL_ACQ_HAPPY_RECALL_HITS}, 
			   # ${TOTAL_ACQ_NEUTRAL_RECALL_HITS}, 
			   # ${TOTAL_ACQ_SAD_RECALL_HITS}, 
			   # ${TOTAL_ACQ_FEARFUL_RECALL_MISSES}, 
			   # ${TOTAL_ACQ_HAPPY_RECALL_MISSES}, 
			   # ${TOTAL_ACQ_NEUTRAL_RECALL_MISSES}, 
			   # ${TOTAL_ACQ_SAD_RECALL_MISSES}, 
			   ${TOTAL_RECALL_MISSES} . "Proportion", 
			   ${TOTAL_RECALL_CORRECT_REJECTIONS} . "Proportion", 
			   ${TOTAL_RECALL_FALSE_ALARMS} . "Proportion", 
			   ${TOTAL_RECALL_OMISSIONS} . "Proportion"
			  )
		  );

  my $fh;
  if ( -f $summaryFilename) {
    open($fh, ">>", $summaryFilename) or die "Cannot open $summaryFilename: $!\n";
  } else {
    open($fh, ">", $summaryFilename) or die "Cannot open $summaryFilename: $!\n";
    print $fh "$header\n";
  }

  if ($DEBUG) {
    print "$header\n";
  }


  ## there are 32 total actors (one per trial) depicting one of 4
  ## face-emotion types
  my $totalTrials=32; 

  # my $totalTrials=0 + 
  #   $recallSummaryHash_ref->{ $TOTAL_RECALL_HITS } + 
  #     $recallSummaryHash_ref->{ $TOTAL_RECALL_MISSES } + 
  # 	$recallSummaryHash_ref->{ $TOTAL_RECALL_CORRECT_REJECTIONS } + 
  # 	  $recallSummaryHash_ref->{ $TOTAL_RECALL_FALSE_ALARMS } + 
  # 	    $recallSummaryHash_ref->{ $TOTAL_RECALL_OMISSIONS };

  ## 32 face-emotions divided by 4 emotion types (fearful, happy,
  ## neutral, and sad). Thus we use 32 / 4 = 6 in the present
  ## calculation.
  ##
  ## Note that this is different from the Pine et al 2004 paper where
  ## they only presented 3 emotion faces and therefroe had 24 / 3 = 8
  ## trials per face-emotion

  my $maxTrialsPerFaceEmotion = $totalTrials / 4;

  # my $totalTrials=0 + 
  #   $recallSummaryHash_ref->{ $TOTAL_RECALL_HITS } + 
  #     $recallSummaryHash_ref->{ $TOTAL_RECALL_MISSES } + 
  # 	$recallSummaryHash_ref->{ $TOTAL_RECALL_CORRECT_REJECTIONS } + 
  # 	  $recallSummaryHash_ref->{ $TOTAL_RECALL_FALSE_ALARMS } + 
  # 	    $recallSummaryHash_ref->{ $TOTAL_RECALL_OMISSIONS };
  
  my $summaryLine=join ("\t", ($subjectId,
			       $recallSummaryHash_ref->{ $TOTAL_RECALL_HITS },
			       $recallSummaryHash_ref->{ $FEARFUL_RECALL_HITS },
			       $recallSummaryHash_ref->{ $HAPPY_RECALL_HITS },
			       $recallSummaryHash_ref->{ $NEUTRAL_RECALL_HITS },
			       $recallSummaryHash_ref->{ $SAD_RECALL_HITS },
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_FEARFUL_RECALL_HITS } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_HAPPY_RECALL_HITS } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_NEUTRAL_RECALL_HITS } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_SAD_RECALL_HITS } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_FEARFUL_RECALL_MISSES } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_HAPPY_RECALL_MISSES } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_NEUTRAL_RECALL_MISSES } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_SAD_RECALL_MISSES } || 0,
			       $recallSummaryHash_ref->{ $TOTAL_RECALL_MISSES },
			       $recallSummaryHash_ref->{ $TOTAL_RECALL_CORRECT_REJECTIONS }, 
			       $recallSummaryHash_ref->{ $TOTAL_RECALL_FALSE_ALARMS },
			       $recallSummaryHash_ref->{ $TOTAL_RECALL_OMISSIONS },
			       $totalTrials,
			       $recallSummaryHash_ref->{ $TOTAL_RECALL_HITS } / $totalTrials,
			       $recallSummaryHash_ref->{ $FEARFUL_RECALL_HITS } / $maxTrialsPerFaceEmotion,
			       $recallSummaryHash_ref->{ $HAPPY_RECALL_HITS } / $maxTrialsPerFaceEmotion,
			       $recallSummaryHash_ref->{ $NEUTRAL_RECALL_HITS } / $maxTrialsPerFaceEmotion,
			       $recallSummaryHash_ref->{ $SAD_RECALL_HITS } / $maxTrialsPerFaceEmotion,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_FEARFUL_RECALL_HITS } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_HAPPY_RECALL_HITS } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_NEUTRAL_RECALL_HITS } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_SAD_RECALL_HITS } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_FEARFUL_RECALL_MISSES } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_HAPPY_RECALL_MISSES } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_NEUTRAL_RECALL_MISSES } || 0,
			       # $recallSummaryHash_ref->{ $TOTAL_ACQ_SAD_RECALL_MISSES } || 0,
			       $recallSummaryHash_ref->{ $TOTAL_RECALL_MISSES } / $totalTrials,
			       $recallSummaryHash_ref->{ $TOTAL_RECALL_CORRECT_REJECTIONS } / $totalTrials, 
			       $recallSummaryHash_ref->{ $TOTAL_RECALL_FALSE_ALARMS } / $totalTrials,
			       $recallSummaryHash_ref->{ $TOTAL_RECALL_OMISSIONS } / $totalTrials
			      )
		       );
  if ($DEBUG) {
    print $summaryLine . "\n";
  }
  print $fh $summaryLine . "\n";
  
  close ($fh);
  
  if ($DEBUG) {
    print "### EXIT printSummaryLine\n";
  }
}

sub normalizeRelativeToFirstStimulusPresentationTime($) {
  my $time=shift;

  if ($firstStimulusPresentedAt == $ILLEGAL_VALUE ) {
    die ("### normalizeRelativeToFirstStimulusPresentationTime - firstStimulusPresentedAt is $ILLEGAL_VALUE. Cannot normalise time ($time) to this base time. Exiting\n");
  }
      
  return ($time - $firstStimulusPresentedAt);
}

sub makeInstructionsInEffectRegressorComponent($) {
  if ($DEBUG) {
    print "### ENTRY makeInstructionsInEffectRegressorComponent\n";
  }

  my $instructionsInEffect=shift;
  my $regressorComponent="";

  if ($instructionsInEffect eq "A") {
    $regressorComponent="1\t0\t0\t0\t0"
  } elsif ($instructionsInEffect eq "B") {
    $regressorComponent="0\t1\t0\t0\t0"
  } elsif ($instructionsInEffect eq "C") {
    $regressorComponent="0\t0\t1\t0\t0"
  } elsif ($instructionsInEffect eq "D") {
    $regressorComponent="0\t0\t0\t1\t0"
  } elsif ($instructionsInEffect eq "R") {
    $regressorComponent="0\t0\t0\t0\t1"
  } else {
    $regressorComponent="U\tU\tU\tU\tU"
  }
  
  if ($DEBUG) {
    print "### EXIT makeInstructionsInEffectRegressorComponent\n";
  }

  return ($regressorComponent);
}

sub makeStimulusRegressorLine($$$$$$$$) {
  if ($DEBUG) {
    print "### ENTRY makeStimulusRegressorLine\n";
  }
  
  my $instructionsInEffect=shift;
  my $targetId=shift;
  my $stimulusType=shift;
  my $onset=shift;
  my $offset=shift;
  my $phase=shift;
  my $facesRecalledHash_ref=shift;
  my $recallSummaryHash_ref=shift;

  if ($DEBUG) {
    print join ( "\n", (
			"### makeStimulusRegressorLine - instructionsInEffect = $instructionsInEffect",
			"### makeStimulusRegressorLine - targetId = $targetId",
			"### makeStimulusRegressorLine - stimulusType = $stimulusType", 
			"### makeStimulusRegressorLine - onset = $onset",
			"### makeStimulusRegressorLine - offset = $offset", 
			"### makeStimulusRegressorLine - phase = $phase"
		       )
	       );
    print "\n";
  }


  ## check the length of the keys in $facesRecalledHash_ref if it's
  ## zero we have a problem and should quit

  my $duration;

  if ($phase == $ACQUISITION) {
    $onset=normalizeRelativeToFirstStimulusPresentationTime($onset);
    $offset=normalizeRelativeToFirstStimulusPresentationTime($offset);

    if ( defined($forceDuration) ) {
      ## $duration = $forceDuration;
      if ($stimulusType==$INSTRUCTIONS) {
	$duration=$INSTRUCTIONS_DURATION;
      } else {
	$duration=$TRIAL_DURATION;
      }
    } else {
      $duration = $offset - $onset;
    }

    # convert to seconds from milliseconds
    $onset=$onset/1000;
    $offset=$offset/1000;
    $duration=$duration/1000;
    
  } elsif ($phase == $RECALL) {
    $onset=9999;
    $offset=9999;
    $duration=9999;
  }

  ## a this point check whether onset is less than zero, if this is
  ## the case we should die immediatly as under all circumstances this
  ## is be a major problem
  if ( $onset < 0 ) {
    die "### makeStimulusRegressorLine - Onset ($onset) is less than zero. Exiting\n";
  }

  ## a this point check whether offset and duration are less than
  ## zero, if this is the case we should die immediatly but only under
  ## certain circumstances as this is not always a major problem
  if ( $offset < 0 && $regressorType == $ONSET_OFFSET_REGRESSOR ) {
    die "### makeStimulusRegressorLine - Offset ($offset) is less than zero. Exiting\n";
  }
  if ( $duration < 0 && $forceDuration == 0) {
    die "### makeStimulusRegressorLine - Duration ($duration) is less than zero. Exiting\n";
  }

  my $line="**** makeStimulusRegressorLine - ERROR default value assigned in makeStimulusRegressorLine. If you can see this in the text file some thing really bad has gone wrong. **** ";
  
  ## this line is useful to test that the sprintf statements are
  ## putting the times in the correct position for each stimulus type
  ## $facesRecalledHash_ref->{ $targetId } = $OMISSION;

  if ($regressorType == $ONSET_DURATION_REGRESSOR ) {
    if ($stimulusType==$FEARFUL) {
      if ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $REMEMBERED ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - FEARFUL REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}
		
	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_FEARFUL_RECALL_HITS } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_FEARFUL_RECALL_HITS } || 0 ) + 1;
	# }

	$line = sprintf("%s\t%0.3f:%0.3f\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			$onset, $duration, # fearful
			$onset, $duration, # fearful remembered
			9999, 9999,        # fearful not remembered
			9999, 9999,        # fearful omission

			9999, 9999, # fixation
			
			9999, 9999, # happy 
			9999, 9999, # happy remembered
			9999, 9999, # happy not remembered
			9999, 9999, # happy omission
			
			9999, 9999, # instructions
			
			9999, 9999, # neutral
			9999, 9999, # neutral remembered 
			9999, 9999, # neutral not remembered 
			9999, 9999, # neutral omission
			
			9999, 9999, # sad
			9999, 9999, # sad remembered 
			9999, 9999, # sad not remembered 
			9999, 9999  # sad omission 
		       );

      } elsif ( defined $facesRecalledHash_ref->{ $targetId } &&  $facesRecalledHash_ref->{ $targetId } eq $NOT_REMEMBERED  ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - FEARFUL NOT REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}

	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_FEARFUL_RECALL_MISSES } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_FEARFUL_RECALL_MISSES } || 0 ) + 1;
	# }

	$line = sprintf("%s\t%0.3f:%0.3f\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			$onset, $duration, # fearful
			9999, 9999,        # fearful remembered
			$onset, $duration, # fearful not remembered
			9999, 9999,        # fearful omission
			
			9999, 9999, # fixation
			
			9999, 9999, # happy 			
			9999, 9999, # happy remembered
			9999, 9999, # happy not remembered
			9999, 9999, # happy omission
			
			9999, 9999, # instructions
			
			9999, 9999, # neutral
			9999, 9999, # neutral remembered 
			9999, 9999, # neutral not remembered 
			9999, 9999, # neutral omission
			
			9999, 9999, # sad 
			9999, 9999, # sad remembered 
			9999, 9999, # sad not remembered 
			9999, 9999  # sad omission 
		       );
	
      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $OMISSION ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - FEARFUL OMISSION: Target ($targetId) is in facesRecalledHash\n";
	}
	
	$line = sprintf("%s\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			$onset, $duration, # fearful
			9999, 9999,        # fearful remembered
			9999, 9999,        # fearful not remembered
			$onset, $duration, # fearful omission

			9999, 9999, # fixation			

			9999, 9999, # happy
			9999, 9999, # happy remembered
			9999, 9999, # happy not remembered
			9999, 9999, # happy omission
			
			9999, 9999, # instructions
			
			9999, 9999, # neutral
			9999, 9999, # neutral remembered 
			9999, 9999, # neutral not remembered 
			9999, 9999, # neutral omission
			
			9999, 9999, # sad 
			9999, 9999, # sad remembered 
			9999, 9999, # sad not remembered 
			9999, 9999  # sad omission 
		       );
	
      } else {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - FEARFUL: No entry for target ($targetId) or it's a false alarm in facesRecalledHash\n";
	}

	$line = sprintf("%s\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			$onset, $duration, # fearful
			9999, 9999,        # fearful remembered
			9999, 9999,        # fearful not remembered
			9999, 9999,        # fearful omission
			
			9999, 9999, # fixation
			
			9999, 9999, # happy 
			9999, 9999, # happy remembered
			9999, 9999, # happy not remembered
			9999, 9999, # happy omission
			
			9999, 9999, # instructions
			
			9999, 9999, # neutral
			9999, 9999, # neutral remembered 
			9999, 9999, # neutral not remembered 
			9999, 9999, # neutral omission
			
			9999, 9999, # sad
			9999, 9999, # sad remembered 
			9999, 9999, # sad not remembered 
			9999, 9999  # sad omission 
		       );
      }
    } elsif ($stimulusType==$FIXATION) {
      $line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
		      makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
		      
		      9999, 9999, # fearful
		      9999, 9999, # fearful remembered
		      9999, 9999, # fearful not remembered
		      9999, 9999, # fearful omission
		      
		      $onset, $duration, # fixation
		      
		      9999, 9999, # happy 
		      9999, 9999, # happy remembered
		      9999, 9999, # happy not remembered
		      9999, 9999, # happy omission
		      
		      9999, 9999, # instructions
		      
		      9999, 9999, # neutral  
		      9999, 9999, # neutral remembered 
		      9999, 9999, # neutral not remembered 
		      9999, 9999, # neutral omission
		      
		      9999, 9999, # sad 
		      9999, 9999, # sad remembered 
		      9999, 9999, # sad not remembered 
		      9999, 9999  # sad omission 
		     );
      
    } elsif ($stimulusType==$HAPPY) {
      if ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $REMEMBERED  ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - HAPPY REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}
	
	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_HAPPY_RECALL_HITS } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_HAPPY_RECALL_HITS } || 0 ) + 1;
	# }
	
	$line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			9999, 9999, # fearful remembered
			9999, 9999, # fearful
			9999, 9999, # fearful not remembered
			9999, 9999, # fearful omission
			
			9999, 9999, # fixation
			
			$onset, $duration, # happy
			$onset, $duration, # happy remembered
			9999, 9999,        # happy not remembered
			9999, 9999,        # happy omission
			
			9999, 9999, # instructions
			
			9999, 9999, # neutral
			9999, 9999, # neutral remembered 
			9999, 9999, # neutral not remembered 
			9999, 9999, # neutral omission
			
			9999, 9999, # sad
			9999, 9999, # sad remembered 
			9999, 9999, # sad not remembered 
			9999, 9999  # sad omission 
		       );
	
      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $NOT_REMEMBERED  ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - HAPPY NOT REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}
	
	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_HAPPY_RECALL_MISSES } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_HAPPY_RECALL_MISSES } || 0 ) + 1;
	# }
	
	$line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			9999, 9999, # fearful
			9999, 9999, # fearful remembered
			9999, 9999, # fearful not remembered
			9999, 9999, # fearful omission
			
			9999, 9999, # fixation
			
			$onset, $duration, # happy
			9999, 9999,        # happy remembered
			$onset, $duration, # happy not remembered
			9999, 9999,        # happy omission
			
			9999, 9999, # instructions
			
			9999, 9999, # neutral
			9999, 9999, # neutral remembered 
			9999, 9999, # neutral not remembered 
			9999, 9999, # neutral omission
			
			9999, 9999, # sad
			9999, 9999, # sad remembered 
			9999, 9999, # sad not remembered 
			9999, 9999  # sad omission 
		       );
	
      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $OMISSION ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - HAPPY OMISSION: Target ($targetId) is in facesRecalledHash\n";
	}
	$line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			9999, 9999, # fearful
			9999, 9999, # fearful remembered			
			9999, 9999, # fearful not remembered
			9999, 9999, # fearful omission
			
			9999, 9999, # fixation
			
			$onset, $duration, # happy
			9999, 9999,        # happy remembered
			9999, 9999,        # happy not remembered
			$onset, $duration, # happy omission
			
			9999, 9999, # instructions
			
			9999, 9999, # neutral
			9999, 9999, # neutral remembered 
			9999, 9999, # neutral not remembered 
			9999, 9999, # neutral omission
			
			9999, 9999, # sad
			9999, 9999, # sad remembered 
			9999, 9999, # sad not remembered 
			9999, 9999  # sad omission 
		       );
	
      } else {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - HAPPY: No entry for target ($targetId) or it's a false alarm in facesRecalledHash\n";
	}
	
	$line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			9999, 9999, # fearful
			9999, 9999, # fearful remembered
			9999, 9999, # fearful not remembered
			9999, 9999, # fearful omission
			
			9999, 9999, # fixation
			
			$onset, $duration, # happy 
			9999, 9999,	   # happy remembered
			9999, 9999,	   # happy not remembered
			9999, 9999,	   # happy omission
			
			9999, 9999, # instructions
			
			9999, 9999, # neutral
			9999, 9999, # neutral remembered 
			9999, 9999, # neutral not remembered 
			9999, 9999, # neutral omission
			
			9999, 9999, # sad
			9999, 9999, # sad remembered 
			9999, 9999, # sad not remembered 
			9999, 9999  # sad omission 
		       );
      }
    } elsif ($stimulusType==$INSTRUCTIONS) {
      $line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
		      makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
		      
		      9999, 9999,	# fearful
		      9999, 9999,	# fearful remembered
		      9999, 9999,	# fearful not remembered
		      9999, 9999,	# fearful omission
		      
		      9999, 9999, # fixation
		      
		      9999, 9999,	# happy 
		      9999, 9999,	# happy remembered
		      9999, 9999,	# happy not remembered
		      9999, 9999,	# happy omission

		      $onset, $duration, # instructions		      

		      9999, 9999,	# neutral 
		      9999, 9999,	# neutral remembered 
		      9999, 9999,	# neutral not remembered 
		      9999, 9999,	# neutral omission
		      
		      9999, 9999,	# sad 
		      9999, 9999,	# sad remembered 
		      9999, 9999,	# sad not remembered 
		      9999, 9999	# sad omission 
		     );

    } elsif ($stimulusType==$NEUTRAL) {
      if ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $REMEMBERED ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - NEUTRAL REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}

	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_NEUTRAL_RECALL_HITS } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_NEUTRAL_RECALL_HITS } || 0 ) + 1;
	# }

	$line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 

			9999, 9999,	    # fearful remembered
			9999, 9999,	    # fearful
			9999, 9999,	    # fearful not remembered
			9999, 9999,	    # fearful omission

			9999, 9999, # fixation
			
			9999, 9999,	    # happy 
			9999, 9999,	    # happy remembered 
			9999, 9999,	    # happy not remembered 
			9999, 9999,	    # happy omission
			
			9999, 9999, # instructions
			
			$onset, $duration,	   # neutral
			$onset, $duration,  # neutral remembered
			9999, 9999,	    # neutral not remembered
			9999, 9999,	    # neutral omission
			
			9999, 9999,	    # sad 
			9999, 9999,	    # sad remembered 
			9999, 9999,	    # sad not remembered 
			9999, 9999	    # sad omission 
		       );

      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $NOT_REMEMBERED ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - NEUTRAL NOT REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}

	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_NEUTRAL_RECALL_MISSES } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_NEUTRAL_RECALL_MISSES } || 0 ) + 1;
	# }

	$line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 

			9999, 9999,	    # fearful
			9999, 9999,	    # fearful remembered
			9999, 9999,	    # fearful not remembered
			9999, 9999,	    # fearful omission

			9999, 9999, # fixation
			
			9999, 9999,	    # happy
			9999, 9999,	    # happy remembered 
			9999, 9999,	    # happy not remembered 
			9999, 9999,	    # happy omission
			
			9999, 9999, # instructions
			
			$onset, $duration,  # neutral
			9999, 9999,	    # neutral remembered
			$onset, $duration,  # neutral not remembered
			9999, 9999,	    # neutral omission
			
			9999, 9999,	    # sad
			9999, 9999,	    # sad remembered 
			9999, 9999,	    # sad not remembered 
			9999, 9999	    # sad omission 
		       );
	
      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $OMISSION ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - NEUTRAL OMISSION: Target ($targetId) is in facesRecalledHash\n";
	}
	$line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 

			9999, 9999,	    # fearful
			9999, 9999,	    # fearful remembered
			9999, 9999,	    # fearful not remembered
			9999, 9999,	    # fearful omission

			9999, 9999, # fixation
			
			9999, 9999,	    # happy
			9999, 9999,	    # happy remembered 
			9999, 9999,	    # happy not remembered 
			9999, 9999,	    # happy omission
			
			9999, 9999, # instructions
			
			$onset, $duration,  # neutral
			9999, 9999,	    # neutral remembered
			9999, 9999,	    # neutral not remembered
			$onset, $duration,  # neutral omission
			
			9999, 9999,	    # sad 
			9999, 9999,	    # sad remembered 
			9999, 9999,	    # sad not remembered 
			9999, 9999	    # sad omission 
		       );
	
      } else { 
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - NEUTRAL: No entry for target ($targetId) or it's a false alarm in facesRecalledHash\n";
	}
	$line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			9999, 9999,	    # fearful
			9999, 9999,	    # fearful remembered
			9999, 9999,	    # fearful not remembered
			9999, 9999,	    # fearful omission
			
			9999, 9999, # fixation

			9999, 9999,	    # happy
			9999, 9999,	    # happy remembered 
			9999, 9999,	    # happy not remembered 
			9999, 9999,	    # happy omission
			
			9999, 9999, # instructions
			
			$onset, $duration,  # neutral 
			9999, 9999,	    # neutral remembered
			9999, 9999,	    # neutral not remembered
			9999, 9999,	    # neutral omission
			
			9999, 9999,	    # sad
			9999, 9999,	    # sad remembered 
			9999, 9999,	    # sad not remembered 
			9999, 9999	    # sad omission 
		       );
      } 
    } elsif ($stimulusType==$SAD) { 
      if ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $REMEMBERED ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - SAD REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}
	
	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_SAD_RECALL_HITS } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_SAD_RECALL_HITS } || 0 ) + 1;
	# }
	
	$line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%0.3f:%0.3f\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 

			9999, 9999,	      # fearful
			9999, 9999,	      # fearful remembered
			9999, 9999,	      # fearful not remembered
			9999, 9999,	      # fearful omission

			9999, 9999, # fixation
			
			9999, 9999,	      # happy 
			9999, 9999,	      # happy remembered 
			9999, 9999,	      # happy not remembered 
			9999, 9999,	      # happy omission
			
			9999, 9999, # instructions

			9999, 9999,	      # neutral
			9999, 9999,	      # neutral remembered 
			9999, 9999, # neutral not remembered 
			9999, 9999, # neutral omission
 			
			$onset, $duration,	     # sad 
			$onset, $duration,	     # sad remembered
			9999, 9999,	     # sad not remembered
			9999, 9999	     # sad omission
		       );

      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $NOT_REMEMBERED ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - SAD NOT REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}

	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_SAD_RECALL_MISSES } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_SAD_RECALL_MISSES } || 0 ) + 1;
	# }
	
	$line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%0.3f:%0.3f\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 

			9999, 9999,	      # fearful
			9999, 9999,	      # fearful remembered
			9999, 9999,	      # fearful not remembered
			9999, 9999,	      # fearful omission

			9999, 9999, # fixation
			
			9999, 9999,	      # happy 
			9999, 9999,	      # happy remembered 
			9999, 9999,	      # happy not remembered 
			9999, 9999,	      # happy omission
			
			9999, 9999, # instructions

			9999, 9999,	      # neutral
			9999, 9999,	      # neutral remembered 
			9999, 9999, # neutral not remembered 
			9999, 9999, # neutral omission
 			
			$onset, $duration,	     # sad 
			9999, 9999,		     # sad remembered
			$onset, $duration,   # sad not remembered
			9999, 9999	     # sad omission
		       );
	
      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $OMISSION ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - SAD OMISSION: Target ($targetId) is in facesRecalledHash\n";
	}
	$line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%0.3f:%0.3f",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 

			9999, 9999,	      # fearful
			9999, 9999,	      # fearful remembered
			9999, 9999,	      # fearful not remembered
			9999, 9999,	      # fearful omission

			9999, 9999, # fixation
			
			9999, 9999,	      # happy 
			9999, 9999,	      # happy remembered 
			9999, 9999,	      # happy not remembered 
			9999, 9999,	      # happy omission
			
			9999, 9999, # instructions

			9999, 9999,	      # neutral
			9999, 9999,	      # neutral remembered 
			9999, 9999, # neutral not remembered 
			9999, 9999, # neutral omission
 			
			$onset, $duration,    # sad 
			9999, 9999,	      # sad remembered
			9999, 9999,	      # sad not remembered
			$onset, $duration     # sad omission
		       );
	
      } else { 
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - SAD: No entry for target ($targetId) or it's a false alarm in facesRecalledHash\n";
	}

	$line = sprintf("%s\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%d:%d\t%0.3f:%0.3f\t%d:%d\t%d:%d\t%d:%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			9999, 9999,	      # fearful
			9999, 9999,	      # fearful remembered
			9999, 9999,	      # fearful not remembered
			9999, 9999,	      # fearful omission
			
			9999, 9999, # fixation

			9999, 9999,	      # happy
			9999, 9999,	      # happy remembered 
			9999, 9999,	      # happy not remembered 
			9999, 9999,	      # happy omission
			
			9999, 9999, # instructions

			9999, 9999,	      # neutral
			9999, 9999,	      # neutral remembered 
			9999, 9999, # neutral not remembered 
			9999, 9999, # neutral omission
 			
			$onset, $duration,    # sad 
			9999, 9999,	      # sad remembered
			9999, 9999,	      # sad not remembered
			9999, 9999	      # sad omission
		       );
      }
    } else {
      die "Cannot make a regressor line of type $stimulusType\n. Exiting";
    }
    ########################################################################################################################################################
    ########################################################################################################################################################
    ########################################################################################################################################################
    ########################################################################################################################################################
    ########################################################################################################################################################
  } elsif ($regressorType == $ONSET_ONLY_REGRESSOR ) {
    if ($stimulusType==$FEARFUL) {
      if ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $REMEMBERED ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - FEARFUL REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}
		
	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_FEARFUL_RECALL_HITS } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_FEARFUL_RECALL_HITS } || 0 ) + 1;
	# }

	$line = sprintf("%s\t%0.3f\t%0.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			$onset,  # fearful
			$onset,  # fearful remembered
			9999,        # fearful not remembered
			9999,        # fearful omission

			9999, # fixation
			
			9999, # happy 
			9999, # happy remembered
			9999, # happy not remembered
			9999, # happy omission
			
			9999, # instructions
			
			9999, # neutral
			9999, # neutral remembered 
			9999, # neutral not remembered 
			9999, # neutral omission
			
			9999, # sad
			9999, # sad remembered 
			9999, # sad not remembered 
			9999  # sad omission 
		       );

      } elsif ( defined $facesRecalledHash_ref->{ $targetId } &&  $facesRecalledHash_ref->{ $targetId } eq $NOT_REMEMBERED  ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - FEARFUL NOT REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}

	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_FEARFUL_RECALL_MISSES } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_FEARFUL_RECALL_MISSES } || 0 ) + 1;
	# }

	$line = sprintf("%s\t%0.3f\t%d\t%0.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			$onset,  # fearful
			9999,        # fearful remembered
			$onset,  # fearful not remembered
			9999,        # fearful omission
			
			9999, # fixation
			
			9999, # happy 			
			9999, # happy remembered
			9999, # happy not remembered
			9999, # happy omission
			
			9999, # instructions
			
			9999, # neutral
			9999, # neutral remembered 
			9999, # neutral not remembered 
			9999, # neutral omission
			
			9999, # sad 
			9999, # sad remembered 
			9999, # sad not remembered 
			9999  # sad omission 
		       );
	
      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $OMISSION ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - FEARFUL OMISSION: Target ($targetId) is in facesRecalledHash\n";
	}
	
	$line = sprintf("%s\t%0.3f\t%d\t%d\t%0.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			$onset,  # fearful
			9999,        # fearful remembered
			9999,        # fearful not remembered
			$onset,  # fearful omission

			9999, # fixation			

			9999, # happy
			9999, # happy remembered
			9999, # happy not remembered
			9999, # happy omission
			
			9999, # instructions
			
			9999, # neutral
			9999, # neutral remembered 
			9999, # neutral not remembered 
			9999, # neutral omission
			
			9999, # sad 
			9999, # sad remembered 
			9999, # sad not remembered 
			9999  # sad omission 
		       );
	
      } else {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - FEARFUL: No entry for target ($targetId) or it's a false alarm in facesRecalledHash\n";
	}

	$line = sprintf("%s\t%0.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			$onset,  # fearful
			9999,        # fearful remembered
			9999,        # fearful not remembered
			9999,        # fearful omission
			
			9999, # fixation
			
			9999, # happy 
			9999, # happy remembered
			9999, # happy not remembered
			9999, # happy omission
			
			9999, # instructions
			
			9999, # neutral
			9999, # neutral remembered 
			9999, # neutral not remembered 
			9999, # neutral omission
			
			9999, # sad
			9999, # sad remembered 
			9999, # sad not remembered 
			9999  # sad omission 
		       );
      }
    } elsif ($stimulusType==$FIXATION) {
      $line = sprintf("%s\t%d\t%d\t%d\t%d\t%0.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
		      makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
		      
		      9999, # fearful
		      9999, # fearful remembered
		      9999, # fearful not remembered
		      9999, # fearful omission
		      
		      $onset,  # fixation
		      
		      9999, # happy 
		      9999, # happy remembered
		      9999, # happy not remembered
		      9999, # happy omission
		      
		      9999, # instructions
		      
		      9999, # neutral  
		      9999, # neutral remembered 
		      9999, # neutral not remembered 
		      9999, # neutral omission
		      
		      9999, # sad 
		      9999, # sad remembered 
		      9999, # sad not remembered 
		      9999  # sad omission 
		     );
      
    } elsif ($stimulusType==$HAPPY) {
      if ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $REMEMBERED  ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - HAPPY REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}
	
	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_HAPPY_RECALL_HITS } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_HAPPY_RECALL_HITS } || 0 ) + 1;
	# }
	
	$line = sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%0.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			9999, # fearful remembered
			9999, # fearful
			9999, # fearful not remembered
			9999, # fearful omission
			
			9999, # fixation
			
			$onset,  # happy
			$onset,  # happy remembered
			9999,        # happy not remembered
			9999,        # happy omission
			
			9999, # instructions
			
			9999, # neutral
			9999, # neutral remembered 
			9999, # neutral not remembered 
			9999, # neutral omission
			
			9999, # sad
			9999, # sad remembered 
			9999, # sad not remembered 
			9999  # sad omission 
		       );
	
      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $NOT_REMEMBERED  ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - HAPPY NOT REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}
	
	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_HAPPY_RECALL_MISSES } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_HAPPY_RECALL_MISSES } || 0 ) + 1;
	# }
	
	$line = sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%d\t%0.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			9999, # fearful
			9999, # fearful remembered
			9999, # fearful not remembered
			9999, # fearful omission
			
			9999, # fixation
			
			$onset,  # happy
			9999,        # happy remembered
			$onset,  # happy not remembered
			9999,        # happy omission
			
			9999, # instructions
			
			9999, # neutral
			9999, # neutral remembered 
			9999, # neutral not remembered 
			9999, # neutral omission
			
			9999, # sad
			9999, # sad remembered 
			9999, # sad not remembered 
			9999  # sad omission 
		       );
	
      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $OMISSION ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - HAPPY OMISSION: Target ($targetId) is in facesRecalledHash\n";
	}
	$line = sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%d\t%d\t%0.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			9999, # fearful
			9999, # fearful remembered			
			9999, # fearful not remembered
			9999, # fearful omission
			
			9999, # fixation
			
			$onset,  # happy
			9999,        # happy remembered
			9999,        # happy not remembered
			$onset,  # happy omission
			
			9999, # instructions
			
			9999, # neutral
			9999, # neutral remembered 
			9999, # neutral not remembered 
			9999, # neutral omission
			
			9999, # sad
			9999, # sad remembered 
			9999, # sad not remembered 
			9999  # sad omission 
		       );
	
      } else {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - HAPPY: No entry for target ($targetId) or it's a false alarm in facesRecalledHash\n";
	}
	
	$line = sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			9999, # fearful
			9999, # fearful remembered
			9999, # fearful not remembered
			9999, # fearful omission
			
			9999, # fixation
			
			$onset,  # happy 
			9999,	   # happy remembered
			9999,	   # happy not remembered
			9999,	   # happy omission
			
			9999, # instructions
			
			9999, # neutral
			9999, # neutral remembered 
			9999, # neutral not remembered 
			9999, # neutral omission
			
			9999, # sad
			9999, # sad remembered 
			9999, # sad not remembered 
			9999  # sad omission 
		       );
      }
    } elsif ($stimulusType==$INSTRUCTIONS) {
      $line = sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
		      makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
		      
		      9999,	# fearful
		      9999,	# fearful remembered
		      9999,	# fearful not remembered
		      9999,	# fearful omission
		      
		      9999, # fixation
		      
		      9999,	# happy 
		      9999,	# happy remembered
		      9999,	# happy not remembered
		      9999,	# happy omission

		      $onset,  # instructions		      

		      9999,	# neutral 
		      9999,	# neutral remembered 
		      9999,	# neutral not remembered 
		      9999,	# neutral omission
		      
		      9999,	# sad 
		      9999,	# sad remembered 
		      9999,	# sad not remembered 
		      9999	# sad omission 
		     );

    } elsif ($stimulusType==$NEUTRAL) {
      if ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $REMEMBERED ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - NEUTRAL REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}

	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_NEUTRAL_RECALL_HITS } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_NEUTRAL_RECALL_HITS } || 0 ) + 1;
	# }

	$line = sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%0.3f\t%d\t%d\t%d\t%d\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 

			9999,	    # fearful remembered
			9999,	    # fearful
			9999,	    # fearful not remembered
			9999,	    # fearful omission

			9999, # fixation
			
			9999,	    # happy 
			9999,	    # happy remembered 
			9999,	    # happy not remembered 
			9999,	    # happy omission
			
			9999, # instructions
			
			$onset, 	   # neutral
			$onset,   # neutral remembered
			9999,	    # neutral not remembered
			9999,	    # neutral omission
			
			9999,	    # sad 
			9999,	    # sad remembered 
			9999,	    # sad not remembered 
			9999	    # sad omission 
		       );

      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $NOT_REMEMBERED ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - NEUTRAL NOT REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}

	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_NEUTRAL_RECALL_MISSES } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_NEUTRAL_RECALL_MISSES } || 0 ) + 1;
	# }

	$line = sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%d\t%0.3f\t%d\t%d\t%d\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 

			9999,	    # fearful
			9999,	    # fearful remembered
			9999,	    # fearful not remembered
			9999,	    # fearful omission

			9999, # fixation
			
			9999,	    # happy
			9999,	    # happy remembered 
			9999,	    # happy not remembered 
			9999,	    # happy omission
			
			9999, # instructions
			
			$onset,   # neutral
			9999,	    # neutral remembered
			$onset,   # neutral not remembered
			9999,	    # neutral omission
			
			9999,	    # sad
			9999,	    # sad remembered 
			9999,	    # sad not remembered 
			9999	    # sad omission 
		       );
	
      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $OMISSION ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - NEUTRAL OMISSION: Target ($targetId) is in facesRecalledHash\n";
	}
	$line = sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%d\t%d\t%0.3f\t%d\t%d\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 

			9999,	    # fearful
			9999,	    # fearful remembered
			9999,	    # fearful not remembered
			9999,	    # fearful omission

			9999, # fixation
			
			9999,	    # happy
			9999,	    # happy remembered 
			9999,	    # happy not remembered 
			9999,	    # happy omission
			
			9999, # instructions
			
			$onset,   # neutral
			9999,	    # neutral remembered
			9999,	    # neutral not remembered
			$onset,   # neutral omission
			
			9999,	    # sad 
			9999,	    # sad remembered 
			9999,	    # sad not remembered 
			9999	    # sad omission 
		       );
	
      } else { 
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - NEUTRAL: No entry for target ($targetId) or it's a false alarm in facesRecalledHash\n";
	}
	$line = sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			9999,	    # fearful
			9999,	    # fearful remembered
			9999,	    # fearful not remembered
			9999,	    # fearful omission
			
			9999, # fixation

			9999,	    # happy
			9999,	    # happy remembered 
			9999,	    # happy not remembered 
			9999,	    # happy omission
			
			9999, # instructions
			
			$onset,   # neutral 
			9999,	    # neutral remembered
			9999,	    # neutral not remembered
			9999,	    # neutral omission
			
			9999,	    # sad
			9999,	    # sad remembered 
			9999,	    # sad not remembered 
			9999	    # sad omission 
		       );
      } 
    } elsif ($stimulusType==$SAD) { 
      if ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $REMEMBERED ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - SAD REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}
	
	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_SAD_RECALL_HITS } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_SAD_RECALL_HITS } || 0 ) + 1;
	# }
	
	$line = sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%0.3f\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 

			9999,	      # fearful
			9999,	      # fearful remembered
			9999,	      # fearful not remembered
			9999,	      # fearful omission

			9999, # fixation
			
			9999,	      # happy 
			9999,	      # happy remembered 
			9999,	      # happy not remembered 
			9999,	      # happy omission
			
			9999, # instructions

			9999,	      # neutral
			9999,	      # neutral remembered 
			9999, # neutral not remembered 
			9999, # neutral omission
 			
			$onset, 	     # sad 
			$onset, 	     # sad remembered
			9999,	     # sad not remembered
			9999	     # sad omission
		       );

      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $NOT_REMEMBERED ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - SAD NOT REMEMBERED: Target ($targetId) is in facesRecalledHash\n";
	}

	# if ( $phase == $ACQUISITION ) {
	#   $recallSummaryHash_ref->{ $TOTAL_ACQ_SAD_RECALL_MISSES } = ( $recallSummaryHash_ref->{ $TOTAL_ACQ_SAD_RECALL_MISSES } || 0 ) + 1;
	# }
	
	$line = sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%d\t%0.3f\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 

			9999,	      # fearful
			9999,	      # fearful remembered
			9999,	      # fearful not remembered
			9999,	      # fearful omission

			9999, # fixation
			
			9999,	      # happy 
			9999,	      # happy remembered 
			9999,	      # happy not remembered 
			9999,	      # happy omission
			
			9999, # instructions

			9999,	      # neutral
			9999,	      # neutral remembered 
			9999, # neutral not remembered 
			9999, # neutral omission
 			
			$onset, 	     # sad 
			9999,		     # sad remembered
			$onset,    # sad not remembered
			9999	     # sad omission
		       );
	
      } elsif ( defined $facesRecalledHash_ref->{ $targetId } && $facesRecalledHash_ref->{ $targetId } eq $OMISSION ) {
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - SAD OMISSION: Target ($targetId) is in facesRecalledHash\n";
	}
	$line = sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%d\t%d\t%0.3f",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 

			9999,	      # fearful
			9999,	      # fearful remembered
			9999,	      # fearful not remembered
			9999,	      # fearful omission

			9999, # fixation
			
			9999,	      # happy 
			9999,	      # happy remembered 
			9999,	      # happy not remembered 
			9999,	      # happy omission
			
			9999, # instructions

			9999,	      # neutral
			9999,	      # neutral remembered 
			9999, # neutral not remembered 
			9999, # neutral omission
 			
			$onset,     # sad 
			9999,	      # sad remembered
			9999,	      # sad not remembered
			$onset,     # sad omission
		       );
	
      } else { 
	if ($DEBUG) {
	  print "### makeStimulusRegressorLine - SAD: No entry for target ($targetId) or it's a false alarm in facesRecalledHash\n";
	}

	$line = sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%0.3f\t%d\t%d\t%d",
			makeInstructionsInEffectRegressorComponent($instructionsInEffect), 
			
			9999,	      # fearful
			9999,	      # fearful remembered
			9999,	      # fearful not remembered
			9999,	      # fearful omission
			
			9999, # fixation

			9999,	      # happy
			9999,	      # happy remembered 
			9999,	      # happy not remembered 
			9999,	      # happy omission
			
			9999, # instructions

			9999,	      # neutral
			9999,	      # neutral remembered 
			9999, # neutral not remembered 
			9999, # neutral omission
 			
			$onset,     # sad 
			9999,	      # sad remembered
			9999,	      # sad not remembered
			9999,	      # sad omission
		       );
      }
    } else {
      die "Cannot make a regressor line of type $stimulusType\n. Exiting";
    }

  } elsif ($regressorType == $ONSET_OFFSET_REGRESSOR ) {
    die "### makeStimulusRegressorLine - Can't make offset regressors. Write the code for that!\n";
  }
  
  if ($DEBUG) {
    print "### EXIT makeStimulusRegressorLine\n";
  }
  
  return ($line);
}

sub extractLogFrames($) {
  # print "extractLogFrames: ### extractLogFrames start\n";

  my $pineDataFile=shift;

  my @logFrameEntries = ();

  my @logFrameEntry = ();
  my $numberOfResets=0;
  while ( <$pineDataFile> ) {
    chomp $_;
    # print "line just read is: $_\n";
    if ( /\*\*\* LogFrame Start \*\*\*/ .. /\*\*\* LogFrame End \*\*\*/ ) {
      #print "extractLogFrames: @@@@@@@@@@ Line: $_\n";
      push @logFrameEntry, $_
    } else {
      # print "extractLogFrames: !!!!!!!!!!!!!!!! Log frame processing ending\n";
      # my $length = $#logFrameEntry + 1;
      # print "extractLogFrames: !!!!!!!!!!!!!!!! Log frame has $length entries\n";
      # print "extractLogFrames: This is entry $numberOfResets\n";
      # for my $i (0 .. $#logFrameEntry) {
      # 	print "extractLogFrames: %%%%%%%%%%%%%% Log line $i: " . $logFrameEntry[$i] . "\n";
      # }

      # print "extractLogFrames: Length: " . @logFrameEntry . "\n";
      # add  the log frame to the array of log frames
      push @logFrameEntries,   [ @logFrameEntry  ] ;

      ## reset the array for the next log entry
      @logFrameEntry = ();
      $numberOfResets=$numberOfResets+1;
    }
  }				## end of while
  
  ## this handles the final LogFrame Entry. This is necessary as after
  ## the end of the file is reached the, the else branch of the if
  ## statement above never executes as there is nothing more ing th file
  ## to satisfy the guard of the while loop

  {				## anonymous block 
				## print "extractLogFrames: !!!!!!!!!!!!!!!! Log frame processing ending\n";
				## my $length = $#logFrameEntry + 1;
				## print "extractLogFrames: !!!!!!!!!!!!!!!! Log frame has $length entries\n";
				## print "extractLogFrames: This is entry $numberOfResets\n";
				## for my $i (0 .. $#logFrameEntry) {
				##   print "extractLogFrames: %%%%%%%%%%%%%% Log line $i: " . $logFrameEntry[$i] . "\n";
				## }
    
				# print "extractLogFrames: Length: " . @logFrameEntry . "\n";
				# add  the log frame to the array of log frames
    push @logFrameEntries,   [ @logFrameEntry  ] ;
    
    ## reset the array for the next log entry
    @logFrameEntry = ();
    $numberOfResets=$numberOfResets+1;
  }

  # print "extractLogFrames: Number of resets: $numberOfResets\n";
  # my $length = $#logFrameEntries + 1;

  #for my $i (0 .. $length ) {
  #  print "@{ $logFrameEntries[$i] }\n";
  #}
  # print "extractLogFrames: Returning $length log frame entries\n";
  # print "extractLogFrames: ### extractLogFrames end\n";

  return (\@logFrameEntries);
}



sub checkSubjectBehavioralDataDirectoryExists($) {
  my $subjectId=shift;
  my $returnValue="0";
  
  my $subjectFolderRoot = "$MDD_ROOT/$subjectId/$subjectId"; # . "_beh";
  my @folderPossibilities= (
			    "$MDD_ROOT/$subjectId/$subjectId" . "_beh", "$MDD_ROOT/$subjectId/$subjectId" . "_Beh"
			   );
  my $oldSubjectId=$subjectId;
  $subjectId =~ s/_[ABCD]//;
  push @folderPossibilities, "$MDD_ROOT/$oldSubjectId/$subjectId" .  "_beh";
  push @folderPossibilities, "$MDD_ROOT/$oldSubjectId/$subjectId" .  "_Beh";
  
    
    foreach my $possibility ( @folderPossibilities ) {
      if ($DEBUG) {
	print "### checkSubjectBehavioralDataDirectoryExists - Trying $possibility\n";
      }
      if ( -d $possibility ) {
	$returnValue=$possibility;
	last;
      }
    }
  
  return $returnValue;
}

sub findSubjectBehavioralDataFile($$$) {
  my $subjectId=shift;
  my $subjectFolder=shift;
  my $phase=shift;

  #my $subjectFolder = "$MDD_ROOT/$subjectId/$subjectId" . "_beh";

  opendir(DIR, $subjectFolder) or die $!;

  my @candidateFiles=();
  if ($phase == $ACQUISITION) {
    @candidateFiles= grep { 
      /fac2sd-.*-\d*\.txt/		  # Begins with a period
	&& -f "$subjectFolder/$_"	  # and is a file
      } readdir(DIR);
  }
  elsif ($phase == $RECALL) {
    @candidateFiles= grep { 
      /rcl2sd-.*-\d*\.txt/		  # Begins with a period
	&& -f "$subjectFolder/$_"	  # and is a file
      } readdir(DIR);
  } else {
    die "Unable to handlle that phase. Please check your command line and ensure you have selected either acuqisition or recall. Exiting\n";
  }
  #if ($DEBUG) {
  print "Candidate Files are (in reverse order): ";
  @candidateFiles = sort {$b cmp $a} @candidateFiles;
  print "@candidateFiles\n";
  
  if ( scalar(@candidateFiles) > 1 ) { 
    print "Found more than one file returning only the first one\n";
  }
  #}
  
  return $candidateFiles[0];
}

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
  {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
  }

# Left trim function to remove leading whitespace
sub ltrim($)
  {
    my $string = shift;
    $string =~ s/^\s+//;
    return $string;
  }

# Right trim function to remove trailing whitespace
sub rtrim($)
  {
    my $string = shift;
    $string =~ s/\s+$//;
    return $string;
  }
