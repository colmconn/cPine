#!/Volumes/opt/local/bin/perl

use strict;
use warnings;

my %myhash=();

## $recallSummaryHash{ $FEARFUL_RECALL_HITS } = $recallSummaryHash{ $FEARFUL_RECALL_HITS } + 1;

$myhash{"key1"} = $myhash{"key1"} + 1;
$myhash{"key2"} = ( $myhash{"key2"} || 0 ) + 1;

print "key1 => " . $myhash{"key1"} . "\n";
print "key2 => " . $myhash{"key2"} . "\n";

$myhash{"key1"} = $myhash{"key1"} + 1;
print "key1 => " . $myhash{"key1"} . "\n";
