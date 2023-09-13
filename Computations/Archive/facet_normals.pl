#!/usr/bin/perl

use strict;
use warnings;

use application "fan";

my $n = $ARGV[0];
print "n: $n\n";
my $input = "SchlaefliFan$n.json";
my $mpha_result = "schlaefli_mc_$n.xz";
print "Files: $input $mpha_result\n";
my $arrangement = load($input);

open my $pipe, "xzcat $mpha_result |";
my $line;
while($line = <$pipe> ){
   print $line;
   my($sig, $rays) = $line =~ m/Signature: (\[[0-9\-,]*\]) Rays: (\[\[.*\]\])/;
   print $sig,"\n";
   $rays = new Matrix<Rational>(eval $rays);
   print $rays;
}

