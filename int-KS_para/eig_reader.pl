#! /usr/bin/perl
#
use strict;
use warnings;

#Name of the eigenvalue output file From QE, in the temp file.  NOT direct QE output!!
my $eig_input = shift @ARGV;

exit unless (-f $eig_input);

open my $fh, '<', $eig_input
   or die " ERROR: Cannot not open $eig_input $!";

my @eigenvalues;
for my $line (<$fh>){

   next if ($line =~ /STEP:/);
   if ($line =~ /Eigenvalues/){
      @eigenvalues = ();
      next;
   }

   push @eigenvalues, $_ for ( split ' ', $line ); 

}
close($fh);


print "$_ \n" for (@eigenvalues);

