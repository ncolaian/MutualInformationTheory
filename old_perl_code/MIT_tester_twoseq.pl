#! /usr/bin/env perl

# This script will create a random file containing a number of sequence and combinations used to test MIT

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Data::Dumper;
use Readonly;

#My Variables
my $help = 0;
my $man = 0;
my $out;
my $seqs;
my $null = 0;
my $real_null;

#Read in the variables from the command line
GetOptions( 'man'   =>  \$man,
            'help'  =>  \$help,
            'out|o=s'   =>  \$out,
            'seqs|s=i'  =>  \$seqs,
            'null|n:i'  =>  \$null,
            'real_n|r:i'    =>  \$real_null,
    ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manal and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

#array for big random
Readonly::Array my @AMINO_ACIDS => qw/
A R N D C Q E G H I L K M F P S T W Y V
    /;

## Main ##
open my $OUT, ">", $out;
open my $OUT2, ">", "$out.2";
#sequences will be 20 amino acids long
for (my $i = 0; $i < $seqs; $i++ ){
    my @array;
    my @array_2;
    
    if ( $real_null ) {
        $array[$real_null-1] = "A";
        @array = ("A") x @array;
        my $changes = 1+int(rand(2));
        my @changes_ar;
        $changes_ar[$changes-1] = 0;
        foreach my $c ( @changes_ar ) {
            $c = int(rand($real_null));
            $array[$c] = $AMINO_ACIDS[rand @AMINO_ACIDS];
        }
        print $OUT ">Seq$i\n" . join("", @array) . "\n";
    }
    
    elsif ( $null ) {
        $array[$null-1] = "A";
        for ( my $i = 0; $i < $null; $i++ ) {
            $array[$i] = $AMINO_ACIDS[rand @AMINO_ACIDS];
        }
        print $OUT ">Seq$i\n" . join("", @array) . "\n";
    }
    else {
        $array[19] = "A";
        $array_2[19] = "A";
        @array = ("A") x @array;
        @array_2 = ("A") x @array;
        
        #amino acids at position 2 and 19 will be linked to change the same time to the same random amino acid
        my $rand = $AMINO_ACIDS[rand @AMINO_ACIDS];
        $array[1] = $rand;
        $array_2[18] = $rand;
        
        #amino acids 4 and 17 will both change every time but to random amino acids
        $array[3] = $AMINO_ACIDS[rand @AMINO_ACIDS];
        $array_2[16] = $AMINO_ACIDS[rand @AMINO_ACIDS];
        
        #amino acids 6 and 15 will change randomly at about 50% of the time to the same amino acid
        if ( rand() > .5 ) {
            $rand = $AMINO_ACIDS[rand @AMINO_ACIDS];
            $array[5] = $rand;
            $array_2[14] = $rand;
        }
        
        #amino acids 8 and 13 will change randomly at about 50% of the time to random amino acids
        if ( rand() > .5 ) {
            $array[7] = $AMINO_ACIDS[rand @AMINO_ACIDS];
            $array_2[12] = $AMINO_ACIDS[rand @AMINO_ACIDS];
        }
        
        #amino acids 10 and 11 will change 50% and to the same amino acid every time
        if ( rand() > .5 ) {
            $array[9] = "N";
            $array_2[10] = "N";
        }
        print $OUT ">Seq$i\n" . join("", @array) . "\n";
        print $OUT2 ">Seq$i\n" . join("", @array_2) . "\n";
    }
}

close($OUT);
close($OUT2);
