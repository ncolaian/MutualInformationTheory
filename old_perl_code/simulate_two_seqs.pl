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
use List::Util qw(sum);

#My Variables
my $help = 0;
my $man = 0;
my $seq_out;
my $num_seqs;
my $ident_1;
my $ident_2;
my $avg_ident_out; #this will be the avg ident 
my $seq_length = 100;

#Read in the variables from the command line
GetOptions( 'man'   =>  \$man,
            'help'  =>  \$help,
            'seq_out|o=s'   =>  \$seq_out,
            'num_seqs|n=i'  =>  \$num_seqs,
            'ident_s1|is=i'  =>  \$ident_1,
            'ident_s2|it=i' =>  \$ident_2,
            'avg_out|a=s'    =>  \$avg_ident_out,
            'seq_length|l:i'    => \$seq_length,
    ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manal and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

#### SETTING EVERYTHING UP ####
#array for big random
Readonly::Array my @AMINO_ACIDS => qw/
A R N D C Q E G H I L K M F P S T W Y V
    /;
    ### BASED ON WAG MODEL ###
Readonly::Array my @AMINO_ACID_FREQ => qw/
0.0866279 0.043972 0.0390894 0.0570451 0.0193078 0.0367281 0.0580589 0.0832518 0.0244313 0.048466 0.086209 0.0620286 0.0195027 0.0384319 0.0457631 0.0695179 0.0610127 0.0143859 0.0352742 0.0708956
    /;
    
    #wag transitions from https://www.ebi.ac.uk/goldman-srv/WAG/wag.dat
    ### ASK CORBIN ABOUT THIS ###
#Readonly::Hash my %TRANSITIONS => (
#    "A" =>  qw/0.551571 0.509848 0.738998 1.027040 0.908598 1.582850 1.416720 0.316954 0.193335 0.397915 0.906265 0.893496 0.210494 1.438550 3.370790 2.121110 0.113133 0.240735 2.006010/,
#    "R" =>  qw/0.509848  0.635346 0.147304 0.528191 3.035500 0.439157 0.584665 2.137150 0.186979 0.497671 5.351420 0.683162 0.102711 0.679489 1.224190 0.554413 1.163920 0.381533 0.251849/,
#    "N" =>  qw/0.738998  0.147304  5.429420 0.265256 1.543640 0.947198 1.125560 3.956290 0.554236 0.131528 3.012010 0.198221 0.0961621 0.195081 3.974230 2.030060 0.0719167 1.086000 0.196246/,
#    "D" =>  qw/1.027040  0.528191  0.265256  0.0302949 0.616783 6.174160 0.865584 0.930676 0.039437 0.0848047 0.479855 0.103754 0.0467304 0.423984 1.071760 0.374866 0.129767 0.325711 0.152335/,
#    "C" =>  qw/0.908598  3.035500  1.543640  0.616783  0.0988179 0.021352 0.306674 0.248972 0.170135 0.384287 0.0740339 0.390482 0.398020 0.109404 1.407660 0.512984 0.717070 0.543833 1.002140/,
#    "Q" =>  qw/1.582850  0.439157  0.947198  6.174160  0.021352  5.469470 0.330052 4.294110 0.113917 0.869489 3.894900 1.545260 0.0999208 0.933372 1.028870 0.857928 0.215737 0.227710 0.301281/,
#    "E" =>  qw/1.416720  0.584665  1.125560  0.865584  0.306674  0.330052  0.567717 0.570025 0.127395 0.154263 2.584430 0.315124 0.0811339 0.682355 0.704939 0.822765 0.156557 0.196303 0.588731/,
#    "G" =>  qw/0.316954  2.137150  3.956290  0.930676  0.248972  4.294110  0.570025  0.249410 0.0304501 0.0613037 0.373558 0.174100 0.049931 0.243570 1.341820 0.225833 0.336983 0.103604 0.187247/,
#    "H" =>  qw/0.193335  0.186979  0.554236  0.039437  0.170135  0.113917  0.127395  0.0304501 0.138190 0.499462 0.890432 0.404141 0.679371 0.696198 0.740169 0.473307 0.262569 3.873440 0.118358/,
#    "I" =>  qw/0.397915  0.497671  0.131528  0.0848047 0.384287  0.869489  0.154263  0.0613037 0.499462  3.170970 0.323832 4.257460 1.059470 0.0999288 0.319440 1.458160 0.212483 0.420170 7.821300/,
#    "L" =>  qw/0.906265  5.351420  3.012010  0.479855  0.0740339 3.894900  2.584430  0.373558  0.890432  0.323832  0.257555 4.854020 2.115170 0.415844 0.344739 0.326622 0.665309 0.398618 1.800340/,
#    "K" =>  qw/0.893496  0.683162  0.198221  0.103754  0.390482  1.545260  0.315124  0.174100  0.404141  4.257460  4.854020  0.934276 0.088836 0.556896 0.967130 1.386980 0.137505 0.133264 0.305434/,
#    "M" =>  qw/0.210494  0.102711  0.0961621 0.0467304 0.398020  0.0999208 0.0811339 0.049931  0.679371  1.059470  2.115170  0.088836  1.190630 0.161444 0.545931 0.171903 1.529640 6.454280 0.649892/,
#    "F" =>  qw//,
#    "P" =>  qw//,
#    "S" =>  (),
#    "T" =>  (),
#    "W" =>  (),
#    "Y" =>  (),
#    "V" =>  ()
#    );
    
my @AA_FRACTIONS;
my $count = 0.000;
for ( my $i = 0; $i < scalar(@AMINO_ACID_FREQ); $i++ ) {
    $count += $AMINO_ACID_FREQ[$i];
    push @AA_FRACTIONS, $count;
}

#my variability across the sequence(conservation decreases across the sequence)
my @CONSERVATION_ARRAY;
for ( my $i = 0; $i < $seq_length; $i++ ) {
    my $rep = int($i/4)+1;
    my @add_array = ($i)x$rep;
    push @CONSERVATION_ARRAY, @add_array;
}

## MAIN ##
#I will be creating paired sequences to perform MIT on
#The sequences will all have an identity of a certain percentage that will be recorded and will
#be very close to what is passed into the script.
#I will output the exact statistic.
# I am looking to see what kind of things I will be able to identify using MIT
# and at what point will I not be able to recognize it.
# I will be performing this at a lot of different levels

#I am for now going to manually set seq_length to 100 for ease of calculations

random_main();

## SUBROUTINES ##

sub random_main {
    my @start_seq1;
    my @start_seq2;
    for ( my $i = 0; $i < $seq_length; $i++ ) {
        push @start_seq1, get_rand_amino_acid();
        push @start_seq2, get_rand_amino_acid();
    }
    my @set_seq1;
    my @set_seq2;
    push @set_seq1, join("",@start_seq1);
    push @set_seq2, join("", @start_seq2);
    
    open my $OUT, ">", "$seq_out.1";
    open my $OUT2, ">", "$seq_out.2";
    print $OUT ">Seq0\n" . join("", @start_seq1) . "\n";
    print $OUT2 ">Seq0\n" . join("", @start_seq2) . "\n";
    #create a sequence for each file
    for ( my $i = 1; $i < $num_seqs; $i++ ) {
        my $seq1_new = create_diff_seq(\@start_seq1, $ident_1);
        my $seq2_new = create_diff_seq(\@start_seq2, $ident_2);
        
        #THIS IS WHERE YOU WOULD ADD SOME MIT changes
        
        push @set_seq1, $seq1_new;
        push @set_seq2, $seq2_new;
        
        print $OUT ">Seq$i\n$seq1_new\n";
        print $OUT2 ">Seq$i\n$seq2_new\n";
    }
    
    close $OUT;
    close $OUT2;
    #get the average sequence identity for the sequences
    my $ident_file_1 = get_avg_seq_identity(\@set_seq1);
    my $ident_file_2 = get_avg_seq_identity(\@set_seq2);
    
    open my $IDENT, ">>", $avg_ident_out;
    print $IDENT "$seq_out.1\t$ident_file_1\n$seq_out.2\t$ident_file_2\n";
    close $IDENT;
}

#based on wag amino acid frequencies
sub get_rand_amino_acid {
    my $rand_num = rand();
    for ( my $i = 0; $i < scalar(@AA_FRACTIONS); $i++ ) {
        if ( $rand_num < $AA_FRACTIONS[$i] ) {
            return($AMINO_ACIDS[$i]);
        }
    }
}

#This will assume a 100 length sequence that has conserved residues and some that are variable
sub create_diff_seq {
    my ($orig_seq_aref, $ident_percentage) = @_;
    my $num_changes = int(((scalar(@$orig_seq_aref) - $ident_percentage)/scalar(@$orig_seq_aref))*100);
    
    my @new_seq = @$orig_seq_aref;
    for ( my $i = 0; $i < $num_changes; $i++ ) {
        #this is where could use transtitions instead of changes
        $new_seq[$CONSERVATION_ARRAY[int(rand(scalar(@CONSERVATION_ARRAY)))]] = get_rand_amino_acid();
    }
    return(join("", @new_seq));
}

sub get_avg_seq_identity {
    my ( $seq_aref ) = @_;
    
    my @ident;
    for ( my $i = 0; $i < scalar(@$seq_aref); $i++ ) {
        for ( my $j = $i+1; $j < scalar(@$seq_aref); $j++ ) {
            push @ident, calc_identity($seq_aref->[$i], $seq_aref->[$j]);
        }
    }
    return(sum(@ident)/scalar(@ident));
}

sub calc_identity {
    my ( $seq1, $seq2 ) = @_;
    my @seq1_a = split //, $seq1;
    my @seq2_a = split //, $seq2;
    
    my $sim = 0;
    for ( my $i = 0; $i < scalar(@seq1_a); $i++ ) {
        if ( $seq1_a[$i] eq $seq2_a[$i] ) {
            $sim++;
        }
    }
    return( $sim/scalar(@seq1_a) );
}
