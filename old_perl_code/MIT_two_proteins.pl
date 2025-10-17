#! /usr/bin/env perl

# This will create a matrix that can be used for MIT
# Need to make sure the genomes where the proteins come from match.

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
my $aa;
my $msa_1;
my $msa_2;
my $out;

#Read in the variables from the command line
GetOptions( 'man'   =>  \$man,
            'help'  =>  \$help,
            'aa'    =>  \$aa,
            'msa_1|mo=s'   =>  \$msa_1,
            'msa_2|mt=s'    =>  \$msa_2,
            'out|o=s'   =>  \$out
    ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manal and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

Readonly::Array my @AMINO_ACIDS => qw/
A R N D C Q E G H I L K M F P S T W Y V
    /;

Readonly::Array my @N_ACIDS => qw/ A T G C /;

main2();

## Subroutines ##

sub main2{
    my $seq_aref_1 = get_sequences($msa_1);
    my $pos_href_1 = get_hash_pos($seq_aref_1);
    my $seq_aref_2 = get_sequences($msa_2);
    my $pos_href_2 = get_hash_pos($seq_aref_2);
    create_matrix($pos_href_1,$pos_href_2);
    perform_MIT_type2();
}

sub get_sequences {
    my ($msa) = @_;
    open my $IN, "<", $msa;
    my $thr_awy = <$IN>;
    
    my @seq;
    my $mm_seq = "";
    while ( <$IN> ) {
        chomp $_;
        if ( $_ =~ /^>/ ) {
            #handle first line
            if ( $mm_seq ne "" ) {
                push @seq, $mm_seq;
                $mm_seq = "";
                next;
            }
            next;
        }
        $mm_seq .= $_;
    }
    #handle last line
    push @seq, $mm_seq;
    close($IN);
    return( \@seq );
}
sub get_hash_pos {
    my ( $seq_aref ) = @_;
    
    my %pos_hash;
    foreach my $seq ( @$seq_aref ) {
        my @bases = split //, $seq;
        
        my $pos = 1;
        foreach my $b ( @bases ) {
            if ( $pos_hash{$pos} ) {
                push @{$pos_hash{$pos}}, $b;
            }
            else {
                $pos_hash{$pos} = [$b];
            }
            $pos++;
        }
    }
    
    return( \%pos_hash );
}

sub create_matrix {
    my ( $pos_href_1, $pos_href_2 ) = @_;
    
    open my $OUT, ">", "$out/mit_mat_full.txt" || die "can't open outfile\n";
    print $OUT "Seq1Pos\tSeq2Pos\t";
    #create a file containing the hash
    my %out_pos;
    my @print_out;
    if ( $aa ) {
        my $posit = 0;
        foreach my $aa ( @AMINO_ACIDS ) {
            foreach my $aa2 ( @AMINO_ACIDS ) {
                my $key = "$aa-$aa2";
                $out_pos{$key} = $posit;
                $posit++;
                push @print_out, $key;
            }
        }
    }
    
    else {
        my $posit = 0;
        foreach my $na ( @N_ACIDS ) {
            foreach my $na2 ( @N_ACIDS ) {
                my $key = "$na-$na2";
                $out_pos{$key} = $posit;
                $posit++;
                push @print_out, $key;
            }
        }
    }
    print $OUT join("\t", @print_out) . "\n";
    
    
    for ( my $i = 1; $i <= keys %$pos_href_1; $i++ ) {
        for ( my $j = 1; $j <= keys %$pos_href_2; $j++) {
            my @array;
            #I changed the 0 to 2 to handle low number of counts
            $array[(keys %out_pos)-1] = 0;
            @array = (0) x @array;
            my $pos = 0;
            foreach my $base ( @{$pos_href_1->{$i}} ) {
                my $href_i_aref = $pos_href_1->{$i};
                my $href_j_aref = $pos_href_2->{$j};
                if ( $href_i_aref->[$pos] eq "-" || $href_j_aref->[$pos] eq "-" || $href_i_aref->[$pos] eq "X" || $href_j_aref->[$pos] eq "X" || $href_j_aref->[$pos] eq "*" || $href_i_aref->[$pos] eq "*") {
                    $pos++;
                    next;
                }
                #should put a try catch block here to handle errors *HANDLED ABOVE*
                my $tag = $href_i_aref->[$pos] . "-" . $href_j_aref->[$pos];
                $array[$out_pos{$tag}]++;
                $pos++;
            }
            print $OUT "$i\t$j\t";
            print $OUT join("\t", @array);
            print $OUT "\n";
        }
    }
    
    close($OUT);
    return(1);
}

#Small Subroutine
sub array_sum {
    my ($aref) = @_;
    my $sum = 0;
    foreach my $part ( @$aref) {
        $sum += $part;
    }
    return $sum;
}

sub perform_MIT_type2 {
    open my $MATRIX, "<", "$out/mit_mat_full.txt";
    open my $OUT, ">", "$out/mit.txt";
    print $OUT "Pos1\tPos2\tEntropy\n";
    
    my $header = <$MATRIX>;
    chomp $header;
    my @header_split = split /\t/, $header;
    
    my %pos;
    for ( my $i = 2; $i < scalar(@header_split); $i++ ) {
        my @mm_array = (split /-/, $header_split[$i]);
        $pos{$i} = \@mm_array;
    }
    
    while ( <$MATRIX> ) {
        chomp $_;
        my @split_line = split /\t/, $_;

        my $pos1 = shift @split_line;
        my $pos2 = shift @split_line;
        
        #start 2 counters for p(ai) and p(bj)
        my %pai;
        my %pbj;
        if ( $aa ) {
            foreach my $base ( @AMINO_ACIDS ) {
                $pai{$base} = 0.00;
                $pbj{$base} = 0.00;
            }
        }
        else {
            foreach my $base ( @N_ACIDS ) {
                $pai{$base} = 0.00;
                $pbj{$base} = 0.00;
            }
        }
        
        #now need to perform MIT
        my $total_comb = array_sum(\@split_line);
        if ( $total_comb == 0 ) {
            next;
        }
        
        my $total_mit = 0.000000;
        #first go through and get the per base probability
        for ( my $i = 0; $i < scalar(@split_line); $i++ ) {
            next if ( $split_line[$i] == 0 );
            my @bases = @{$pos{$i+2}}; #go through the paired amino acids
            
            my $comb_prob = $split_line[$i]/$total_comb;
            $pai{$bases[0]} += $comb_prob;
            $pbj{$bases[1]} += $comb_prob;
        }
        for ( my $i = 0; $i < scalar(@split_line); $i++ ) {
            next if ( $split_line[$i] == 0 );
            
            my @bases = @{$pos{$i+2}}; #go through the paired amino acids
            
            my $comb_prob = $split_line[$i]/$total_comb;
            #perform the math
            $total_mit += $comb_prob*log( $comb_prob/($pai{$bases[0]} * $pbj{$bases[1]}) );
        }
        print $OUT "$pos1\t$pos2\t$total_mit\n";
    }
    close $MATRIX;
    close $OUT;
}
