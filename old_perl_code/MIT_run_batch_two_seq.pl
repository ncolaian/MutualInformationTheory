#! /usr/bin/env perl

# This script will run MIT on all of the protein combinations

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Data::Dumper;

#My Variables
my $help = 0;
my $man = 0;
my $seq_file_ext;
my $out_dir;
my $file_dir;
my $mit_path;

#Read in the variables from the command line
GetOptions( 'man|m'   =>  \$man,
            'help|h'  =>  \$help,
            'file_ext|e=s' => \$seq_file_ext,
            'out_dir|o=s'  => \$out_dir,
            'file_path|f=s'   => \$file_dir,
            'mit_path|m=s' => \$mit_path,
            ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manual and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

## Main ##
make_ordered_files();
#run_MIT();

## Subroutines ##
sub run_MIT {
   my $files = `ls $out_dir/ordered_files`;
   `mkdir $out_dir/out_files`;
   `mkdir $out_dir/mit_res`;
   my @files_a = split /\n/, $files;
   
   my @comparisons;
   for ( my $i = 0; $i<scalar(@files_a); $i++ ){
    for ( my $j = $i+1; $j < scalar(@files_a); $j++ ) {
      #get the correct names
      my @first_path = (split /\//, $files_a[$i]);
      my @second_path = (split /\//, $files_a[$j]);
      my $first = $first_path[scalar(@first_path)-1];
      my $second = $second_path[scalar(@second_path)-1];
      
      $first = (split /\./, $first)[0] . "_" . (split /\./, $first)[1];
      $second = (split /\./, $second)[0] . "_" . (split /\./, $second)[1];
      #put together names
      my $combined = $first . "-" . $second;
      #push to holder
      push @comparisons, $combined;
      `mkdir $out_dir/mit_res/$combined`;
      my $file_one = $files_a[$i];
      my $file_two = $files_a[$j];
      chomp $file_one;
      chomp $file_two;
      `sbatch -o $out_dir/out_files/$combined.out --wrap="perl $mit_path/MIT_two_proteins.pl -aa -mo $out_dir/ordered_files/$file_one -mt $out_dir/ordered_files/$file_two -o $out_dir/mit_res/$combined"`;
    }
   }
   
   return(1)
}
sub make_ordered_files {
   my $files = `ls $file_dir/*$seq_file_ext`;
   my @files_a = split /\n/, $files;
   open(my $order_IN, "<", $files_a[0]) or die "First file doesn't exist\n";
   `mkdir $out_dir/ordered_files`;
   
   while ( <$order_IN>) {
    if ($_ !~ /^>/) {
        next;
    }
    chomp $_;
    my @grep_genome_a = split /-/, $_;
    my $grep_genome = $grep_genome_a[scalar(@grep_genome_a)-1];
    
    my @genome;
    my @grep_print;
    my $print_seq = 1;
    foreach my $n (@files_a) {
      my @N = split /\//, $n;
      my $n_new = $N[scalar(@N)-1];
      my $grep = `grep -A 1 -- "-$grep_genome" $n`;
      if ($grep !~ />/) {
        $print_seq=0;
      }
      else{
         push @genome, $n_new;
         push @grep_print, $grep;
      }
    }
    if ($print_seq) {
        for ( my $i = 0; $i< scalar(@genome); $i++) {
         open(my $OUT, ">>", "$out_dir/ordered_files/$genome[$i]") or die "can't add grep seqs\n";
         print $OUT $grep_print[$i];
         close $OUT;
        }
    }
    
   }
   close $order_IN;
   return(1);
}



__END__
=head1 TITLE



=head1 VERSION



=head1 INCLUDED MODULES

Getopt::Long;
Pod::Usage;
Carp;
Readonly;
Path::Class;
Data::Dumper;
Log::Log4perl qw(:easy);
Log::Log4perl::CommandLine qw(:all);

=head1 INHERIT

=head1 SYNOPSIS

=head1 PARAMETERS

=head1 CONFIGURATION AND ENVIRONMENT



=head1 DEPENDENCIES


    
=head1 INCOMPATIBILITIES

    None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests	
	
=head1 AUTHOR

Nicholas Colaianni
contact via C<< <ncolaian@live.unc.edu> >>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2018, Nicholas Colaianni
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut