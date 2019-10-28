#!/usr/bin/perl -w

###########################################
# PAPAC_diveristy.pl
#
# Wei Shao, Ph.D.
# Advanced Biomedical Computing Science/HIV Replication and Diversity Program
# Frederick National Labatory for Cancer Research
# Frederick, MD, USA
# 02-20-09
#
# Usage: PAPAC_diveristy.pl <Input_blast_result file>
############################################

use strict;
use locale;

my ($blast_result_file, $ref_length, $output_file, $good_seq_line, @mutation_line, $sequence_number, $compared_seq_pairs);
my ($A, $C, $G, $T, $N, $diff_each_position, $total_diff, $deletion, $insertion, $total_nt, $adusted_total_nt);
my ($wt_nt, $average_seq_length, $P_distance, $mutation_each_position, $wt_each_position, %nt_position, $index);
my ($base_number, $site_pairs, $diverstiy_pos, $total_diversity, $total_bases, $pos_with_gaps, $noGap_length);
my ($position); 
my ($JC_distance, $K80_distance, $F84_distance, $TN93_distance);
my ($total_A, $total_C, $total_G, $total_T, $transition, $transition_AG, $transition_CT, $transversion);
my ($Pi_A, $Pi_C, $Pi_G, $Pi_T, $Pi_Y, $Pi_R, $transition_rate, $transition_AG_rate, $transition_CT_rate, $transversion_rate);
my ($ref_position);

$blast_result_file  = $ARGV[0];
$output_file = $blast_result_file."_diveristy.txt";

print "Usage:  PAPAC_diveristy.pl <input blast result file>\n";

%nt_position = ( 'A' => 1, 
                 'C' => 2,
                'G' => 3,
		 'T' => 4,);

$transition = 0;
$transversion = 0;
$transition_CT = 0;
$transition_AG = 0;


open (IN, "$blast_result_file");
	
open (OUT, ">$output_file");

print OUT "Diversity of $blast_result_file\n";	
	
print  "reference\tA\tC\tG\tT\tN\tDeletions\ta\tc\tg\tt\tn\twt base #\tdiff\tsite pairs\tdiveristy\n"; 

my $found_ref_line = 0;
$ref_position = 0;
while (<IN>) {
    chomp;
    if ($_ =~ /Good sequences/) {
	$good_seq_line = $_ ;
	$good_seq_line =~ /(\w+):\s+(\d+)/;
    
	$sequence_number = $2; 
	$compared_seq_pairs = ($sequence_number*($sequence_number - 1) )/2;  ## For Mega diversity calculation.
    }
    $found_ref_line = 1 if ($_ =~ /^reference/);
    if ($found_ref_line) {
	$ref_position++;
    }
    
    next unless ($_ =~ /^\w\d+/);
    
    print "&&&& $_\n";
    $ref_length++;

    next if ($_ =~ /No mutation found/);
   
    @mutation_line = split(/\t/, $_);
 
    if (scalar@mutation_line<2) {

        $diff_each_position = 0;
        $mutation_each_position = 0;
        $mutation_line[1] = 0;
        $mutation_line[2] = 0;
        $mutation_line[3] = 0;
        $mutation_line[4] = 0;
        $mutation_line[5] = 0;
        $mutation_line[6] = 0;

    }

    $wt_nt = $mutation_line[0];
    $wt_nt =~ s/\d+//;
    $wt_nt = uc($wt_nt);
   
    $index = $nt_position{$wt_nt};
 
    $mutation_each_position = $mutation_line[1] + $mutation_line[2] + $mutation_line[3] + $mutation_line[4] + $mutation_line[5] + $mutation_line[6]; ## $mutation_line[6] is deletion. 

    $wt_each_position = $sequence_number - $mutation_each_position;
   
    $mutation_line[$index] = $wt_each_position;

    print "1 @mutation_line\n";
    print "2 $wt_nt = $mutation_line[0]\n";
    print "3 $index = $nt_position{$wt_nt}\n";
   
 print "4 $mutation_line[$index] = $wt_each_position\n";


    $total_A += $mutation_line[1];
    $total_C += $mutation_line[2];
    $total_G += $mutation_line[3];
    $total_T += $mutation_line[4];


    $diff_each_position =  $mutation_line[1]*( $mutation_line[2] + $mutation_line[3] + $mutation_line[4]) ## Not include $mutation_line[5] which is an ambigous N.
                            + $mutation_line[2]*($mutation_line[3] + $mutation_line[4])
                            + $mutation_line[3]*($mutation_line[4]);


    $transition = $transition + ( ($mutation_line[1]*$mutation_line[3]) + ($mutation_line[2]*$mutation_line[4]) );  ## transition = A*G + C*T

    $transition_CT = $transition_CT + ($mutation_line[2]*$mutation_line[4]);

    $transition_AG = $transition_AG + ($mutation_line[1]*$mutation_line[3]);

    $transversion = $transversion + ($mutation_line[1]*$mutation_line[2] + $mutation_line[3]*$mutation_line[4] 
				     +  $mutation_line[2]*$mutation_line[3] +  $mutation_line[1]*$mutation_line[4]); # A*C + G*T + C*G + A*T

    $base_number = $sequence_number - $mutation_line[5] - $mutation_line[6];  ## the sequences without deletion at this site and without N
    
    $site_pairs = ($base_number*($base_number - 1))/2;
    
    next if (  $site_pairs == 0);
    $diverstiy_pos = sprintf("%.6f", ($diff_each_position/$site_pairs));

    $total_diff += $diff_each_position;

    $total_diversity += $diverstiy_pos;

    $deletion += $mutation_line[6];  ### for the method of Pairwise Deletion of gaps only, no need for Complete Deletion of gaps. 
    $total_bases += $base_number; ### for the method of Pairwise Deletion of gaps only, no need for Complete Deletion of gaps.

}


    $Pi_A = $total_A/($ref_length*$sequence_number);  
    $Pi_C = $total_C/($ref_length*$sequence_number);
    $Pi_G = $total_G/($ref_length*$sequence_number );
    $Pi_T = $total_T/($ref_length*$sequence_number);


    $Pi_Y = $Pi_C + $Pi_T;
    $Pi_R = $Pi_A + $Pi_G;


   $transition_AG_rate = $transition_AG/($ref_length*$compared_seq_pairs);
   $transition_CT_rate = $transition_CT/($ref_length*$compared_seq_pairs);

   $transition_rate = $transition/($ref_length*$compared_seq_pairs );  

    $transversion_rate = $transversion/($ref_length*$compared_seq_pairs );


$P_distance = sprintf ("%.6f", ($total_diversity/$ref_length));  

$JC_distance = JCdistance ($P_distance);

$K80_distance = K80_distance ($transition_rate, $transversion_rate ); 

$TN93_distance = TN93_distance ($Pi_A, $Pi_C, $Pi_G, $Pi_T, $Pi_Y, $Pi_R, $transition_AG_rate, $transition_CT_rate, $transversion_rate);

print OUT "P-distance: $P_distance\t JC: $JC_distance\t K80: $K80_distance\t TN93: $TN93_distance\n";
#print OUT "diversity by the method of Pairwise  Deletion of gaps/missing data: $P_distance\n";
print "P-distance: $P_distance\t JC: $JC_distance\t K80: $K80_distance\t TN93: $TN93_distance\n";
#print "diversity by the method of Pairwise  Deletion of gaps/missing data: $P_distance \tMega pairwise deletion diversity\t$Mega_diversity\n";


#################### sub routines ####################

sub JCdistance {
   my  $p_d = shift;
   my ($jc_d);

   $jc_d = sprintf("%.5f", -(3/4)*log(1- (4*$p_d)/3) );  ## log is e based
   
   
   return $jc_d;

}

sub  K80_distance {
    my ($s, $v) = @_;
     my $k80_d;

     $k80_d = sprintf ("%.5f", -(0.5)*log(1- 2*$s - $v) - ((0.25)*log(1 - 2*$v)) );
     
     return $k80_d;
 }


sub TN93_distance {
    my ($pi_a, $pi_c, $pi_g, $pi_t, $pi_y, $pi_r, $s1, $s2, $v) = @_;
    my ($a1, $a2, $b, $tn93_dist);

########## The following numbers are from Zihen Yang's Computational Molecular Evolution for testing purpose ###############
#    $pi_a = 0.3265;
#    $pi_c = 0.2605;
#    $pi_g = 0.1946;
#    $pi_t = 0.2184;
#    $pi_y = $pi_c + $pi_t;
#    $pi_r = $pi_a + $pi_g;
#    $s1 = 0.05591;
#    $s2 = 0.03270;
#    $v = 0.00633;
############################################################################################################################
    
    $a1 = -log(1 - ($pi_y*$s1)/(2*$pi_t*$pi_c) - $v/(2*$pi_y) );
    $a2 = -log(1 - ($pi_r*$s2)/(2*$pi_a*$pi_g) - $v/(2*$pi_r) );
    $b = -log(1 - $v/(2*$pi_y*$pi_r) );
    
    $tn93_dist = sprintf("%.5f", ( (2*$pi_t*$pi_c)/$pi_y)*($a1 - $pi_r*$b) + ((2*$pi_a*$pi_g)/$pi_r)*($a2 - $pi_y*$b) + 2*$pi_y*$pi_r*$b);
 
    return $tn93_dist;

}

