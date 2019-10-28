#!/usr/bin/perl
################################################
# blast_to_mutation_table.pl
#
# Wei Shao, Ph.D.
# Advanced Biomedical Computing Science/HIV Replication and Diversity Programm
# Frederick National Laboratory for Cancer Research
# Frederick, Maryland, USA
#
# 10-25-2019
# 
# This script parses  blast files to calculate the numbers of nucleotides at each position of a sample. 
# The reference sequence is used for aligning sample sequences only. Its nucleotide frequencies are not used in PAPNC script.  
# The output mutation table is the input file for PAPNC_diversity.pl
#
# Note: In the blast file, the query is sample sequence while subject is reference sequence. 
# Note: Reversed sequences are filtered here because our experiments do not sequence from the backward direction. 
#
# Usage: blast_to_mutation_table.pl <reference_sequence> <blast_file>
#####################################################

use locale;
use strict;

my ($input_ref, $input_blast, $output_file, @blastblock, $s, $i);
my ($original_ref, $original_ref_length);
my ($blastresult, @alignblock, @query_subj_lines, @refseqline,@sampleline);
my ($rest, $ref_ID, $reflength, $each_query, @each_compare, $score, $identity, @Iline);
my ($Qlength, $Identity, $sampleID, @aligned_lines);
my ($line,@single_refline, @single_sampleline, $first_ref, $last_ref, $first_sample);
my ($last_sample, $ref_seq, $sample_seq, $firstdiff_length);
my (@mutations, $mutations_ref, @all_wt_mutations);
my ($pos_mut, $total_number, $good_sequences, $read_no_error);

$input_ref = $ARGV[0];
$input_blast = $ARGV[1];
$output_file = "mutation_table_".$input_blast.".xls";

$read_no_error = 0;

print " **** blast_to_mutation_table.pl <reference_sequence> <blast_file> ****\n";

open (RIN, "$input_ref") || die "cannot open $input_ref\n";


<RIN>;

$original_ref = <RIN>;

$original_ref_length = length($original_ref);
$total_number = 0;

$/ = "\nQuery=\s+";

open (IN, "$input_blast") || die "cannot open $input_blast\n";;
while (<IN>) {
   
    chomp;
    @blastblock = split(/Query=/, $_);

}

### Now look for mutations in each alignment block.
foreach $blastresult (@blastblock) {
    
  next if ($blastresult =~ /^BLASTN/);
  $total_number++;

    @query_subj_lines = split (/>/, $blastresult);

    @refseqline =();
    @sampleline = ();


    ($sampleID, $rest) = split(/(\s+Length=)/, $query_subj_lines[0]);   ### 0703-2013 for blast+

    $sampleID =~ s/\s+//g;
    $sampleID =~ s/\($b//;



    for ($s = 1; $s < scalar@query_subj_lines; $s++) { ## only one query block.
        $each_query = $query_subj_lines[$s]; 

        @each_compare = split (/Strand=/,$each_query);  ### modified for blastPlus
   
        ($ref_ID, $reflength, $score, $rest, $identity) = split (/\n/, $each_compare[0]);
     
	$reflength =~ s/Length=//g; ### 0703-2013 for blast+
        $reflength =~ s/\s+//g;
        $identity =~ s/Identities =//;
        $identity =~ s/\(/#/g;
        $identity =~ s/%/#/;
        $identity =~ s/\//#/;     
        @Iline = split (/#/, $identity);
        $Qlength = $Iline[0];  ## aligned sample sequence length
        $Identity = $Iline[2]; 
  	       
	@aligned_lines = split(/\n/, $each_compare[1]);
        shift@aligned_lines;
	    
	foreach $line (@aligned_lines) {
	    
	   if ($line =~ /Query/) { ### 0703-2013 for blast+
		
	       $line =~ s/Query\s+//; ### 0703-2013 for blast+
               @single_sampleline = split(/\s+/, $line);
	       @sampleline =(@sampleline, @single_sampleline);
	       @single_sampleline=();
		    
	   }
	   if ($line =~ /Sbjct/) {
	       $line =~ s/Sbjct\s+//; ### 0703-2013 for blast+

               @single_refline = split(/\s+/, $line);
	       @refseqline = (@refseqline,@single_refline); 
	       @single_refline =();		    
	   }
		
	} ## foreach $line (@aligned_lines)
 
    $first_ref = shift(@refseqline);
			
    $last_ref = pop(@refseqline);
    $first_sample = shift(@sampleline);
    $last_sample = pop(@sampleline);
    
    $ref_seq = "@refseqline";
    $sample_seq = "@sampleline";   		
    $ref_seq =~ s/\s+|\d+//g;
    $sample_seq =~ s/\s+|\d+//g;

    $ref_seq = uc($ref_seq);
    $sample_seq = uc($sample_seq);
 
    next if ($last_ref < $first_ref); ## There should be no reversed seq here.	
 #   next if ($original_ref_length - ($last_ref - $first_ref) >20); ## Skip very short ones. e.g. ref20-----------ref149
 #   next if ($original_ref_length - $last_ref >=10); ## if too much truncation, discard it. e.g. ref1-----ref149

    $good_sequences++;

     if (($Qlength == $reflength) && ($Identity == 100)) {

	    $read_no_error++;

	    last;
     }
               
    ## initialization for each sample.

    @mutations = ();  

    $firstdiff_length = $first_ref -1; ## the length of ref fragment that is not aligned.

     $mutations_ref =ScanMutations($ref_seq, $sample_seq, $firstdiff_length);
     @mutations = @{$mutations_ref};    ## mutations for each reads	 

     @all_wt_mutations = (@all_wt_mutations, @mutations); ## if ($ref_ID =~ /wt/);  ## Add mutations from all reads. 
		              
     last; 
   } 
	
} ## foreach $blastresult; An alignment has been processed here. 
   
$pos_mut  = AddMutations (\@all_wt_mutations, $good_sequences);

open (OUT, ">$output_file");

print OUT "Mutations in $input_blast\n\n";

print OUT "Total sequences in this group: $total_number\n";
print OUT "Good sequences in this group: $good_sequences\n";
print OUT "reads identical to reference: $read_no_error\n";
print OUT "Blank lines indicate no mutations found at the sites\n";
print OUT "====================================\n";
print OUT "Compared to wild type reference:\n";

PrintMutations ($original_ref, $pos_mut);


##################################################################
# sub ScanMutations
###################################################################

sub ScanMutations {
   
 my ($ref_seq1, $sample_seq1, $length_diff) = @_;

    my ($i1, $m, @mutations1, $refbase, $samplebase, $ref_length);

    my $insertion1 = 0;


    $ref_length = length($ref_seq1);
  
    @mutations1 = ();
 

     for ($i1=0; $i1 < $ref_length; $i1++) {
         
	    $refbase = substr($ref_seq1, $i1, 1);
	    
	    $samplebase = substr($sample_seq1, $i1, 1);
	    
	    $insertion1++ if ($refbase eq "-"); ## This is for mapping mutations to ref
	   
	    $m = ($i1 + $length_diff)- $insertion1; ## $length_diff is length of truncated $ref_seq1 compared to original ref.
      
          

	    ## $m is used to map a $samplebase to refseq position. Discount all inserts.
          

	    if ($samplebase  ne $refbase){             	
              
		if ($refbase eq "-") {
		    push @mutations1, $m."+".lc($samplebase); ## This is insertion.
		   
		}
		else { 

		     push @mutations1, $m."+".$samplebase;   ### finds mutation here.

		}
	    }
            
	    $refbase ="";
	    $samplebase ="";
        
	}
 return (\@mutations1);
   
}


##################################################################################################
# sub AddMutations
# This subroutin converts arrays into hashes so that mutations at a postion is aligned together
##################################################################################################

sub AddMutations {
    my ($pos_mut_arrayref, $good_seq) = @_;

    my ($pos, $mutation, @pos_mut, %pos_mut);
    my ($mFS, $pFS, $A, $C, $G, $T, $i2, $del, $nt, $N);
    my ($a, $c, $g, $t, $n);  ## insertions
    my ($point_mut_freq, $indel_mut_freq);
    
    @pos_mut = @{$pos_mut_arrayref};
    %pos_mut =();
    $mFS = 0;
    $pFS = 0;
    
   foreach $pos_mut(sort @pos_mut) {
       
       ($pos, $mutation) = split (/\+/, $pos_mut);
       $mFS++ if ($mutation eq "-");
      
       $pFS++ if ($mutation =~ /[acgt]/);

       if ($pos_mut{$pos}) {
          $pos_mut{$pos} = $pos_mut{$pos}.$mutation;
     	   
   
       }
       else {
          
          $pos_mut{$pos} = $mutation;
       }
 
    } 
   
    foreach $pos (keys %pos_mut) {  
	$A = 0;
	$C = 0;
	$G = 0;
	$T = 0;
        $N = 0;
	$del = 0;
        $a = 0;
        $c = 0;
        $g = 0;
        $t = 0;
        $n = 0;
        for ($i2 = 0; $i2 <length$pos_mut{$pos}; $i2++) {
            $nt = substr($pos_mut{$pos}, $i2,1);
            $A++ if ($nt eq "A");
            $C++ if ($nt eq "C");
            $G++ if ($nt eq "G");
            $T++ if ($nt eq "T");
            $N++ if ($nt eq "N");
            $del++ if ($nt eq "-");
            $a++ if ($nt eq "a");
	    $c++ if ($nt eq "c");
	    $g++ if ($nt eq "g");
	    $t++ if ($nt eq "t");
            $n++ if ($nt eq "n");
	    
        }
        
	$point_mut_freq = sprintf ("%.2f", ( ($A + $C + $G + $T +$N)/$good_seq) *100);
        $indel_mut_freq = sprintf ("%.2f", ( ($del + $a + $c + $g+ $t + $n)/$good_seq) *100);
        $pos_mut{$pos} =$A."\t".$C."\t".$G."\t".$T."\t".$N."\t".$del."\t".$a."\t".$c."\t".$g."\t".$t."\t".$n. "\t".$point_mut_freq. "\t".$indel_mut_freq;

    }
    
    return (\%pos_mut);
}


#######################
# sub PrintMutations
######################

sub PrintMutations { 

my ($ref, $posmut) = @_;
 
    my ($j, $c, %posmut);



    %posmut = %{$posmut};
       
    print OUT "\t\tMutations\t\t\t\tDeletions\t\tInsertions\n";
    print OUT "reference\tA\tC\tG\tT\tN\tDeletions\ta\tc\tg\tt\tn\tpoint mutation %\tindel %\n";
   
   for ($j=0; $j<(length$ref)-1; $j++) {
     if ($posmut{$j}) { 
         print OUT substr($ref, $j, 1), $j+1, "\t$posmut{$j}\n";

     }
     else {
         print OUT substr($ref, $j,1),  $j+1, "\n";
	
     }
   } 

    print OUT "********************************\n";
}
