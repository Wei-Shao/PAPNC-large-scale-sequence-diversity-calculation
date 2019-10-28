# PAPNC-large-scale-sequence-diversity-calcution
This tool is for the calculation of large scale sequences generated with next generation sequencing methods. It was created for HIV sequences. However, it can be used for other sequences too. 
-------------------------------------------------------------------
Steps to calculate sequence diversities with PAPNC_diverity.pl:

(1)	Select a reference sequence and name it, for example, HIV1B_shortRT.fasta. The reference sequence must be the same length of  sample sequences (e.g. HIV1B_test.fasta).

(2)	Build a query database with the reference sequence with the name, for example, HIV1B_shortRTDB.
If you use blast makeblastdb:
makeblastdb -in HIV1B_shortRT.fasta -dbtype nucl -out HIV1B_shortRTDB

(3)	Blast the sample sequences to the reference sequence database
For example:
blastn -task blastn -dust no -query HIV1B_test.fasta -db HIV1B_shortRTDB -out HIV1B_test.blast

(4)	Construct mutation table:
blast_to_mutation_table.pl  <HIV1B_shortRT.fasta>  <HIV1B_test.blast>

(5)	Run PAPNC_diversity.pl on the mutation table file to calculate diversity
PAPNC_diversity.pl   <HIV1B_shortRT_mutation_table>


Please cite the following paper if you use these scripts
Shao W, Kearney MF, Boltz VF, Spindler JE, Mellors JW, Maldarelli F, Coffin JM. 2014. PAPNC, a novel method to calculate nucleotide diversity from large scale next generation sequencing data. J Virol Methods, 203:73-80.
