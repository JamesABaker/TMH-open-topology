#DGpred
##Description
Calculate the hydrophobicity (deltaG values) of protein sequences or peptides.

The web-server for dgpred is accessible at http://dgpred.cbr.su.se.

##Usage
The main script for calculating the hydrophobicity is `myscanDG.pl`
This script is modified from the original script `DGpred_multi_scan.pl`
It can be used to calculate DG values by using either single sequence
information or multiple sequence information. In the latter case, the path of 
NCBI blastpgp should be set.

Examples:
In the subfolder "test"

* calculate DG values for the test sequence file 'test.fa' using single sequence information

    `$ ../myscanDG.pl test.fa > test.dgvalue.single.txt`

* calculate DG values for the test sequence file 'test.fa' based on multiple sequence alignment
<pre>
  the blastpgp program should be located at /path/to/blast/bin/blastpgp
  /path/to/blastdb, e.g. /data/blastdb/nr
  /path/to/scampi-msa, is the scampi-msa from the modhmm package
</pre>

    `$ ../myscanDG.pl -multi test.fa -ncbidir /path/to/blast -blastdb /path/to/blastdb -scampipath /path/to/scampi-msa  > test.dgvalue.multi.txt`

* calculate DG values from the prf file 

    `$ ../myscanDG.pl test.prf > test.dgvalue.by.multi.prf.txt`

##Other scripts

* analyze.pl    analyze the delta G values of given a peptide
e.g.

    `$ ../analyze.pl  YIYLGGAILAEVIGTTLMKF -o rst_pep.html`

* calc_dG.pl    calculate the delta G values given a list of peptides
e.g.

    `$ ../calc_dG.pl fragfile.txt -o rst_calc.html`

* dG_predTM.pl  predict the most hydrophobic regions given a protein sequence
e.g.

    `$ ../dG_predTM.pl test.fa -o rst_predTM.html`


