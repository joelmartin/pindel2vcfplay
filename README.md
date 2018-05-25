# pindel2vcfplay

pindel2vcf for messy genomes with lots of contigs/scaffolds, and normal genomes too.

timing;  
1. kitaake - 12 chromosomes followed by 1300 scaffolds 389 mb
*  0.6.3 - 56 minutes   
*  0.6.0 - 5 minutes   
*  this - 30 seconds  

2. panicum - 9 chromosomes followed by 8400 scaffolds  537 mb
   result files pre-grepped for ChrID lines ( speedup larger for non-pre-grepped )
*   0.6.3  killed after 2 days, may try again for completeness
*   0.6.0  22 hours 46 minutes
*   this   41 minutes

3. clostridium - 1 contig, 3.5 mb
* 0.6.3 2 seconds
* 0.6.0 1 second
* this  1 second

changes are:  
revert to getline from the get loop  
short-cut out of the more intense 'is this a summary line' check if the line can't be one  
use fasta index file (.fai) ( pindel already requires one ) instead of walking through whole fasta.  
index the pindel output during first pass over it, avoids re-reading whole thing for every new contig
 
 output identical, 95% of code identical but moved around so I could follow it.
 some cruft added to get xcode to shut up and because I'm terrible.
