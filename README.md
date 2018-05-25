# pindel2vcfplay

pindel2vcf for messy genomes with lots of contigs/scaffolds, and normal genomes too.

timing;
  kitaake - 12 chromosomes followed by 1300 scaffolds
   v 0.6.3  56 minutes
   v 0.6.0   5 minutes
   v this   30 seconds

  panicum - 9 chromosomes followed by 8400 scaffolds
   result files pre-grepped for ChrID lines ( speedup more extreme for non-pre-grepped )
   v 0.6.3  killed after 2 days, may try again for completeness
   v 0.6.0  22 hours 46 minutes
   v this   41 minutes

  clostridium - 1 contig, 3.5mb
   v 0.6.3 2 seconds
   v 0.6.0 1 second
   v this  1 second

changes are:
  revert to getline from the get loop
  short-cut out of the more intense 'is this a summary line' check if the line can't be one
  use fasta index file (.fai) ( pindel already requires one ) instead of walking through whole fasta.
 
 output identical, 95% of code identical but moved around so I could follow it.
 some cruft added to get xcode to shut up and because I'm terrible.
