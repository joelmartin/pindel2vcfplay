# pindel2vcfplay

pindel2vcf for plant 30,000+ contig genomes, and normal genomes but most important for those

on a rice( kitaake v3 ) with just 58 contigs, speedup is
3300 seconds ( pindel2vcf 0.6.3 )
340  seconds ( pindel2vcf 0.6.0 - before the slowdowns were added )
128 seconds  ( this version )


changes are:
  revert to getline from the get loop
  short-cut out of the more intense 'is this a summary line' check if the line can't be one
  use fasta index file (.fai) ( pindel already requires one ) instead of walking through whole fasta.
 
 output identical, 95% of code identical but moved around so I could follow it.
 some cruft added to get xcode to shut up and because I'm terrible.
