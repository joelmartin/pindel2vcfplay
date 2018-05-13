#ifndef PINDEL2VCF_H
#define PINDEL2VCF_H

#include <string>

/* the global parameter g_par stores the values of all parameters set by the user; the Parameter class, in contrast, is for user-friendly IO. */
struct ParameterSettings {
   std::string reference;
   std::string referenceName;
   std::string referenceDate;
   std::string pindelfile;
   std::string pindelroot;
   std::string vcffile;
   std::string chromosome;
   int windowSize;
   int MinCoverage;

   double HetCutoff;
   double HomCutoff;
   int minsize;
   int maxsize;
   bool bothstrands;
   int minsuppSamples;
   int minsuppReads;
   int maxSuppReads;
   int regionStart;
   int regionEnd;
   int maxInterRepeatNo;
   int maxInterRepeatLength;
   int maxPostRepeatNo;
   int maxPostRepeatLength;
   bool onlyBalancedSamples;
   int minimumStrandSupport;
   int compactOutput;
   bool showHelp;
   bool somatic;
   bool gatkCompatible;
} g_par;

void createComplement( const std::string& dna, std::string& complement );

#endif

