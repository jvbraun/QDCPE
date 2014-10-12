
//
// File: SeqAnal.h
//
// Purpose: This file defines the sequence analyzer process for analysis
//    of a given FASTA sequence file.
//
//	Author: Ron Braun
//
// History:
//		07/30/98	RKB	Created.
//

#ifndef SEQANALH
#define SEQANALH

/////////////////////////////
// class SEQUENCE_ANALYZER //
/////////////////////////////

// this object performs the analysis of a FASTA sequence file

// TBD this is done crudely for the moment, until it is clearer how
// the analyzer function will be used.  it may be useful to let the
// simulators define sequence analyzers for their use.  for the moment,
// it seems more efficient to let the sims handle their own analyses,
// since they do some optimizations that might not be generally
// supportable through this mechanism.

class SEQUENCE_ANALYZER {
   // filename for fasta dna sequence
   // TBD eventually pull off a local copy of this instead of relying
   // on the outside context to maintain the string allocation.
   char *sequenceFileName;

   // max number of change points
   int maxCPs;

   // minimum segment length
   long minSegmentLength;

   // algorithm used to analyze sequence
   ALGORITHMS algorithmType;

   // power factor for computing C_N
   double cnPower;

public:
   /////////////////
   // constructors

   SEQUENCE_ANALYZER ( char *fName, int numCPs, long minLen,
                       ALGORITHMS algType, double cn ) {
      // save state
      sequenceFileName = fName;
      maxCPs = numCPs;
      minSegmentLength = minLen;
      algorithmType = algType;
      cnPower = cn;

   };

   ////////////
   // methods

   // perform the analysis of the sequence
   void AnalyzeSequence ( void );

};

#endif

