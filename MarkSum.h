
//
// File: MarkSum.h
//
// Purpose: This file defines a summarizer used to build a Markov chain
// model for some class of interest.  This process takes as input some DNA
// sequences and creates a 7th-order Markov chain approximation for the
// sequences.  This summarization may be incorporated into a scoring function
// for use by the sequence analyzer.
//
//	Author: Ron Braun
//
// History:
//		07/30/98	RKB	Created.
//

#ifndef MARKSUMH
#define MARKSUMH

class MARKOV_CHAIN_SUMMARIZER {
   DNA_SEQUENCE_DATA trainingSets[2];  // training sets for the markov
                                       // approximation
public:
   // constructors

   // methods
   SCORING_FUNCTION *PerformMarkovApproximation ( void );
      // performs the 7th-order markov chain approximation on the
      // local training sets.  a scoring function is returned by this
      // procedure.

   // destructors

};

#endif

