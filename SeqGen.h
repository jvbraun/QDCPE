
//
// File: SeqGen.h
//
// Purpose: This file defines a generator to create test sequences for
//    characterizing the finite-sample properties of different analysis
//    functions.
//
//	Author: Ron Braun
//
// History:
//		07/31/98	RKB	Created.
//

#ifndef DNASEQDTH
#define DNASEQDTH

#include "DNASeqDt.h"

class SEQUENCE_GENERATOR {
public:
   // constructors

   // methods
   DNA_SEQUENCE_DATA *GenerateTestSequence ( int r );
      // generates a test sequence

   // destructors

};

#endif

