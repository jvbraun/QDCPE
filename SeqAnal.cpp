
//
// File: SeqAnal.cpp
//
// Purpose: This file defines the sequence analyzer process for analysis
//    of a given FASTA sequence file.
//
//	Author: Ron Braun
//
// History:
//		10/25/98 RKB   Created.
//

#include <math.h>
#include <stdlib.h>

#include "AnalAlg.h"
#include "SeqAnal.h"

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

/////////////////
// constructors

// defined inline

////////////
// methods

// perform the analysis of the sequence
void SEQUENCE_ANALYZER :: AnalyzeSequence ( void ) {
   // algorithm used to analyze sequence
   ANALYSIS_ALGORITHM *analAlg;

   // interface to desired FASTA sequence
   FASTA_SEQUENCE_FILE seqFile ( sequenceFileName );

   if ( seqFile.OpenDNASequenceFile ( ) == FALSE ) {
      fflush ( stdin );
      getchar ( );
      exit ( 1 );
   };

   // create sequence from bacteriophage datafile
   DATAFILE_SEQUENCE sequence ( &seqFile );

   // create analysis algorithm to be used on sequence
   switch ( algorithmType ) {
      case BRUTE_FORCE_FAST_ALGORITHM :
         analAlg = new BRUTE_FORCE_FAST ( &sequence, maxCPs );
      break;

      case AUGER_LAWRENCE_FAST_ALGORITHM :
         analAlg = new AUGER_LAWRENCE_FAST ( &sequence, maxCPs );
      break;

      case AUGER_LAWRENCE_ALGORITHM :
         analAlg = new AUGER_LAWRENCE ( &sequence, maxCPs );
      break;

      case MULTI_START_FAST_ALGORITHM :
         analAlg = new MULTI_START_FAST ( &sequence, maxCPs, 100 );
      break;

      case MULTI_START_ALGORITHM :
         analAlg = new MULTI_START_FAST ( &sequence, maxCPs, 100 );
      break;

   case ORDERED_SPLIT_ALGORITHM :
     analAlg = new ORDERED_SPLIT ( &sequence, maxCPs );
   break;

      default :
         printf ( "No algorithm specified!\n" );
         fflush ( stdin );
         getchar ( );
         exit ( 1 );

   };

   // perform analysis
   analAlg->Analyze ( minSegmentLength );

   printf ( "Results of datafile sequence analysis:\n" );
   analAlg->analysisResults->Print ( );

   seqFile.CloseDNASequenceFile ( );

   // delete dynamics
   delete analAlg;

};

