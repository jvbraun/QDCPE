
//
// File: DNASeqDt.h
//
// Purpose: This file defines the interface to a DNA sequence data file.
// This file must be in the FASTA format.  This object encapsulates the data
// file and provides a consistent interface for the analyzer and other
// interested components.
//
//	Author: Ron Braun
//
// History:
//		07/30/98	RKB	Created.
//

#ifndef DNASEQDTH
#define DNASEQDTH

#include <stdio.h>

#include "GenFile.h"
#include "Globals.h"
#include "Matrix.h"
#include "Random.h"


/////////////////////////////
// class DNA_SEQUENCE_DATA //
/////////////////////////////

// this general class represents a DNA sequence of a given length

// dna sequence data
class DNA_SEQUENCE_DATA {
public:
   // length of sequence
   long length;
   // sequence data.  zeroeth element is not used, to conform to
   // indexing inherited from original fortran code.
   long *data;
   // number of elements in sequence alphabet
   long alphaLength;

   /////////////////
	// constructors

   // create sequence of length N.  sequence data is not populated!
   DNA_SEQUENCE_DATA ( long N );

   ////////////
   // methods

   // returns quasi-deviance of subsequence (l1, l2] of this dna sequence.
   // if l1 and l2 are not specified, the quasi-deviance of the entire
   // sequence is returned.
   double QuasiDeviance ( long l1 = 0, long l2 = 0 );

   // computes the mean of subsequence (l1, l2] of this dna sequence and
   // stores it in the mu parameter.  if l1 and l2 are not specified, the
   // mean of the entire sequence is computed and stored.
   void ComputeMean ( double mu[5], long l1 = 0, long l2 = 0 );

   ////////////////
   // destructors

   // destroy this sequence
   ~DNA_SEQUENCE_DATA ( void );

};


////////////////////////////////////////////////////////
// subclass MULTINOMIAL_SEQUENCE of DNA_SEQUENCE_DATA //
////////////////////////////////////////////////////////

class MULTINOMIAL_SEQUENCE : public DNA_SEQUENCE_DATA {
   // max permitted number of change points
   long rMax;

   // utility functions
   long mlndev ( double p[] );
   void CreateNewMu0 ( double jumpSize );
   void CreateNewTau0 ( long rate,
                        double minRatio );

public:
   // actual number of change points for sequence
   long rTrue;
   // mean function characterizing multinomial data
   MATRIX<double> *mu0;
   // change points for multinomial data
   double *tau0;

   /////////////////
   // constructors

   // create multinomial sequence using length and max permitted number
   // of change points.  the sequence is not populated by the constructor.
   MULTINOMIAL_SEQUENCE ( long len, long numMax ) :
      DNA_SEQUENCE_DATA ( len ) {
      // save max number of change points
      rMax = numMax;

      // create muo and tau characterization matrices
      mu0 = new MATRIX<double> ( numMax+1, alphaLength );
      tau0 = new double[numMax+2];

   };

   ////////////
   // methods

   // takes the mean step-function defined by TAU0 and MU0 and generates
   // a vector of multinomial variates following those means.  this vector
   // is used to populate the sequence data.
   void GenerateNewSequence ( long rate,
                              double jumpSize,
                              long minSegmentLength );

   ////////////////
   // destructors

   // destroy mu0 and tau0 matrices
   ~MULTINOMIAL_SEQUENCE ( void ) {
      delete mu0;
      delete[] tau0;
   };

};


/////////////////////////////////////////////////////
// subclass DATAFILE_SEQUENCE of DNA_SEQUENCE_DATA //
/////////////////////////////////////////////////////

class DATAFILE_SEQUENCE : public DNA_SEQUENCE_DATA {
   FASTA_SEQUENCE_FILE *sequenceFile;

public:
   /////////////////
   // constructors

   // create datafile sequence from a given datafile object
   DATAFILE_SEQUENCE ( FASTA_SEQUENCE_FILE *seqFile );

   ////////////
   // methods

   // converts the fasta format sequence into internal sequence
   // representation
   void GenerateNewSequence ( void );

   ////////////////
   // destructors

};

#endif

