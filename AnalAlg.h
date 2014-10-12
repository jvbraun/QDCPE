
//////////////////////////////////////////////////////////////////////////////
//
// File: AnalAlg.h
//
// Purpose: This file defines analysis algorithms to be applied to some
//    sequence for segmentation.  There may be two flavors of a given
//    algorithm, a general definition and a FAST definition.  The FAST
//    version works for small sequences by building a scoring matrix for
//    use by the algorithm.  Sequences that are too large for a scoring
//    matrix (which is of a size length**2) must use the general algorithm.
//
//	Author: Ron Braun
//
// History:
//		07/30/98	RKB	Created.
//
//////////////////////////////////////////////////////////////////////////////

#ifndef ANALALGH
#define ANALALGH

#include "AnalRes.h"
#include "DNASeqDt.h"
#include "Globals.h"
#include "Matrix.h"


//////////////////////////////
// class ANALYSIS_ALGORITHM //
//////////////////////////////

// general class to organize all analysis algorithms.

class ANALYSIS_ALGORITHM {
   // annotate whether we are responsible for deallocating
   // the results object
   boolean deallocateResults;

protected:
   // sequence to be analyzed by algorithm
   DNA_SEQUENCE_DATA *sequence;

public:
   // results of analysis
   ANALYSIS_RESULTS *analysisResults;

   /////////////////
   // constructors
   //

   // create new general algorithm, building a new
   // analysis results object in the process.  results
   // object will be deallocated when this object is
   // deleted.
   ANALYSIS_ALGORITHM ( DNA_SEQUENCE_DATA *seq,
                        long numCPs ) {
      // save sequence to be analyzed
      sequence = seq;

      // create new results object for algorithm
      analysisResults = new ANALYSIS_RESULTS ( numCPs );

      // we are responsible for deallocating this object
      deallocateResults = TRUE;

   };

   // create a new general algorithm, using a given
   // analysis results object to store the results.
   // caller has responsibility to deallocate results
   // object.
   ANALYSIS_ALGORITHM ( DNA_SEQUENCE_DATA *seq,
                        ANALYSIS_RESULTS *analRes ) {
      // save sequence to be analyzed
      sequence = seq;

      // save pointer to results object to be used
      analysisResults = analRes;

      // caller must deallocate results object
      deallocateResults = FALSE;

   };

   ////////////
   // methods
   //

   // permit reuse of the existing results object if desired.  use of
   // this option requires caller to assume responsibility for the
   // deallocation of the results object.
   void UseNewResults ( ANALYSIS_RESULTS *analRes );

   // analysis function, specified for each algorithm class.  this
   // implements the actual algorithm to perform segmentation.
   virtual void Analyze ( long minSegmentSize ) { };

   ////////////////
   // destructors
   //

   // delete allocated results arrays
   ~ANALYSIS_ALGORITHM ( void ) {
      if ( deallocateResults )
         delete analysisResults;

   };

};


////////////////////////
// class SCORE_MATRIX //
////////////////////////

// this class implements a matrix which stores scoring values for all
// subsequences of the sequence to be analyzed.  this support efficient
// algorithm development, but because the matrix is of size length**2,
// only shorter sequences can be analyzed with this technique.

class SCORE_MATRIX {
   // annotates whether or not original score matrix has been deallocated
   boolean deallocatedSM;

protected:
   // matrix of segment scores.  a whole matrix is defined to favor
   // speed over space, though only half of it is used.  the data
   // is populated at matrix creation.
   MATRIX<double> *scoreMatrix;

   // length of sequence which generated this scoring function
   long length;

public:
   /////////////////
   // constructors
   //

   // create an empty segment matrix given a sequence of some length.
   // the data is not populated yet!  it should be populated with a
   // given sequence of the same length before the scoring matrix is
   // used computationally.
   SCORE_MATRIX ( long seqLength );

   ////////////
   // methods
   //

   // access an element of the segment matrix
   double &operator () ( int x, int y ) {
      return scoreMatrix->Data ( x, y );
   };

   // generate the segment scores using the quasi-deviance method.  this
   // should probably be generalized for other methodologies, when it is
   // obvious how to do so.
   void ScoreQuasiDeviance ( DNA_SEQUENCE_DATA *Y );

   // use someone else's score matrix in place of our own.  warning: if
   // the score matrix passed in here becomes deallocated by its owner
   // then we lose use of the matrix too!
   void UseAlternateScoreMatrix ( MATRIX<double> *alternateSM );

   // fetch this score matrix
   MATRIX<double> *GetScoreMatrix ( void ) {
      return scoreMatrix;

   };

   ////////////////
   // destructors
   //

   // delete matrix
   ~SCORE_MATRIX ( void );

};


/////////////////////////////////////////////////////
// subclass BRUTE_FORCE_FAST of ANALYSIS_ALGORITHM //
/////////////////////////////////////////////////////

// this algorithm performs brute-force search for optimized segmentation
// as per original fortran code.  it has been generalized to work for
// an arbitrary number of change points using recursion, rather than the
// hard-coded fortran solution.  though obsoleted by Auger-Lawrence, it is
// retained as a benchmark.  see J. Braun's thesis research for context
// of original algorithm.

class BRUTE_FORCE_FAST : public ANALYSIS_ALGORITHM,
                         public SCORE_MATRIX {
   //////////////////////
   // utility functions
   //

   // initialize the loop iterators and kick of the recursion.  used to
   // implement the nested for-structures of the original fortran in
   // a generalized manner.
   void FindChangePoints ( int numCPs, long minSegmentSize );

   // this recursive function iterates down another level into an arbitrary
   // length nested for-structure.  this generalizes the hardcoded loops of
   // the original fortran.
   void ForIterate ( int recurseStep, int *loopIterators,
                     int numCPs, long minSegmentSize );

public:
   /////////////////
   // constructors
   //

   // create instance of brute force algorithm, creating local
   // results object
   BRUTE_FORCE_FAST ( DNA_SEQUENCE_DATA *seq, long numCPs ) :
         ANALYSIS_ALGORITHM ( seq, numCPs ),
         SCORE_MATRIX ( seq->length ) {
      // class constructors completely handle construction

   };

   // create instance of brute force algorithm, using given
   // results object
   BRUTE_FORCE_FAST ( DNA_SEQUENCE_DATA *seq, ANALYSIS_RESULTS *analRes ) :
         ANALYSIS_ALGORITHM ( seq, analRes ),
         SCORE_MATRIX ( seq->length ) {
      // class constructors completely handle construction

   };

   ////////////
   // methods
   //

   // perform analysis of segment matrix
   void Analyze ( long minSegmentSize );

   ////////////////
   // destructors
   //

};


////////////////////////////////////////////////////////
// subclass AUGER_LAWRENCE_FAST of ANALYSIS_ALGORITHM //
////////////////////////////////////////////////////////

// see Auger-Lawrence paper entitled "Algorithms for the Optimal
// Identification of Segment Neighborhoods" for the specification of
// this algorithm.

class AUGER_LAWRENCE_FAST : public ANALYSIS_ALGORITHM,
                            public SCORE_MATRIX {

public:
   /////////////////
   // constructors
   //

   // create instance of auger-lawrence algorithm, creating local
   // results object
   AUGER_LAWRENCE_FAST ( DNA_SEQUENCE_DATA *seq, long numCPs ) :
         ANALYSIS_ALGORITHM ( seq, numCPs ),
         SCORE_MATRIX ( seq->length ) {
      // class constructors completely handle construction

   };

   // create instance of auger-lawrence algorithm, using given
   // results object
   AUGER_LAWRENCE_FAST ( DNA_SEQUENCE_DATA *seq, ANALYSIS_RESULTS *analRes ) :
         ANALYSIS_ALGORITHM ( seq, analRes ),
         SCORE_MATRIX ( seq->length ) {
      // class constructors completely handle construction

   };

   ////////////
   // methods
   //

   // perform analysis of segment matrix
   void Analyze ( long minSegmentSize );

   ////////////////
   // destructors
   //

};


///////////////////////////////////////////////////
// subclass AUGER_LAWRENCE of ANALYSIS_ALGORITHM //
///////////////////////////////////////////////////

// see Auger-Lawrence paper entitled "Algorithms for the Optimal
// Identification of Segment Neighborhoods" for the specification of
// this algorithm.

class AUGER_LAWRENCE : public ANALYSIS_ALGORITHM {

public:
   /////////////////
   // constructors
   //

   // create instance of auger-lawrence algorithm, creating local
   // results object
   AUGER_LAWRENCE ( DNA_SEQUENCE_DATA *seq, long numCPs ) :
         ANALYSIS_ALGORITHM ( seq, numCPs ) {
      // class constructor completely handles construction

   };

   // create instance of auger-lawrence algorithm, using given
   // results object
   AUGER_LAWRENCE ( DNA_SEQUENCE_DATA *seq, ANALYSIS_RESULTS *analRes ) :
      ANALYSIS_ALGORITHM ( seq, analRes ) {
      // class constructor completely handles construction

   };

   ////////////
   // methods
   //

   // perform analysis of segment matrix
   void Analyze ( long minSegmentSize );

   ////////////////
   // destructors
   //

};


/////////////////////////////////////////////////////
// subclass MULTI_START_FAST of ANALYSIS_ALGORITHM //
/////////////////////////////////////////////////////

// see Braun, Braun, and Mueller paper entitled "Multiple Change-Point
// Fitting Via Quasi-Likelihood, With Applications To DNA Sequence
// Segmentation" for specification of this algorithm.

class MULTI_START_FAST : public ANALYSIS_ALGORITHM,
                         public SCORE_MATRIX {
   // number of random CP configurations to try in minimizing
   // the desired metric
   int numSamples;

   //////////////////////
   // utility functions
   //

   // generates a random segmentation of a sequence.  the segmentation
   // must obey minimal segment size constraints.
   void GenerateNewGrid ( int *grid, int gridSize,
                          RANDOM_NUMBER_GENERATOR &rndGen,
                          long minSegmentSize );

   // performs the segmentation analysis for a single number of change
   // points.  the scoring matrix needs to have been computed prior to
   // this point!
   void AnalyzeSingleCP ( long minSegmentSize, int numCPs );

public:
   /////////////////
   // constructors
   //

   // create instance of multi-start algorithm, creating local
   // results object
   MULTI_START_FAST ( DNA_SEQUENCE_DATA *seq, long numCPs, int numReps ) :
         ANALYSIS_ALGORITHM ( seq, numCPs ),
         SCORE_MATRIX ( seq->length ) {
      // annotate number of sample repititions used by algorithm
      numSamples = numReps;

   };

   // create instance of multi-start algorithm, using given
   // results object
   MULTI_START_FAST ( DNA_SEQUENCE_DATA *seq, ANALYSIS_RESULTS *analRes,
                      int numReps ) :
         ANALYSIS_ALGORITHM ( seq, analRes ),
         SCORE_MATRIX ( seq->length ) {
      // annotate number of sample repititions used by algorithm
      numSamples = numReps;

   };

   ////////////
   // methods
   //

   // perform complete analysis of sequence using scoring matrix
   void Analyze ( long minSegmentSize );

   // perform an analysis of sequence using multi-start method for only
   // one change point number
   void AnalyzeCP ( long minSegmentSize, int cpNum, boolean suppressScore );

   ////////////////
   // destructors
   //

};


////////////////////////////////////////////////
// subclass MULTI_START of ANALYSIS_ALGORITHM //
////////////////////////////////////////////////

// see Braun, Braun, and Mueller paper entitled "Multiple Change-Point
// Fitting Via Quasi-Likelihood, With Applications To DNA Sequence
// Segmentation" for specification of this algorithm.

class MULTI_START : public ANALYSIS_ALGORITHM {
   // number of random CP configurations to try in minimizing
   // the desired metric
   int numSamples;

   //////////////////////
   // utility functions
   //

   // generates a random segmentation of a sequence
   void GenerateNewGrid ( int *grid, int gridSize,
                          RANDOM_NUMBER_GENERATOR &rndGen,
                          long minSegmentSize );

   // performs the segmentation analysis for a single number of change
   // points
   void AnalyzeSingleCP ( long minSegmentSize, int numCPs );

public:
   /////////////////
   // constructors
   //

   // create instance of multi-start algorithm, creating local
   // results object
   MULTI_START ( DNA_SEQUENCE_DATA *seq, long numCPs, int numReps ) :
         ANALYSIS_ALGORITHM ( seq, numCPs ) {
      // annotate number of sample repititions used by algorithm
      numSamples = numReps;

   };

   // create instance of multi-start algorithm, using given
   // results object
   MULTI_START ( DNA_SEQUENCE_DATA *seq, ANALYSIS_RESULTS *analRes,
                 int numReps ) :
         ANALYSIS_ALGORITHM ( seq, analRes ) {
      // annotate number of sample repititions used by algorithm
      numSamples = numReps;

   };

   ////////////
   // methods
   //

   // perform complete analysis of sequence
   void Analyze ( long minSegmentSize );

   // perform an analysis of sequence using multi-start method for only
   // one change point number
   void AnalyzeCP ( long minSegmentSize, int cpNum );

   ////////////////
   // destructors
   //
};

//////////////////////////////////////////////////
// subclass ORDERED_SPLIT of ANALYSIS_ALGORITHM //
//////////////////////////////////////////////////

// see e-mail from jerry for definition of this algorithm.

class ORDERED_SPLIT : public ANALYSIS_ALGORITHM {

public:
   /////////////////
   // constructors
   //

   // create instance of ordered-split algorithm, creating local
   // results object
   ORDERED_SPLIT ( DNA_SEQUENCE_DATA *seq, long numCPs ) :
         ANALYSIS_ALGORITHM ( seq, numCPs ) {
      // class constructor completely handles construction

   };

   ////////////
   // methods
   //

   // perform analysis of segment matrix
   void Analyze ( long minSegmentSize );

   ////////////////
   // destructors
   //

private:
   struct TRACKING_NODE {
      long firstElem;   // first element of subsequence
      long lastElem;    // last element of subsequence
      long cpElem;      // best change point for this subsequence
      double qdDelta;   // change in qd if best cp is used
   };

   // utility functions
   void findBestChangePointAndQdDelta(TRACKING_NODE &node,
                                      double currentQd);


};

#endif

