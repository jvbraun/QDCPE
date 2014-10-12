
//////////////////////////////////////////////////////////////////////////////
//
// File: AnalAlg.cpp
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
//		09/21/98	RKB	Created.
//
//////////////////////////////////////////////////////////////////////////////

#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// TBD
#include <time.h>

#include "AnalAlg.h"
#include "ScoreFun.h"


//////////////////////////////
// class ANALYSIS_ALGORITHM //
//////////////////////////////

// general class to organize all sequence analysis algorithms.

/////////////////
// constructors
//

// defined inline in spec

////////////
// methods
//

// permit reuse of the existing results object if desired.  use of
// this option requires caller to assume responsibility for the
// deallocation of the results object.
void ANALYSIS_ALGORITHM :: UseNewResults ( ANALYSIS_RESULTS *analRes ) {
   // if we are responsible for old results, free them up
   if ( deallocateResults ) {
      delete analysisResults;

      // we are no longer responsible for the results
      deallocateResults = FALSE;

   };

   // save new results object
   analysisResults = analRes;

};

////////////////
// destructors
//

// defined inline in spec


////////////////////////
// class SCORE_MATRIX //
////////////////////////

// this class implements a matrix which stores scoring values for all
// subsequences of the sequence to be analyzed.  this support efficient
// algorithm development, but because the matrix is of size length**2,
// only shorter sequences can be analyzed with this technique.

/////////////////
// constructors
//

// create an empty segment matrix given a sequence Y.  the data is
// not populated yet!
SCORE_MATRIX :: SCORE_MATRIX ( long seqLength ) {
   // annotate length
   length = seqLength;

   // allocate new matrix of necessary size
   scoreMatrix = new MATRIX<double> ( length, length );

   // matrix is initially allocated
   deallocatedSM = FALSE;

};

////////////
// methods
//

// generate the segment scores using the quasi-deviance method.  this
// should probably be generalized for other methodologies, when it is
// obvious how to do so.
void SCORE_MATRIX :: ScoreQuasiDeviance ( DNA_SEQUENCE_DATA *Y ) {

   // sanity check to ensure score matrix is compatible with sequence
   if ( length != Y->length ) {
      printf ( "Error!  Scoring matrix and sequence size mismatch!\n" );

      fflush ( stdin );
      getchar ( );
      exit ( 1 );

   };

   double mu[5];
   int muCount[5];
   int numElems;
   double q;

   // populate data by computing segment scores.  the start of the
   // segment is l1+1, so that l1 represents the change-point estimator.

   // build column by column
   for ( int column = length; column >= 1; column-- ) {
      // initialize column state
      numElems = 0;
      for ( int k = 1; k <= 4; k++ ) {
         mu[k] = 0.0;
         muCount[k] = 0;
      };

      for ( int row = column; row >= 1; row-- ) {
         // update state with new element added
         mu[Y->data[row]]++;
         muCount[Y->data[row]]++;
         numElems++;

         // recompute quasi-deviance with new state
         q = 0.0;
         for ( int k = 1; k <= 4; k++ ) {
            if ( mu[k] ) q += muCount[k] * log ( mu[k] / numElems );
         };

         scoreMatrix->Data ( row, column ) = -2 * q;

      };

   };

};

// use someone else's score matrix in place of our own.  warning: if
// the score matrix passed in here becomes deallocated by its owner
// then we lose use of the matrix too!
void SCORE_MATRIX :: UseAlternateScoreMatrix ( MATRIX<double> *alternateSM ) {

   // if original score matrix has not been deallocated, do so
   if ( ! deallocatedSM ) {
//JVB      free ( scoreMatrix );
         delete scoreMatrix;

      // annotate that original matrix has been deallocated
      deallocatedSM = TRUE;

   };

   // annotate new matrix
   scoreMatrix = alternateSM;

};

////////////////
// destructors
//

// delete matrix

SCORE_MATRIX :: ~SCORE_MATRIX ( void ) {
   // if the score matrix has not been already deallocated,
   // do so now
   if ( ! deallocatedSM )
//JVB      free ( scoreMatrix );
        delete scoreMatrix;

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

//////////////////////
// utility functions
//

// this recursive function iterates down another level into an arbitrary
// length nested for-structure.  this generalizes the hardcoded loops of
// the original fortran.
void BRUTE_FORCE_FAST :: ForIterate ( int recurseStep, int *loopIterators,
                                      int numCPs, long minSegmentSize ) {

   // convenience aliases
   CHANGE_POINTS &changePoints = *(analysisResults->changePointMatrix);
   double *fitMetric = analysisResults->fitMetric;

   // loop index for this level of nesting
   int loopIndex = numCPs - recurseStep + 1;

   // iterate on the curStep level of the nested loop
   for ( int segPos = loopIterators[loopIndex-1] + minSegmentSize;
         segPos <= length - recurseStep * minSegmentSize;
         segPos++ ) {
      // annotate this nested loop index
      loopIterators[loopIndex] = segPos;

      // recursive base case
      if ( recurseStep != 1 ) {
         // recurse into next nested loop
         ForIterate ( recurseStep-1, loopIterators, numCPs, minSegmentSize );

      } else {
         // compute final metric value
         double q = 0;

         for ( int cpNum = 1; cpNum <= numCPs; cpNum++ )
            q += (*scoreMatrix) ( loopIterators[cpNum-1] + 1,
                                  loopIterators[cpNum] );

         q += (*scoreMatrix) ( loopIterators[numCPs] + 1, length );

         // see if we have a new minimum
         if ( q <= fitMetric[numCPs] ) {
            // new minimum, save it
            fitMetric[numCPs] = q;

            // update change points
            for ( int cpNum2 = 1; cpNum2 <= numCPs; cpNum2++ )
               changePoints ( cpNum2, numCPs ) = loopIterators[cpNum2];

         };

      };

   };

};

// initialize the loop iterators and kick of the recursion.  used to
// implement the nested for-structures of the original fortran in
// a generalized manner.
void BRUTE_FORCE_FAST :: FindChangePoints ( int numCPs,
                                            long minSegmentSize ) {

   // convenience renames
   double *fitMetric = analysisResults->fitMetric;

   // create nested for-loop iterators
   int *loopIterators = new int[numCPs+1];

   loopIterators[0] = 0;

   // initialize fit metric
   fitMetric[numCPs] = fitMetric[0];

   // start recursion
   ForIterate ( numCPs, loopIterators, numCPs, minSegmentSize );

   // normalize the fit-metric
   analysisResults->NormalizeFitMetric ( sequence->length );

   delete[] loopIterators;

};

/////////////////
// constructors
//

// defined inline in spec

////////////
// methods
//

// perform brute-force segmentation of sequence as per original
// fortran code
void BRUTE_FORCE_FAST :: Analyze ( long minSegmentSize ) {

   // ensure that segmentation is possible -- come to a grinding
   // halt if not
   if ( ( analysisResults->maxCPs + 1 ) * minSegmentSize > length ) {
      // can't perform this segmentation!
      printf ( "Error!  Impossible to perform segmentation with %d change\n",
               analysisResults->maxCPs );
      printf ( "points, min. seg. size of %d and sequence of length %d.\n",
               minSegmentSize, length );
      fflush ( stdin );
      getchar ( );
      exit ( 1 );

   };

   // convenience renames
   CHANGE_POINTS &changePoints = *(analysisResults->changePointMatrix);
   double *fitMetric = analysisResults->fitMetric;

   // fill in score matrix for current sequence
   ScoreQuasiDeviance ( sequence );

   // no need to initialize results matrix, since algorithm
   // doesn't care

   // no change-points.  the segment matrix is already set for this one.
   fitMetric[0] = (*scoreMatrix) ( 1, length );

   // find change points up to max number
   for ( int cpNum = 1; cpNum <= changePoints.numCPs; cpNum++ )
      FindChangePoints ( cpNum, minSegmentSize );

   // normalize the fit-metric
   analysisResults->NormalizeFitMetric ( sequence->length );

};

////////////////
// destructors
//

// defined inline in spec


////////////////////////////////////////////////////////
// subclass AUGER_LAWRENCE_FAST of ANALYSIS_ALGORITHM //
////////////////////////////////////////////////////////

// see Auger-Lawrence paper entitled "Algorithms for the Optimal
// Identification of Segment Neighborhoods" for the specification of
// this algorithm.

/////////////////
// constructors
//

// defined inline in spec

////////////
// methods
//

// perform analysis of segment using a segment matrix.  this optimizes
// the performance of the algorithm, but the method may not be used for
// sequences of large size.  function is split off from optimized function
// to avoid switching logic and thus to increase performance.
// TBD may wish to constrain function by sequence length check.  currently,
// assume caller knows what they are doing.
void AUGER_LAWRENCE_FAST :: Analyze ( long minSegmentSize ) {

   // ensure that segmentation is possible -- come to a grinding
   // halt if not
   if ( ( analysisResults->maxCPs + 1 ) * minSegmentSize > length ) {
      // can't perform this segmentation!
      printf ( "Error!  Impossible to perform segmentation with %d change\n",
               analysisResults->maxCPs );
      printf ( "points, min. seg. size of %d and sequence of length %d.\n",
               minSegmentSize, length );
      fflush ( stdin );
      getchar ( );
      exit ( 1 );

   };

   // convenience aliases
   CHANGE_POINTS &changePoints = *(analysisResults->changePointMatrix);
   double *fitMetric = analysisResults->fitMetric;

   // create support matrices of number of segments by length
   MATRIX<double> c ( analysisResults->maxCPs+1, length );
   MATRIX<long> cpos ( analysisResults->maxCPs+1, length );

   // fill in score matrix for current sequence
   ScoreQuasiDeviance ( sequence );

   // initialize results
   analysisResults->Initialize ( length );

   // see Auger-Lawrence paper for algorithm details.

   // step 1 already performed by scoring function:
   // populate for q = 1 from segment matrix
   for ( int segNum = 1; segNum <= length; segNum++ )
      c ( 1, segNum ) = (*scoreMatrix) ( 1, segNum );

   // step 2 iterate on remaining segment lengths:
   // populate for remaining q's.
   // q, c, j, v conform to the algorithm specification.
   for ( int q = 2; q <= analysisResults->maxCPs + 1; q++ ) {
      for ( int j = q * minSegmentSize;
            j <= length;
            j++ ) {

         // assume first element is the minimized metric
         c ( q, j ) = c ( q-1, (q-1)*minSegmentSize ) +
                      (*scoreMatrix) ( (q-1)*minSegmentSize+1, j );

         cpos ( q, j ) = (q - 1) * minSegmentSize;

         // check remaining elements to see if a more optimal solution
         // exists
         for ( int v = (q - 1) * minSegmentSize + 1;
               v <= j - minSegmentSize;
               v++ ) {
            if ( c ( q, j ) >= c ( q-1, v ) + (*scoreMatrix) ( v+1, j ) ) {
                c ( q, j ) = c ( q-1, v ) + (*scoreMatrix) ( v+1, j );
                cpos ( q, j ) = v;
            };
         };
      };
   };

   // populate change point data
   for ( int cpNum = analysisResults->maxCPs; cpNum >= 0; cpNum-- ) {
      changePoints ( cpNum, 0 ) = 0;
      changePoints ( cpNum+1, cpNum ) = length;
      for ( int cpVal = cpNum; cpVal >= 1; cpVal-- ) {
         changePoints ( cpVal, cpNum ) =
            cpos ( cpVal+1, changePoints ( cpVal+1, cpNum ) );
      };
   };

   // populate fitmetric (qd) data
   for ( int cpNum2 = 1; cpNum2 <= analysisResults->maxCPs + 1; cpNum2++ )
      fitMetric[cpNum2-1] = c ( cpNum2, length );

   // normalize the fit-metric
   analysisResults->NormalizeFitMetric ( sequence->length );

};

////////////////
// destructors
//

// defined inline in spec


///////////////////////////////////////////////////
// subclass AUGER_LAWRENCE of ANALYSIS_ALGORITHM //
///////////////////////////////////////////////////

// see Auger-Lawrence paper entitled "Algorithms for the Optimal
// Identification of Segment Neighborhoods" for the specification of
// this algorithm.

/////////////////
// constructors
//

// defined inline in spec

////////////
// methods
//

// perform analysis of sequence without using a segment matrix.  this
// is less efficient time-wise but permits analysis of arbitrary length
// sequences.  performance of the algorithm is a function of the number
// of change points with respect to the optimized version above.  function
// is split off from optimized function to avoid switching logic and thus
// to increase performance.
void AUGER_LAWRENCE :: Analyze ( long minSegmentSize ) {

   // ensure that segmentation is possible -- come to a grinding
   // halt if not
   if ( ( analysisResults->maxCPs + 1 ) * minSegmentSize >
        sequence->length ) {
      // can't perform this segmentation!
      printf ( "Error!  Impossible to perform segmentation with %d change\n",
               analysisResults->maxCPs );
      printf ( "points, min. seg. size of %d and sequence of length %d.\n",
               minSegmentSize, sequence->length );
      fflush ( stdin );
      getchar ( );
      exit ( 1 );

   };

   // convenience aliases
   CHANGE_POINTS &changePoints = *(analysisResults->changePointMatrix);
   double *fitMetric = analysisResults->fitMetric;

   // create support matrices of number of segments by length
   MATRIX<double> c ( analysisResults->maxCPs+1, sequence->length );
   MATRIX<long> cpos ( analysisResults->maxCPs+1, sequence->length );

   // current column of scoring matrix
   double *currentColumn = new double[sequence->length+1];

   // state for computing current column
   double mu[5];
   int muCount[5];
   int numElems;
   double sum;

// TBD
//char str[9];

   // initialize results
   analysisResults->Initialize ( sequence->length );

// TBD
//_strtime(str);
//printf("%s Segmenting using 1 change point...\n", str);

   // see Auger-Lawrence paper for algorithm details.

   // step 1 : populate for q = 1 from segment matrix; build this
   // dynamically for complete sequence

   // initialize column state
   numElems = 0;
   for ( int muNum = 1; muNum <= 4; muNum++ ) {
      mu[muNum] = 0.0;
      muCount[muNum] = 0;
   };

   // traverse complete sequence
   for ( int seqPos = 1; seqPos <= sequence->length; seqPos++ ) {
      // update state with new element added
      mu[sequence->data[seqPos]]++;
      muCount[sequence->data[seqPos]]++;
      numElems++;

      // compute quasi-deviance with new state
      sum = 0.0;
      for ( int muNum2 = 1; muNum2 <= 4; muNum2++ ) {
         if ( mu[muNum2] ) sum += muCount[muNum2] *
                                  log ( mu[muNum2] / numElems );
      };

      c ( 1, seqPos ) = -2 * sum;

   };

   // step 2 iterate on remaining segment lengths:
   // populate for remaining q's.
   // q, c, j, v conform to the algorithm specification.
   for ( int q = 2; q <= analysisResults->maxCPs + 1; q++ ) {

// TBD
//_strtime(str);
//printf("%s Segmenting using %d change points...\n", str, q);

      for ( int j = q * minSegmentSize;
            j <= sequence->length;
            j++ ) {
         // initialize column state
         numElems = 0;
         for ( int k2 = 1; k2 <= 4; k2++ ) {
            mu[k2] = 0.0;
            muCount[k2] = 0;
         };

         // optimally compute jth column of segment matrix
         for ( int row = j;
               row >= ( q - 1 ) * minSegmentSize;
               row-- ) {
            // update state with new element added
            mu[sequence->data[row]]++;
            muCount[sequence->data[row]]++;
            numElems++;

            // compute quasi-deviance with new state
            sum = 0.0;
            for ( int muPos3 = 1; muPos3 <= 4; muPos3++ ) {
               if ( mu[muPos3] ) sum += muCount[muPos3] *
                                        log ( mu[muPos3] / numElems );
            };

            currentColumn[row] = -2 * sum;

         };

         // assume first element is the minimized metric
         c ( q, j ) = c ( q-1, (q-1)*minSegmentSize ) +
                      currentColumn[(q-1)*minSegmentSize+1];

         cpos ( q, j ) = (q - 1) * minSegmentSize;

         // check remaining elements to see if a more optimal solution
         // exists
         for ( int v = (q - 1) * minSegmentSize + 1;
               v <= j - minSegmentSize;
               v++ ) {
            if ( c ( q, j ) >= c ( q-1, v ) + currentColumn[v+1] ) {
               // found new minimum
               c ( q, j ) = c ( q-1, v ) + currentColumn[v+1];
               cpos ( q, j ) = v;
            };
         };
      };
   };

   // populate change point data
   for ( int cpNum = analysisResults->maxCPs; cpNum >= 0; cpNum-- ) {
      changePoints ( cpNum, 0 ) = 0;
      changePoints ( cpNum+1, cpNum ) = sequence->length;
      for ( int cpVal = cpNum; cpVal >= 1; cpVal-- ) {
         changePoints ( cpVal, cpNum ) =
            cpos ( cpVal+1, changePoints ( cpVal+1, cpNum ) );
      };
   };

   // populate qd data
   for ( int cpNum2 = 1; cpNum2 <= analysisResults->maxCPs + 1; cpNum2++ )
      fitMetric[cpNum2-1] = c ( cpNum2, sequence->length );

   // normalize the fit-metric
   analysisResults->NormalizeFitMetric ( sequence->length );

   delete[] currentColumn;

};

////////////////
// destructors
//

// defined inline in spec


/////////////////////////////////////////////////////
// subclass MULTI_START_FAST of ANALYSIS_ALGORITHM //
/////////////////////////////////////////////////////

// see Braun, Braun, and Mueller paper entitled "Multiple Change-Point
// Fitting Via Quasi-Likelihood, With Applications To DNA Sequence
// Segmentation" for specification of this algorithm.

//////////////////////
// utility functions
//

// generates a random segmentation of a sequence.  the segmentation
// must obey minimal segment size constraints.
void MULTI_START_FAST :: GenerateNewGrid ( int *grid, int gridSize,
                                           RANDOM_NUMBER_GENERATOR &rndGen,
                                           long minSegmentSize ) {

   // number of generated grid elements thus far
   int numElems = 0;

   // prep grid
   grid[0] = 0;
   grid[1] = length;

   // we are finished if gridSize is 0
   if ( gridSize == 0 ) return;

   // define allowable interval for new CPs.  this interval will shrink
   // dynamically to optimally restrict point selection to regions of the
   // sequence not already covered by a previously selected point.
   int high = length - minSegmentSize;
   int highIndex = 1;
   int low = minSegmentSize;
   int lowIndex = 0;

   do {
      // generate a new random number in allowable interval
      int nextNum = (int) rndGen.GetRandomNumber ( ) * ( high - low + 1 ) + low;

      // search for proper location of new number
      for ( int gridPos = 1; gridPos <= numElems+1; gridPos++ ) {
         if ( nextNum < grid[gridPos] ) {
            // found the proper place for num, but is it valid?
            if ( ( nextNum >= grid[gridPos-1] + minSegmentSize ) &&
                 ( nextNum <= grid[gridPos] - minSegmentSize ) ) {
               // valid number indeed, insert it into the list
               numElems++;

               // shift successors up the list
               for ( int gridPos2 = numElems + 1;
                     gridPos2 > gridPos;
                     gridPos2-- ) {
                  grid[gridPos2] = grid[gridPos2-1];

               };

               // the highIndex shifted too
               highIndex++;

               grid[gridPos] = nextNum;

               // see if this number shrinks the allowable interval low
               for ( int lowPos = lowIndex+1;
                     lowPos <= numElems+1;
                     lowPos++ ) {
                  if ( grid[lowPos] < low + minSegmentSize ) {
                     // we can increase the lower bound
                     low = grid[lowPos] + minSegmentSize;
                     lowIndex = lowPos;

                  } else {
                     break;

                  };

               };

               // see if this number shrinks the allowable interval high
               for ( int hiPos = highIndex-1; hiPos >= 0; hiPos-- ) {
                  if ( grid[hiPos] > high - minSegmentSize ) {
                     // we can decrease the upper bound
                     high = grid[hiPos] - minSegmentSize;
                     highIndex = hiPos;

                  } else {
                     break;

                  };

               };

               // if adding new points is physically impossible, discard
               // this set of change points and begin anew
               if ( ( numElems < gridSize ) &&
                    ( high < low ) ) {
                  // reset number of generated grid elements thus far
                  numElems = 0;

                  // reprep grid
                  grid[0] = 0;
                  grid[1] = length;

                  // reset allowable interval for new CPs
                  high = length - minSegmentSize;
                  highIndex = 1;
                  low = minSegmentSize;
                  lowIndex = 0;

               };

            };    // if num is valid...

            // number has been processed, stop looking for proper location
            break;

         } else if ( nextNum == grid[gridPos] ) {
            // this number is useless, stop looking for proper location
            break;

         };

      };    // search for proper location of number

   // continue if we've found enough grid points
   } while ( numElems < gridSize );

};

// performs the segmentation analysis for a single number of change points.
// the scoring matrix needs to have been computed prior to this point!
void MULTI_START_FAST :: AnalyzeSingleCP ( long minSegmentSize,
                                           int numCPs ) {

   // ensure that segmentation is possible -- come to a grinding
   // halt if not
   if ( ( numCPs + 1 ) * minSegmentSize > length ) {
      // can't perform this segmentation!
      printf ( "Error!  Impossible to perform segmentation with %d change\n",
               numCPs );
      printf ( "points, min. seg. size of %d and sequence of length %d.\n",
               minSegmentSize, length );
      fflush ( stdin );
      getchar ( );
      exit ( 1 );

   };

   // convenience aliases
   CHANGE_POINTS &changePoints = *(analysisResults->changePointMatrix);
   double *fitMetric = analysisResults->fitMetric;

   // grid of change points
   int *grid = new int[numCPs+2];

   fitMetric[numCPs] = fitMetric[0];

   // run desired number of samples
   for ( int sample = 1; sample <= numSamples; sample++ ) {
      // generate a valid grid
      GenerateNewGrid ( grid, numCPs, randGen, minSegmentSize );

      // working vars
      double qd, qdTemp;

      // float movable change points
      for ( int cpNum = 1; cpNum <= numCPs; cpNum++ ) {
         // find new change point in interval
         qd = fitMetric[0];

         for ( int segPos = grid[cpNum-1] + minSegmentSize;
               segPos <= grid[cpNum+1] - minSegmentSize;
               segPos++ ) {
            qdTemp = (*scoreMatrix) ( grid[cpNum-1]+1, segPos ) +
                     (*scoreMatrix) ( segPos+1, grid[cpNum+1] );

            // check for new minimum
            if ( qdTemp <= qd ) {
               qd = qdTemp;
               grid[cpNum] = segPos;
            };

         };

      };

      // determine final qd
      qdTemp = 0;

      for ( int cpNum2 = 1; cpNum2 <= numCPs; cpNum2++ )
         qdTemp += (*scoreMatrix) ( grid[cpNum2-1]+1, grid[cpNum2] );

      qdTemp += (*scoreMatrix) ( grid[numCPs]+1, grid[numCPs+1] );

      // is this the new minimum configuration?
      if ( qdTemp <= fitMetric[numCPs] ) {
         fitMetric[numCPs] = qdTemp;

         // record change points
         for ( int cpNum3 = 1; cpNum3 <= numCPs; cpNum3++ )
            changePoints ( cpNum3, numCPs ) = grid[cpNum3];

      };

   };    // end of for sample

   delete[] grid;

};

/////////////////
// constructors
//

// defined inline in spec

////////////
// methods
//

// perform complete analysis of segment using multi-start method
void MULTI_START_FAST :: Analyze ( long minSegmentSize ) {

   // convenience aliases
   double *fitMetric = analysisResults->fitMetric;

   // fill in score matrix for current sequence
   ScoreQuasiDeviance ( sequence );

// TBD remove this when all is debugged
// local random number generator
//RANDOM_NUMBER_GENERATOR rndGen ( -2000 );

   // initialize results
   analysisResults->Initialize ( length );

   // no change-points.  the segment matrix is already set for this one.
   fitMetric[0] = (*scoreMatrix) ( 1, length );

   // find change points up to max number
   for ( int q = 1; q <= analysisResults->maxCPs; q++ ) {
      AnalyzeSingleCP ( minSegmentSize, q );

   };    // end for q

   // normalize the fit-metric
   analysisResults->NormalizeFitMetric ( sequence->length );

};

// perform an analysis of segment using multi-start method for only
// one change point number
void MULTI_START_FAST :: AnalyzeCP ( long minSegmentSize,
                                     int cpNum,
                                     boolean suppressScore ) {

   // convenience aliases
   double *fitMetric = analysisResults->fitMetric;

   // fill in score matrix for current sequence, if requested
   if ( ! suppressScore ) {
      ScoreQuasiDeviance ( sequence );

   };

   // initialize results
   analysisResults->Initialize ( length );

   // no change-points.  the segment matrix is already set for this one.
   fitMetric[0] = (*scoreMatrix) ( 1, length );

   // find change point for desired change point number
   AnalyzeSingleCP ( minSegmentSize, cpNum );

   // normalize the fit-metric
   analysisResults->NormalizeFitMetric ( sequence->length );

};

////////////////
// destructors
//

// defined inline in spec


////////////////////////////////////////////////
// subclass MULTI_START of ANALYSIS_ALGORITHM //
////////////////////////////////////////////////

// see Braun, Braun, and Mueller paper entitled "Multiple Change-Point
// Fitting Via Quasi-Likelihood, With Applications To DNA Sequence
// Segmentation" for specification of this algorithm.

//////////////////////
// utility functions
//

// generates a random segmentation of a sequence.  the segmentation
// must obey minimal segment size constraints.
void MULTI_START :: GenerateNewGrid ( int *grid, int gridSize,
                                      RANDOM_NUMBER_GENERATOR &rndGen,
                                      long minSegmentSize ) {

   // number of generated grid elements thus far
   int numElems = 0;

   // prep grid
   grid[0] = 0;
   grid[1] = sequence->length;

   // we are finished if gridSize is 0
   if ( gridSize == 0 ) return;

   // define allowable interval for new CPs.  this interval will shrink
   // dynamically to optimally restrict point selection to regions of the
   // sequence not already covered by a previously selected point.
   int high = sequence->length - minSegmentSize;
   int highIndex = 1;
   int low = minSegmentSize;
   int lowIndex = 0;

   do {
      // generate a new random number in allowable interval
      int nextNum = (int) rndGen.GetRandomNumber ( ) * ( high - low + 1 ) + low;

      // search for proper location of new number
      for ( int gridPos = 1; gridPos <= numElems+1; gridPos++ ) {
         if ( nextNum < grid[gridPos] ) {
            // found the proper place for num, but is it valid?
            if ( ( nextNum >= grid[gridPos-1] + minSegmentSize ) &&
                 ( nextNum <= grid[gridPos] - minSegmentSize ) ) {
               // valid number indeed, insert it into the list
               numElems++;

               // shift successors up the list
               for ( int gridPos2 = numElems + 1;
                     gridPos2 > gridPos;
                     gridPos2-- ) {
                  grid[gridPos2] = grid[gridPos2-1];

               };

               // the highIndex shifted too
               highIndex++;

               grid[gridPos] = nextNum;

               // see if this number shrinks the allowable interval low
               for ( int lowPos = lowIndex+1;
                     lowPos <= numElems+1;
                     lowPos++ ) {
                  if ( grid[lowPos] < low + minSegmentSize ) {
                     // we can increase the lower bound
                     low = grid[lowPos] + minSegmentSize;
                     lowIndex = lowPos;

                  } else {
                     break;

                  };

               };

               // see if this number shrinks the allowable interval high
               for ( int hiPos = highIndex-1; hiPos >= 0; hiPos-- ) {
                  if ( grid[hiPos] > high - minSegmentSize ) {
                     // we can decrease the upper bound
                     high = grid[hiPos] - minSegmentSize;
                     highIndex = hiPos;

                  } else {
                     break;

                  };

               };

               // if adding new points is physically impossible, discard
               // this set of change points and begin anew
               if ( ( numElems < gridSize ) &&
                    ( high < low ) ) {
                  // reset number of generated grid elements thus far
                  numElems = 0;

                  // reprep grid
                  grid[0] = 0;
                  grid[1] = sequence->length;

                  // reset allowable interval for new CPs
                  high = sequence->length - minSegmentSize;
                  highIndex = 1;
                  low = minSegmentSize;
                  lowIndex = 0;

               };

            };    // if num is valid...

            // number has been processed, stop looking for proper location
            break;

         } else if ( nextNum == grid[gridPos] ) {
            // this number is useless, stop looking for proper location
            break;

         };

      };    // search for proper location of number

   // continue if we've found enough grid points
   } while ( numElems < gridSize );

};

// performs the segmentation analysis for a single number of change points.
void MULTI_START :: AnalyzeSingleCP ( long minSegmentSize,
                                      int numCPs ) {

   // ensure that segmentation is possible -- come to a grinding
   // halt if not
   if ( ( numCPs + 1 ) * minSegmentSize > sequence->length ) {
      // can't perform this segmentation!
      printf ( "Error!  Impossible to perform segmentation with %d change\n",
               numCPs );
      printf ( "points, min. seg. size of %d and sequence of length %d.\n",
               minSegmentSize, sequence->length );
      fflush ( stdin );
      getchar ( );
      exit ( 1 );

   };

   // convenience aliases
   CHANGE_POINTS &changePoints = *(analysisResults->changePointMatrix);
   double *fitMetric = analysisResults->fitMetric;

   // grid of change points
   int *grid = new int[numCPs+2];

   fitMetric[numCPs] = fitMetric[0];

   // row and column of virtual score matrix
   double *currentRow    = new double[sequence->length+1];
   double *currentColumn = new double[sequence->length+1];

   // state for computing current row and column of virtual score matrix
   double mu[5];
   int muCount[5];
   int numElems;
   double sum;

   // run desired number of samples
   for ( int sample = 1; sample <= numSamples; sample++ ) {
      // generate a valid grid
      GenerateNewGrid ( grid, numCPs, randGen, minSegmentSize );

      // working vars
      double qd, qdTemp;

      // float movable change points
      for ( int cpNum = 1; cpNum <= numCPs; cpNum++ ) {
         // find new change point in interval
         qd = fitMetric[0];

         // compute virtual score matrix column
         // initialize column state
         numElems = 0;
         for ( int muPos = 1; muPos <= 4; muPos++ ) {
            mu[muPos] = 0.0;
            muCount[muPos] = 0;
         };

         // optimally compute appropriate column of virtual score matrix
         for ( int row = grid[cpNum+1];
               row >= grid[cpNum-1] + minSegmentSize;
               row-- ) {
            // update state with new element added
            mu[sequence->data[row]]++;
            muCount[sequence->data[row]]++;
            numElems++;

            // compute quasi-deviance with new state
            sum = 0.0;
            for ( int muPos1 = 1; muPos1 <= 4; muPos1++ ) {
               if ( mu[muPos1] ) sum += muCount[muPos1] *
                                        log ( mu[muPos1] / numElems );
            };

            currentColumn[row] = -2 * sum;

         };

         // compute virtual score matrix row
         // initialize row state
         numElems = 0;
         for ( int muPos2 = 1; muPos2 <= 4; muPos2++ ) {
            mu[muPos2] = 0.0;
            muCount[muPos2] = 0;
         };

         // optimally compute appropriate row of virtual score matrix
         for ( int column = grid[cpNum-1]+1;
               column <= grid[cpNum+1] - minSegmentSize;
               column++ ) {
            // update state with new element added
            mu[sequence->data[column]]++;
            muCount[sequence->data[column]]++;
            numElems++;

            // compute quasi-deviance with new state
            sum = 0.0;
            for ( int muPos3 = 1; muPos3 <= 4; muPos3++ ) {
               if ( mu[muPos3] ) sum += muCount[muPos3] *
                                        log ( mu[muPos3] / numElems );
            };

            currentRow[column] = -2 * sum;

         };

         // float current change point and minimize
         for ( int segPos = grid[cpNum-1] + minSegmentSize;
               segPos <= grid[cpNum+1] - minSegmentSize;
               segPos++ ) {
            qdTemp = currentRow[segPos] + currentColumn[segPos+1];

            // check for new minimum
            if ( qdTemp <= qd ) {
               qd = qdTemp;
               grid[cpNum] = segPos;
            };

         };

      };

      // determine final qd
      qdTemp = 0;

      for ( int cpNum2 = 1; cpNum2 <= numCPs; cpNum2++ )
         qdTemp += sequence->QuasiDeviance ( grid[cpNum2-1],
                                             grid[cpNum2] );

      qdTemp += sequence->QuasiDeviance ( grid[numCPs], grid[numCPs+1] );

      // is this the new minimum configuration?
      if ( qdTemp <= fitMetric[numCPs] ) {
         fitMetric[numCPs] = qdTemp;

         // record change points
         for ( int cpNum3 = 1; cpNum3 <= numCPs; cpNum3++ )
            changePoints ( cpNum3, numCPs ) = grid[cpNum3];

      };

   };    // end of for sample

   delete[] grid;
   delete[] currentRow;
   delete[] currentColumn;

};

/////////////////
// constructors
//

// defined inline in spec

////////////
// methods
//

// perform complete analysis of segment using multi-start method
void MULTI_START :: Analyze ( long minSegmentSize ) {

   // convenience aliases
   double *fitMetric = analysisResults->fitMetric;

// TBD remove this when all is debugged
// local random number generator
//RANDOM_NUMBER_GENERATOR rndGen ( -2000 );

   // initialize results
   analysisResults->Initialize ( sequence->length );

   // no change-points
   fitMetric[0] = sequence->QuasiDeviance ( 0, sequence->length );

   // find change points up to max number
   for ( int q = 1; q <= analysisResults->maxCPs; q++ ) {
      AnalyzeSingleCP ( minSegmentSize, q );

   };    // end for q

   // normalize the fit-metric
   analysisResults->NormalizeFitMetric ( sequence->length );

};

// perform an analysis of segment using multi-start method for only
// one change point number
void MULTI_START :: AnalyzeCP ( long minSegmentSize,
                                int cpNum ) {

   // convenience aliases
   double *fitMetric = analysisResults->fitMetric;

   // initialize results
   analysisResults->Initialize ( sequence->length );

   // no change-points.  the segment matrix is already set for this one.
   fitMetric[0] = sequence->QuasiDeviance ( 0, sequence->length );

   // find change point for desired change point number
   AnalyzeSingleCP ( minSegmentSize, cpNum );

   // normalize the fit-metric
   analysisResults->NormalizeFitMetric ( sequence->length );

};

////////////////
// destructors
//

// defined inline in spec



//////////////////////////////////////////////////
// subclass ORDERED_SPLIT of ANALYSIS_ALGORITHM //
//////////////////////////////////////////////////

// see jerry's e-mail for specification of this algorithm.

//////////////////////
// utility functions
//
void ORDERED_SPLIT::findBestChangePointAndQdDelta
   (TRACKING_NODE &node, double currentQd)
{
//JVB  oops!  the change is based only on the current segment
//JVB  rather than the whole sequence... no use to pass in currentQd...
     currentQd = sequence->QuasiDeviance( node.firstElem-1, node.lastElem );
     double minQd = currentQd;
//JVB   double minQd = currentQd;
//JVB
   // first check to see whether the sequence has been
   // split as far as it will go already.  If so, then
   // assign a negative value to the change -- analytically,
   // the drop in quasi-deviance must always be non-negative,
   // so that any other drop will look greater in comparison --
   // the point is to disallow this and this should probably
   // be handled in the calling routine
   if ( node.lastElem - node.firstElem <= 1  ) {
     node.qdDelta = -1.0;           // an arbitrary negative value
     node.cpElem = node.firstElem;  //
   }
   else {

   // by default, make the next segment only 1 (minimum)
   // observations wide so that even if nothing really is
   // changing when the observations are computed, at least
   // *something* plausible will be returned
     node.qdDelta = 0.0;
     node.cpElem = node.firstElem+1;
//JVB

   // TBD incorporate minSequenceSize in here -- you'll have
   // to pass it in as a parameter.  check these bounds for errors!
   // remember that QuasiDeviance works on (low, high].
   // check this logic, it is slightly tricky!!
   for (int seqIndex = node.firstElem+1;
//JVB        seqIndex < node.lastElem - 1;
        seqIndex <= node.lastElem - 1;
//JVB
        seqIndex++) {
      double newQd =
         // check logic!!
         sequence->QuasiDeviance(node.firstElem-1, seqIndex) +
         sequence->QuasiDeviance(seqIndex, node.lastElem);
      if (newQd < minQd) {
         minQd = newQd;
         node.cpElem = seqIndex;
      }
      // should probably test to ensure that node.cpElem gets set
      // at some point.  is it provable that a minQd less than
      // currentQd will always be found?
   }
   node.qdDelta = currentQd - minQd;
 } // else

//JVB Debug printf("Testing %d to %d, result %d with %lf\n",node.firstElem,node.lastElem,node.cpElem,node.qdDelta);

}

/////////////////
// constructors
//

// defined inline in spec

////////////
// methods
//

// perform complete analysis of segment using ordered-split method
void ORDERED_SPLIT :: Analyze ( long minSegmentSize ) {

   // ensure that segmentation is possible -- come to a grinding
   // halt if not.  this could be pulled out of here and added to
   // the ANALYSIS_ALGORITHM class as a (non-virtual) method called
   // checkSegmentation(long minSeg).  each algorithm that wants to
   // make this check just needs to call the method, e.g.
   // checkSegmentation(minSegmentSize) as its first action.
   // note that analysisResults and sequence are part of the
   // ANALYSIS_ALGORITHM and thus don't need to be passed in.
   if ((analysisResults->maxCPs + 1) * minSegmentSize >
       sequence->length) {
      // can't perform this segmentation!
      printf ("Error!  Impossible to perform segmentation with %d change\n",
              analysisResults->maxCPs);
      printf ("points, min. seg. size of %d and sequence of length %d.\n",
              minSegmentSize, sequence->length);
      fflush (stdin);
      getchar ();
      exit (1);

   };

   // convenience aliases
   CHANGE_POINTS &changePoints = *(analysisResults->changePointMatrix);
   double *fitMetric = analysisResults->fitMetric;

   // algorithm goes here.  it is natural to incorporate the
   // fitMetric and changePoints into the results object as
   // we go along.

   // bookkeeping array to track current options for new
   // change points
   TRACKING_NODE *trackingArray =
      new TRACKING_NODE[analysisResults->maxCPs+1];
   // next open node in array
   long nextNode = 0;

   // compute qd[0], the qd of the whole sequence
   fitMetric[0] = sequence->QuasiDeviance (0, sequence->length);

   // initialize limits of changePoints matrix
   for (int cp = 0; cp <= analysisResults->maxCPs; cp++) {
//JVB      changePoints ( cp, 0 ) = 0;
//JVB      changePoints ( cp+1, cp ) = sequence->length;
      changePoints ( 0, cp ) = 0;   // the cpth column
      changePoints ( cp+1, cp ) = sequence->length; // the cpth column
   };

   // create first tracking array node
   trackingArray[nextNode].firstElem = 1;
   trackingArray[nextNode].lastElem = sequence->length;

   // find best change point and qd delta for first node
   findBestChangePointAndQdDelta(trackingArray[nextNode],
                                 fitMetric[0]);
//JVB Debug printf("Quasi-deviance %lf\n",fitMetric[0]);
   nextNode++;

   for (int cpNum = 1; cpNum <= analysisResults->maxCPs; cpNum++) {
      // pick largest qd delta in tracking array
      int bestNode = 0;
      double bestQdDelta = trackingArray[bestNode].qdDelta;
//JVB Debug printf("cpNum %d\n",cpNum);
      for (int node = bestNode+1; node < nextNode; node++) {
//JVB Debug printf("Node %d, bestqdDelta %lf  contender %lf\n",node,bestQdDelta,trackingArray[node].qdDelta);
         if (trackingArray[node].qdDelta > bestQdDelta) {
            bestNode = node;
            bestQdDelta = trackingArray[node].qdDelta;
         }
      }

      // save current qd as last qd - qd delta
      fitMetric[cpNum] = fitMetric[cpNum-1] - bestQdDelta;

//JVB Debug printf("Quasi-deviance %lf\n", fitMetric[cpNum]);
      // copy rightmost tracking nodes up one position
      for (int shift = nextNode; shift > bestNode+1; shift--)
         trackingArray[shift] = trackingArray[shift-1];

      nextNode++;

      // split node into two and compute best cp and qd delta
      // for both of the new subsequences
      trackingArray[bestNode+1] = trackingArray[bestNode];
      trackingArray[bestNode].lastElem =
         trackingArray[bestNode].cpElem;
      findBestChangePointAndQdDelta(trackingArray[bestNode],
                                    fitMetric[cpNum]);
      trackingArray[bestNode+1].firstElem =
         trackingArray[bestNode+1].cpElem;
      findBestChangePointAndQdDelta(trackingArray[bestNode+1],
                                    fitMetric[cpNum]);

      // record new changePoints sequence from tracking array
      for (int cpIndex = 1; cpIndex <= cpNum; cpIndex++)
//JVB         changePoints ( cpNum, cpIndex ) =
         changePoints ( cpIndex, cpNum ) =      // that is, the cpNumth column
            trackingArray[cpIndex-1].lastElem;

   }

   delete[] trackingArray;

   // normalize the fit-metric
   analysisResults->NormalizeFitMetric (sequence->length);

};

////////////////
// destructors
//

// defined inline in spec

