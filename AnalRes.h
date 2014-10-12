
//
// File: AnalRes.h
//
// Purpose: This file defines The results of a sequence analysis, output
// by the sequence analysis process.
//
//	Author: Ron Braun
//
// History:
//		07/30/98	RKB	Created.
//

#ifndef ANALRESH
#define ANALRESH

#include <stdio.h>

#include "Matrix.h"


/////////////////////////
// class CHANGE_POINTS //
/////////////////////////

class CHANGE_POINTS {
   // list of change points
   MATRIX<long> *CPMatrix;

public:
   // upper bound on number of desired change points
   long numCPs;

   /////////////////
   // constructors
   CHANGE_POINTS ( long num );

   ////////////
   // methods

   // initialize the matrix of change-point estimators.
   // the matrix contains the estimates (in the N scale).
   void Initialize ( long N );

   // access a change point location using convenient notation
   long &operator () ( int x, int y ) {
      return CPMatrix->Data ( x, y );
   };

   // access a change point location using brute notation
   long &Data ( int x, int y ) {
      return CPMatrix->Data ( x, y );
   };

   // print segmentation for num change points
   void Print ( long num );
   void Print ( FILE *fd, long num );

   ////////////////
   // destructors
   ~CHANGE_POINTS ( void );

};


////////////////////////////
// class ANALYSIS_RESULTS //
////////////////////////////

class ANALYSIS_RESULTS {
public:
   // change point segmentation derived from input score matrix
   CHANGE_POINTS *changePointMatrix;

   // fit metric specifier related to change point.  this may be a
   // measure like quasi-deviance or some other fitting function.
   double *fitMetric;

   // upper bound on number of change points to be considered
   long maxCPs;

   /////////////////
   // constructors

   ANALYSIS_RESULTS ( long numCPs ) {
      // annotate max number of change points
      maxCPs = numCPs;

      // create CP matrix for results of this algorithm
      changePointMatrix = new CHANGE_POINTS ( numCPs );

      // create fit metric array
      fitMetric = new double[numCPs+1];

   };

   ////////////
   // methods

   // initialize results arrays
   void Initialize ( long lastElement );

   // normalize the fit metric
   void NormalizeFitMetric ( long norm );

   // print results in nice format
   void Print ( void );

   ////////////////
   // destructors
   ~ANALYSIS_RESULTS ( void ) {
      delete changePointMatrix;
      delete[] fitMetric;

   };

};

#endif

