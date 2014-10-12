
//
// File: AnalRes.cpp
//
// Purpose: This file defines the results of a sequence analysis, output
// by the sequence analysis process.
//
//	Author: Ron Braun
//
// History:
//		09/29/98	RKB	Created.
//

#include "AnalRes.h"


/////////////////////////
// class CHANGE_POINTS //
/////////////////////////

/////////////////
// constructors

CHANGE_POINTS :: CHANGE_POINTS ( long num ) {
   // annotate upper bound on number of change points
   numCPs = num;

   // allocate new CP matrix of necessary size
   CPMatrix = new MATRIX<long> ( num+1, num+1 );

};

////////////
// methods

// initialize the matrix of change-point estimators.
// the matrix contains the estimates (in the N scale).
void CHANGE_POINTS :: Initialize ( long N ) {
   for ( int r = 0; r <= numCPs; r++ ) {
      for ( int s = 0; s <= r+1; s++ ) {
         CPMatrix->Data ( s, r ) = 0;
      };

      CPMatrix->Data ( r+1, r ) = N;

   };

};

// print segmentation for num change points
void CHANGE_POINTS :: Print ( long num ) {
   printf ( "CP(%d) : ( ", num );

   for ( int i = 0; i <= num+1; i++ )
      printf ( "%d ", CPMatrix->Data ( i, num ) );

   printf ( ")\n" );

};

// print segmentation for num change points to a file
void CHANGE_POINTS :: Print ( FILE *fd, long num ) {
   fprintf ( fd, "   ( ", num );

   for ( int i = 0; i <= num+1; i++ )
      fprintf ( fd, "%d ", CPMatrix->Data ( i, num ) );

   fprintf ( fd, ")\n" );

};

////////////////
// destructors

CHANGE_POINTS :: ~CHANGE_POINTS ( void ) {
//JVB   free ( CPMatrix );
delete CPMatrix;

};


////////////////////////////
// class ANALYSIS_RESULTS //
////////////////////////////

////////////
// methods

// initialize results arrays
void ANALYSIS_RESULTS :: Initialize ( long lastElement ) {
   // initialize change points using the last element number
   changePointMatrix->Initialize ( lastElement );

   // init fitMetric
   for ( int i = 0; i <= maxCPs; i++ )
      fitMetric[i] = 0;

};

// normalize the fit metric
void ANALYSIS_RESULTS :: NormalizeFitMetric ( long norm ) {
   for ( int r = 0; r <= maxCPs; r++ )
      fitMetric[r] = fitMetric[r] / norm;

};

// print results in nice format
void ANALYSIS_RESULTS :: Print ( void ) {
   printf ( "\nSegmented sequence using up to %d change points.\n\n",
            maxCPs );

   for ( int i = 1; i <= maxCPs; i++ ) {
      changePointMatrix->Print ( i );
      printf ( "Normalized fit metric: %8.6f\n\n", fitMetric[i] );
   };

   printf ( "\n" );

};

