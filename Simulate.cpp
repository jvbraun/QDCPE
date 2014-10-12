
//
// File: Simulate.cpp
//
// Purpose: This file defines the simulator used to investigate the
//    finite-sample properties of an analysis algorithm.  The simulator
//    generates multinomial test sequences and applies the analysis process
//    to the sequences.  The results are collected and summarized.
//
//	Author: Ron Braun
//
// History:
//		09/24/98	RKB	Created.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AnalAlg.h"
#include "DNASeqDt.h" 
#include "Globals.h"
#include "Random.h"
#include "ScoreFun.h"
#include "Simulate.h"


/////////////////////
// class SIMULATOR //
/////////////////////

//////////////////////
// utility functions

//**********************************************************************
//*  Calculate the IMSE of the estimator.
//**********************************************************************

void GET_IMSE ( double &IMSE, double &WIMSE, DNA_SEQUENCE_DATA &y,
                long RTRUE, double tau0[], MATRIX<double> &mu0,
                long RHAT, CHANGE_POINTS &CP ) {
   double mu[101][5], muhat[5];
   double V[5];

//*  First we need the estimated mean function by observation, mu_i.
   for ( int r1 = 1; r1 <= RHAT + 1; r1++ ) {
      y.ComputeMean ( muhat, CP ( r1-1, RHAT ), CP ( r1, RHAT ) );
      for ( int i1 = CP ( r1-1, RHAT ) + 1;
            i1 <= CP ( r1, RHAT );
            i1++ ) {
         for ( int j1 = 1; j1 <= 4; j1++ ) {
            mu[i1][j1] = muhat[j1];
         };

      };

   };

   IMSE = 0.0;
   WIMSE = 0.0;

   for ( int r = 1; r <= RTRUE + 1; r++ ) {
      for ( int j2 = 1; j2 <= 4; j2++ ) {
         V[j2] = mu0 ( r, j2 ) * ( 1 - mu0 ( r, j2 ) );
      };

      for ( int i = int ( y.length * tau0[r-1] ) + 1;
            i <= int ( y.length * tau0[r] );
            i++ ) {
         for ( int j = 1; j <= 4; j++ ) {
            IMSE = IMSE + ( mu[i][j] - mu0 ( r, j ) ) *
                          ( mu[i][j] - mu0 ( r, j ) );
            WIMSE = WIMSE + ( mu[i][j] - mu0 ( r, j ) ) *
                            ( mu[i][j] - mu0 ( r, j ) ) / V[j];

         };

      };

   };

   IMSE = IMSE / y.length;
   WIMSE = WIMSE / y.length;

};

//**********************************************************************
//*  Calculate the Schwarz criterion.  The
//*  penalty is n^1/4.  Note that the quasi-deviance is already
//*  assumed normalized in this picture.
//**********************************************************************

long schwarz ( double *qd, long N, long RMAX, double penalty ) {
   double sc[32];
   long RHAT;

   for ( int r1 = 0; r1 <= RMAX; r1++ ) {
      sc[r1] = 0.5 * N * log ( qd[r1] );
      sc[r1] = sc[r1] + r1 * penalty;
   };

   RHAT = 0;

   for ( int r = 1; r <= RMAX; r++ ) {
      if ( sc[r] < sc[RHAT] )
         RHAT = r;

   };

   return RHAT;

};

//*********************************************************************
//*  Calculate the uniformity statistic.
//*********************************************************************

void GET_U ( double &U, long RTRUE, double TAU0[32] ) {

   U = 0.0;
   if ( RTRUE > 0 ) {
      for ( int r = 1; r <= RTRUE + 1; r++ ) {
         U = U + pow ( ( TAU0[r] - TAU0[r-1] ) - 1 / ( RTRUE + 1 ), 2 );
      };
   };

};

//*********************************************************************
//*  Calculate the SSJ.
//*********************************************************************

double GET_SSJ ( long RTRUE, MATRIX<double> *MU0 ) {
   // computed value of SSJ
   double result = 0.0;

   for ( int r = 1; r <= RTRUE; r++ ) {
      for ( int j = 1; j <= 4; j++ ) {
         result = result + pow ( (*MU0) ( r+1, j ) - (*MU0) ( r, j ), 2 );

      };
   };

   // return computed SSJ
   return result;

};

//*********************************************************************
//*  Get the estimated mean parameters.
//*********************************************************************

MATRIX<double> *GET_MUH ( DNA_SEQUENCE_DATA &y, long RHAT,
               CHANGE_POINTS &CP ) {
   double MUtemp[5];

   MATRIX<double> *MUH = new MATRIX<double> ( RHAT+1, y.alphaLength );

//*  We need the estimated mean function.
   for ( int r = 1; r <= RHAT + 1; r++ ) {
      y.ComputeMean ( MUtemp, CP ( r-1, RHAT ), CP ( r, RHAT ) );
      for ( int j = 1; j <= 4; j++ ) {
         MUH->Data ( r, j ) = MUtemp[j];
      };
   };

   return MUH;

};


/////////////////////////////////////
// subclass BENCHMARK of SIMULATOR //
/////////////////////////////////////

// this simulator corresponds to the original simulation code developed
// in fortran by J. Braun.  it is retained as a benchmark for use in code
// development, since it generates a known set of outputs.

////////////
// methods

// perform simulation scenario for benchmark
void BENCHMARK :: Simulate ( void ) {

//*  General set-up.  We do not expect to use sample sizes larger
//*  than N=100; also it seems safe that we select less than 30
//*  change-points.
   long N, NS, RMAX, RATE;
   long i, r, s, RHAT, sim, ind;
   long sim_N[10];
   double U, SSJ, IMSE, WIMSE;
   long SWAP[32];
   MATRIX<double> *MUH;
   double sim_qd[10][10], sim_p2[10][10];
   char sample[8] = " abcdef";
   char filename[20];

   // parameter setup

   // specify which algorithm to use; this stuff will be tied to the
   // command line eventually, so leave it crude now...
//   ALGORITHMS algorithm = BRUTE_FORCE_FAST_ALGORITHM;
//   ALGORITHMS algorithm = AUGER_LAWRENCE_FAST_ALGORITHM;
//   ALGORITHMS algorithm = AUGER_LAWRENCE_ALGORITHM;
//   ALGORITHMS algorithm = MULTI_START_FAST_ALGORITHM;
   //   ALGORITHMS algorithm = MULTI_START_ALGORITHM;
   ALGORITHMS algorithm = ORDERED_SPLIT_ALGORITHM;
   ANALYSIS_ALGORITHM *analAlg;

//*  Initialize the output files.

   FILE *fd[10];

   fd[1] = fopen ( "final4_1.dat", "wt" );
   fd[2] = fopen ( "final4_2.dat", "wt" );
   fd[3] = fopen ( "final4_3.dat", "wt" );
   fd[4] = fopen ( "final4_4.dat", "wt" );

   for ( i = 1; i <= 4; i++ ) {
      if ( fd[i] == NULL ) {
         printf ( "Cannot open a data file.\n" );
         fflush ( stdin );
         getchar ( );
         exit ( 1 );

      };
   };

//*------------------ SIMULATION CONTROL ------------------------------*
//*  N is the data size, NS is the number of simulations we will perform.
   NS = 1000;
   RMAX = 5;

//*  We roll through the different sample sizes and rates.
// TBD
   for ( N = 50; N <= 100; N += 10 ) {
      // create empty sequence of length N with at most
      // RMAX change points
      MULTINOMIAL_SEQUENCE y ( N, RMAX );

      // create analysis algorithm to be used on sequence
      // TBD make this a utility
      switch ( algorithm ) {
         case BRUTE_FORCE_FAST_ALGORITHM :
            analAlg = new BRUTE_FORCE_FAST ( &y, RMAX );
         break;

         case AUGER_LAWRENCE_FAST_ALGORITHM :
            analAlg = new AUGER_LAWRENCE_FAST ( &y, RMAX );
         break;

         case AUGER_LAWRENCE_ALGORITHM :
            analAlg = new AUGER_LAWRENCE ( &y, RMAX );
         break;

         case MULTI_START_FAST_ALGORITHM :
            analAlg = new MULTI_START_FAST ( &y, RMAX, 100 );
         break;

         case MULTI_START_ALGORITHM :
            analAlg = new MULTI_START ( &y, RMAX, 100 );
         break;

         case ORDERED_SPLIT_ALGORITHM :
            analAlg = new ORDERED_SPLIT ( &y, RMAX );
         break;

         default :
            printf ( "No algorithm specified!\n" );
            fflush ( stdin );
            getchar ( );
            exit ( 1 );

      };

//*  Initialize the simulation statistics matrices.
      for ( r = 0; r <= RMAX; r++ ) {
         sim_N[r] = 0;
         for ( s = 0; s <= RMAX; s++ ) {
            sim_qd[r][s] = 0;
            sim_p2[r][s] = 0;
         };
      };

      for ( RATE = 1; RATE <= 5; RATE += 2 ) {

         // let us SIMULATE...
         for ( sim = 1; sim <= NS; sim++ ) {

            // simulate some data from the multinomial distribution
            // with our mean function
            y.GenerateNewSequence ( RATE, .2, (long) pow ( y.length, 0.5 ) );

            // perform analysis
            analAlg->Analyze ( (long) pow ( y.length, 0.5 ) );

//*  Now we are finished estimating change-points.  Let us go
//*  on to other exciting tasks!

            // define some useful aliases
            double *fitMetric = analAlg->analysisResults->fitMetric;
            CHANGE_POINTS &changePoints =
               *( analAlg->analysisResults->changePointMatrix );

//*  Let us set down some statistics which describe the true
//*  situation.  They are the true SSJ and the statistic U.
            SSJ = GET_SSJ ( y.rTrue, y.mu0 );
            GET_U ( U, y.rTrue, y.tau0 );

            fprintf ( fd[1], " %1d %3.2f %6.2f\n", y.rTrue, U, SSJ );

//*  Now suppose we knew exactly where the change-points were.
//*  How would the estimation of IMSE, WIMSE, and SSJ proceed?
//*  In some sense this provides a sensible baseline to measure
//*  the effects of estimating the number of change-points.

//*  We will have to swap the appropriate column of CP in and out.
            for ( int r1 = 1; r1<= y.rTrue; r1++ ) {
               SWAP[r1] = changePoints ( r1, y.rTrue );
               changePoints ( r1, y.rTrue ) =
                  long ( y.length * y.tau0[r1] );
            };

            GET_IMSE ( IMSE, WIMSE, y, y.rTrue, y.tau0, *(y.mu0),
                       y.rTrue, changePoints );
            MUH = GET_MUH ( y, y.rTrue, changePoints );
            SSJ = GET_SSJ ( y.rTrue, MUH );
            delete MUH;

            for ( int r2 = 1; r2 <= y.rTrue; r2++ ) {
               changePoints ( r2, y.rTrue ) = SWAP[r2];
            };

            fprintf ( fd[2], " %6.2f %6.2f %6.2f\n", IMSE, WIMSE, SSJ );

//*  We would also like to see how the procedure does with the
//*  known number of change-points.
            GET_IMSE ( IMSE, WIMSE, y, y.rTrue, y.tau0, *(y.mu0),
                       y.rTrue, changePoints );
            MUH = GET_MUH ( y, y.rTrue, changePoints );
            SSJ = GET_SSJ ( y.rTrue, MUH );
            delete MUH;

            fprintf ( fd[3], " %6.2f %6.2f %6.2f\n", IMSE, WIMSE, SSJ );

//*  Let us see how the procedure does with the estimated number
//*  of change-points.
            RHAT = schwarz ( fitMetric, y.length, RMAX, pow ( N, 0.25 ) );

            GET_IMSE ( IMSE, WIMSE, y, y.rTrue, y.tau0, *(y.mu0),
                       RHAT, changePoints );
            MUH = GET_MUH ( y, RHAT, changePoints );
            SSJ = GET_SSJ ( RHAT, MUH );
            delete MUH;

            fprintf ( fd[4], " %1d %6.4f %6.4f %6.2f\n", RHAT,IMSE, WIMSE, SSJ );

//*  Finally, we will jot down the deviance results to this point.
            sim_N[y.rTrue] = sim_N[y.rTrue] + 1;
            for ( int r3 = 0; r3 <= RMAX; r3++ ) {
               sim_qd [y.rTrue][r3] = sim_qd[y.rTrue][r3] + fitMetric[r3];
            };


//*------------------ SIMULATION CONTROL ------------------------------*
         };    // for sim

      };    // for RATE

//*  Now write down the simulation statistics.
      i = y.length / 10 - 4;

      strcpy ( filename, "final4_" );
      strncat ( filename, &sample[i], 1 );
      strcat ( filename, ".qd" );

      fd[9] = fopen ( filename, "wt" );

      if ( fd[9] == NULL )
         printf ( "Cannot open file %s\n", filename );

      for ( int r4 = 0; r4 <= RMAX; r4++ ) {
         if ( sim_N[r4] > 0) {
            for ( int s1 = 0; s1 <= RMAX; s1++ ) {
               sim_qd[r4][s1] = sim_qd[r4][s1] / sim_N[r4];
               sim_p2[r4][s1] = sim_p2[r4][s1] / sim_N[r4];
            };
         };

         fprintf ( fd[9], " %1d ", r4 );
         for ( ind = 0; ind <= RMAX; ind++ )
            fprintf ( fd[9], "%8.4f ", sim_qd[r4][ind] );
         fprintf ( fd[9], "%4d\n", sim_N[r4] );

      };

      fclose ( fd[9] );

      // clean up dynamics
      delete analAlg;

   };    // for N

//*--------------------------------------------------------------------*

   for ( int i1 = 1; i1 <= 4; i1++ )
      fclose ( fd[i1] );

};


///////////////////////////////////
// subclass PHASE_I of SIMULATOR //
///////////////////////////////////

// this simulator corresponds to the first series of simulations described
// in section 4 or "Multiple Change-Point Fitting Via Quasi-Likelihood,
// With Application To DNA Sequence Seqmentation".

/////////////////
// constructors
//

////////////
// methods
//

// perform phase I simulation
void PHASE_I :: Simulate ( void ) {

   /////////////////////////////
   // simulation parameter setup

   // varying sequence size
   int const SEQ_SIZE_LOW  = 100;
   int const SEQ_SIZE_HIGH = 1000;
   int const SEQ_SIZE_STEP = 100;

   // varying number of change points
   int const NUM_CPS_LOW  = 1;
   int const NUM_CPS_HIGH = 10;
   int const NUM_CPS_STEP = 1;

   // varying jump sizes
   double const JUMP_SIZE_LOW  = 0.2;
   double const JUMP_SIZE_HIGH = 0.8;
   double const JUMP_SIZE_STEP = 0.3;

   // number of replications on a given cell
   int const NUM_CELL_REPS = 10;

   // max number of tries for multi-start to generate proper change
   // points before aborting attempt
//JVB   int const MULTI_START_THRESHHOLD = 10000;

   // minimum sub-sequence length in segmentation
   int const MIN_SEGMENT_SIZE = 1;

   ///////////////////////
   // begin the simulation

   // open file for simulation results
   FILE *fd = fopen ( "PhaseI.dat", "wt" );
//JVB Debug FILE *fd = stdout;

   // make sure it opened
   if ( fd == NULL ) {
      printf ( "Cannot open data file for phase I results.\n" );
      fflush ( stdin );
      getchar ( );
      exit ( 1 );

   };

   // vary sequence sizes
   for ( int seqSize = SEQ_SIZE_LOW;
         seqSize <= SEQ_SIZE_HIGH;
         seqSize += SEQ_SIZE_STEP ) {

      // create empty sequence of length seqSize with at most
      // the maximum number of change points
      MULTINOMIAL_SEQUENCE sequence ( seqSize, NUM_CPS_HIGH );

      // TBD it seems superfluous to pass in the number of CPs to
      // both the sequence and the algorithm.  investigate whether
      // this is necessary and remove the offending parameter if not.

      // vary number of change points
      for ( int numCPs = NUM_CPS_LOW;
            numCPs <= NUM_CPS_HIGH;
            numCPs += NUM_CPS_STEP ) {

         // create auger-lawrence algorithms for this sequence
         AUGER_LAWRENCE_FAST augLawAlg  ( &sequence, numCPs );

//JVB         // create multi-start algorithm for this sequence.  we will
//JVB         // run it for one iteration each time.
//JVB         MULTI_START_FAST multiStartAlg ( &sequence, numCPs, 1 );
         // create ordered split algorithm for this sequence.
         ORDERED_SPLIT orderedSplitAlg ( &sequence, numCPs );
//JVB

         // define some convenience aliases
         CHANGE_POINTS &augLawCPs =
                          *(augLawAlg.analysisResults->changePointMatrix);
//JVB         CHANGE_POINTS &multiStartCPs =
//JVB                          *(multiStartAlg.analysisResults->changePointMatrix);
         CHANGE_POINTS &orderedSplitCPs =
                          *(orderedSplitAlg.analysisResults->changePointMatrix);
//JVB

         // vary jump sizes
         for ( double jumpSize = JUMP_SIZE_LOW;
               jumpSize <= JUMP_SIZE_HIGH;
               jumpSize += JUMP_SIZE_STEP ) {

            // run prescribed number of replications for multi-start
            for ( int numReps = 1;
                  numReps <= NUM_CELL_REPS;
                  numReps++ ) {

               // generate test sequence around the current model
               sequence.GenerateNewSequence ( numCPs, jumpSize,
                                              MIN_SEGMENT_SIZE );

               // perform auger-lawrence gold-standard analysis
               augLawAlg.Analyze ( MIN_SEGMENT_SIZE );

               // write results and actual segmentation for this test
               fprintf ( fd, " %6d %4d %6.2f\n", seqSize, numCPs, jumpSize );
               augLawCPs.Print ( fd, numCPs );

//JVB               // since we already did the work to generate the score
//JVB               // matrix for this sequence, let the multi-start algorithm
//JVB               // use the same score matrix
//JVB               multiStartAlg.UseAlternateScoreMatrix
//JVB                  ( augLawAlg.GetScoreMatrix ( ) );
//JVB
//JVB               // perform multi-start until we find a match or reach the
//JVB               // multi-start threshhold for each replication
//JVB               int numTrys = 0;
//JVB
//JVB               for ( int multiStartNum = 1;
//JVB                     multiStartNum <= MULTI_START_THRESHHOLD;
//JVB                     multiStartNum++ ) {
//JVB                  // run the multi-start algorithm once, suppressing the
//JVB                  // rescoring of the sequence since we are using the
//JVB                  // auger lawrence score matrix already computed
//JVB                  multiStartAlg.AnalyzeCP ( MIN_SEGMENT_SIZE, numCPs, TRUE );
                  orderedSplitAlg.Analyze ( MIN_SEGMENT_SIZE );

//JVB

//JVB                  // compare change points with that of auger-lawrence
//JVB                  // and stop if we find a match
//JVB                  boolean cpsMatch = TRUE;
//JVB
//JVB                  for ( int cpPos = 1; cpPos <= numCPs; cpPos++ ) {
//JVB                     // see if this change point differs
//JVB                     if ( augLawCPs ( cpPos, numCPs ) !=
//JVB                          multiStartCPs ( cpPos, numCPs ) ) {
//JVB                          orderedSplitCPs ( cpPos, numCPs ) ) {
//JVB                        // match failed
//JVB                        cpsMatch = FALSE;
//JVB
//JVB                        // cease trying
//JVB                        break;
//JVB
//JVB                     };
//JVB
//JVB                  };    // end of for cpsMatch
//JVB
//JVB                  // quit if we have a match
//JVB                  if ( cpsMatch ) {
//JVB                     // annotate number of reps needed
//JVB                     numTrys = multiStartNum;
//JVB
//JVB                     // dump out of multi-start loop
//JVB                     break;
//JVB
//JVB                  };
//JVB
//JVB               };    // end of for multiStartNum

               // write results and actual segmentation for this test
//JVB               fprintf ( fd, " %6d %4d %6.2f %6d ", seqSize, numCPs,
//JVB                         jumpSize, numTrys );

//JVB            augLawCPs.Print ( fd, numCPs );
//JVB
               orderedSplitCPs.Print (fd, numCPs );
//JVB

               // flush the file for real-time tracking
               fflush ( fd );

            };    // end of for numReps

         };    // end of for jumpsize

      };    // end of for numCPs

   };    // end of for seqSize

   // close simulation results file
   fclose ( fd );

};

////////////////
// destructors
//


///////////////////////////////////
// subclass PHASE_II of SIMULATOR //
///////////////////////////////////

// this simulator corresponds to the second series of simulations described
// in section 4 or "Multiple Change-Point Fitting Via Quasi-Likelihood,
// With Application To DNA Sequence Seqmentation".

/////////////////
// constructors
//

////////////
// methods
//

// perform phase I simulation
void PHASE_II :: Simulate ( void ) {

   /////////////////////////////
   // simulation parameter setup

   // varying number of change points
   int const NUM_CPS_LOW  = 0;
//JVB   int const NUM_CPS_HIGH = 50;
   int const NUM_CPS_HIGH = 20;
   int const NUM_CPS_STEP = 1;

   // varying jump sizes
   double const JUMP_SIZE_LOW  = 0.2;
   double const JUMP_SIZE_HIGH = 0.8;
   double const JUMP_SIZE_STEP = 0.3;

   // number of replications on a given cell
   int const NUM_CELL_REPS = 2;

   // number of times multi-start should be executed in minimizing
   // the fit metric
   int const MULTI_START_REPS = 100;

   // minimum sub-sequence length in segmentation
   int const MIN_SEGMENT_SIZE = 1;

   ///////////////////////
   // begin the simulation

   // open file for simulation results
   FILE *fd = fopen ( "PhaseII.dat", "wt" );

   // make sure it opened
   if ( fd == NULL ) {
      printf ( "Cannot open data file for phase II results.\n" );
      fflush ( stdin );
      getchar ( );
      exit ( 1 );

   };

   // define sequence lengths o~f interest
//JVB   int const NUM_SEQ_LENS = 6;
//JVB
//JVB   int seqSizes[NUM_SEQ_LENS+1] =
//JVB          { 0, 1000, 5000, 10000, 50000, 100000, 500000 };
   int const NUM_SEQ_LENS = 3;
   int seqSizes[NUM_SEQ_LENS+1] =
          { 0, 1000, 10000, 100000 };

//JVB ADD THE FOLLOWING LINES TO //JVB
   int const NUM_CP_LEVELS = 5;
   int CPLevels[NUM_CP_LEVELS+1] =
          { 0, 1, 5, 10, 20 };
//JVB

   // vary sequence sizes
   for ( int seqPos = 1; seqPos <= NUM_SEQ_LENS; seqPos++ ) {

      // get sequence length
      int seqSize = seqSizes[seqPos];

      // create empty sequence of length seqSize with at most
      // the maximum number of change points
      MULTINOMIAL_SEQUENCE sequence ( seqSize, NUM_CPS_HIGH );

      // TBD it seems superfluous to pass in the number of CPs to
      // both the sequence and the algorithm.  investigate whether
      // this is necessary and remove the offending parameter if not.

      // vary true number of change points
//JVB      for ( int numCPs = NUM_CPS_LOW;
//JVB            numCPs <= NUM_CPS_HIGH;
//JVB            numCPs += NUM_CPS_STEP ) {
//JVB SUBSTITUTE THE FOLLOWING 4 LINES
      for ( int CPPos = 1; CPPos <= NUM_CP_LEVELS; CPPos++) {

         // get true number of change-points
         int numCPs = CPLevels[CPPos];

         // create multi-start algorithm for this sequence
//JVB         MULTI_START multiStartAlg ( &sequence, numCPs,
//JVB                                     MULTI_START_REPS );
//JVB         AUGER_LAWRENCE_FAST AugLawAlg ( &sequence, NUM_CPS_HIGH );
//JVB         AUGER_LAWRENCE AugLawAlg ( &sequence, NUM_CPS_HIGH );
         ORDERED_SPLIT OrdSplitAlg ( &sequence, NUM_CPS_HIGH );
//JVB

         // define some convenience aliases
//JVB         double *fitMetric = multiStartAlg.analysisResults->fitMetric;
//JVB         double *fitMetric = AugLawAlg.analysisResults->fitMetric;
         double *fitMetric = OrdSplitAlg.analysisResults->fitMetric;
//JVB

         // vary jump sizes
         for ( double jumpSize = JUMP_SIZE_LOW;
               jumpSize <= JUMP_SIZE_HIGH;
               jumpSize += JUMP_SIZE_STEP ) {

            // run prescribed number of replications for multi-start
            for ( int numReps = 1;
                  numReps <= NUM_CELL_REPS;
                  numReps++ ) {

               // generate test sequence around the current model
               sequence.GenerateNewSequence ( numCPs, jumpSize,
                                              MIN_SEGMENT_SIZE );

               // run the multi-start algorithm
//JVB               multiStartAlg.AnalyzeCP ( MIN_SEGMENT_SIZE, numCPs );
//JVB               AugLawAlg.Analyze ( MIN_SEGMENT_SIZE );
               OrdSplitAlg.Analyze ( MIN_SEGMENT_SIZE );
//JVB

               // write results for this test
//JVB               fprintf ( fd, " %6d %4d %6.2f %10.6f \n", seqSize, numCPs,
//JVB                         jumpSize, fitMetric[numCPs] );
               fprintf ( fd, " %6d %4d %6.2f ", seqSize, numCPs, jumpSize);
               for (int numMetricResults = NUM_CPS_LOW;
                    numMetricResults <= NUM_CPS_HIGH;
                    numMetricResults++ ) {
                 fprintf ( fd, "%10.6f ", fitMetric[numMetricResults] );
               }
               fprintf (fd, "\n");
//JVB

               // flush the file for real-time tracking
               fflush ( fd );

            };    // end of for numReps

         };    // end of for jumpsize

      };    // end of for numCPs

   };    // end of for seqSize

   // close simulation results file
   fclose ( fd );

};

////////////////
// destructors
//


