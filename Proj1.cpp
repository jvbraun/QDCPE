
//
// File: Proj1.cpp
//
// Purpose: This file defines the control code for the PROJ1 application.
//    This provides the glue functionality to tie together control of the
//    other software components.  This item also includes the interface to
//    the application in toto.
//
//	Author: Ron Braun
//
// History:
//		08/06/98	RKB	Created.
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
#include "SeqAnal.h"
#include "Simulate.h"


void main ( ) {
   // the following specifications may be tied to the command line
   // or some other interface later.  for the moment, we'll specify
   // control crudely here.

   // choose task for proj1
//   TASKS taskType = SIMULATION_TASK;
   TASKS taskType = DATA_ANALYSIS_TASK;

   // for SIMULATION, choose simulator type
//   SIMULATIONS simulationType = BENCHMARK_SIMULATION;
   SIMULATIONS simulationType = PHASE_I_SIMULATION;
//   SIMULATIONS simulationType = PHASE_II_SIMULATION;

   // for DATA_ANALYSIS, choose which algorithm to use
//   ALGORITHMS algorithmType = BRUTE_FORCE_FAST_ALGORITHM;
// ALGORITHMS algorithmType = AUGER_LAWRENCE_FAST_ALGORITHM;
    ALGORITHMS algorithmType = AUGER_LAWRENCE_ALGORITHM;
//   ALGORITHMS algorithmType = MULTI_START_FAST_ALGORITHM;
//   ALGORITHMS algorithmType = MULTI_START_ALGORITHM;
//   ALGORITHMS algorithmType = ORDERED_SPLIT_ALGORITHM;

   // for DATA_ANALYSIS, set the following parameters
   char seqFileName[15] = "lamcg.fas";
   int maxCPs = 50;
   long minLength = 1;
   double cnPower = 0.20;

   switch ( taskType ) {
      case SIMULATION_TASK :
         // simulator to be executed
         SIMULATOR *simulation;

         // instantiate appropriate simulator
         switch ( simulationType ) {
            case BENCHMARK_SIMULATION :
               simulation = new BENCHMARK;
            break;

            case PHASE_I_SIMULATION :
               simulation = new PHASE_I;
            break;

            case PHASE_II_SIMULATION :
               simulation = new PHASE_II;
            break;

            default :
               printf ( "Simulation Type is not supported!\n" );
               fflush ( stdin );
               getchar ( );
               exit ( 1 );

         };

         // perform simulation
         simulation->Simulate ( );

         // delete dynamics
         delete simulation;
      break;

      case DATA_ANALYSIS_TASK :
         SEQUENCE_ANALYZER *seqAnal;

         // instantiate analyzer
         seqAnal = new SEQUENCE_ANALYZER ( seqFileName, maxCPs, minLength,
                                           algorithmType, cnPower );

         // analyze!
         seqAnal->AnalyzeSequence ( );

         // delete dynamics
         delete seqAnal;
      break;

      default :
         printf ( "My existence has no purpose!\n" );
         fflush ( stdin );
         getchar ( );
         exit ( 1 );

   };

   printf ( "Processing finished...\n" );
   fflush ( stdin );
   getchar ( );

};

