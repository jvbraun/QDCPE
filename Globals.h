
//
// File: Globals.h
//
// Purpose: This file defines global types, constants, variables, and
//          utility functions.
//
//	Author: Ron Braun
//
// History:
//		08/06/98	RKB	Created.
//

#ifndef GLOBALSH
#define GLOBALSH

#define MIN(a,b) ((a < b)) ? (a) : (b)
#define MAX(a,b) ((a > b)) ? (a) : (b)

//////////
// types
//

#ifdef __cplusplus
typedef bool boolean;
#endif

#ifndef __cplusplus
// useful boolean type
enum boolean { FALSE, TRUE };
#endif /* __cplusplus */

#include "Random.h"

// available program functions
enum TASKS { SIMULATION_TASK,
             DATA_ANALYSIS_TASK };

// available analysis algorithms
enum ALGORITHMS { BRUTE_FORCE_FAST_ALGORITHM,
                  AUGER_LAWRENCE_FAST_ALGORITHM,
                  AUGER_LAWRENCE_ALGORITHM,
                  MULTI_START_FAST_ALGORITHM,
                  MULTI_START_ALGORITHM,
		  ORDERED_SPLIT_ALGORITHM };

// available simulation scenarios
enum SIMULATIONS { BENCHMARK_SIMULATION,
                   PHASE_I_SIMULATION,
                   PHASE_II_SIMULATION };

//////////////
// constants
//

static const char *DNA_ALPHABET = "aAcCgGtT";

//////////////
// variables
//

// create a random number generator.  we'll use one generator for
// the whole project so that we can reproduce a given run.
static RANDOM_NUMBER_GENERATOR randGen ( -2000 );

// utility functions
//

#endif

