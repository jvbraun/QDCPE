
//
// File: Random.h
//
// Purpose: This file defines the random number generator and utilities
//    used by the simulator.
//
//	Author: Ron Braun
//
// History:
//		09/23/98	RKB	Created.
//

#ifndef RANDOMH
#define RANDOMH

////////////////////////////
// class RANDOM_GENERATOR //
////////////////////////////

class RANDOM_NUMBER_GENERATOR {
   // initial seed for generator
   long seed;

   long idum2;
   long iy;
   long iv[32+1];

public:
   /////////////////
   // constructors

   // initialize generator to seed
   RANDOM_NUMBER_GENERATOR ( long idum ) {
      // save seed
      seed = idum;

      idum2 = 123456789;
      iy = 0;

   };

   ////////////
   // methods

   // generate another random number
   double GetRandomNumber ( void );

   // does something spooky but useful
   double ExpDev ( void );

   ////////////////
   // destructors

};

#endif

