
//
// File: Random.cpp
//
// Purpose: This file defines the random number generator and utilities
//    used by the simulator.
//
//	Author: Ron Braun
//
// History:
//		09/23/98	RKB	Created.
//

#include <math.h>
#include <stdlib.h>

#include "Random.h"
#include "Globals.h"


///////////////////////////////////
// class RANDOM_NUMBER_GENERATOR //
///////////////////////////////////

////////////
// methods

//**********************************************************************
//*  Numerical Recipes software, used in the spirit of academic research.
//*  This function is a renamed version of ran2 from Numerical Recipes.
//**********************************************************************

double RANDOM_NUMBER_GENERATOR :: GetRandomNumber ( void ) {
   const long IM1 = 2147483563,
              IM2 = 2147483399;
   const double AM = 1.0 / IM1;
   const long IMM1 = IM1 - 1,
              IA1 = 40014,
              IA2 = 40692,
              IQ1 = 53668,
              IQ2 = 52774,
              IR1 = 12211,
              IR2 = 3791,
              NTAB = 32,
              NDIV = 1 + IMM1 / NTAB;
   const double EPS = 1.2e-7,
                RNMX = 1.0 - EPS;

// removed the following statics to allow for multiple independent
// random number generators. in actual application, these should
// probably be added back in to ensure randomness across all generators.
//   static long idum2 = 123456789;
//   static long iy = 0;
//   static long iv[NTAB+1];

   long j, k;

   if ( seed <= 0) {
      seed = MAX ( -seed, long ( 1 ) );
      idum2 = seed;
      for ( j = NTAB + 8; j >= 1; j-- ) {
         k = seed / IQ1;
         seed = IA1 * ( seed - k * IQ1 ) - k * IR1;
         if ( seed < 0 )
            seed = seed + IM1;
         if ( j <= NTAB )
            iv[j] = seed;
      };

      iy = iv[1];

   };

   k = seed / IQ1;
   seed = IA1 * ( seed - k * IQ1 ) - k * IR1;
   if ( seed < 0 )
      seed = seed + IM1;
   k = idum2 / IQ2;
   idum2 = IA2 * ( idum2 - k * IQ2 ) - k * IR2;
   if ( idum2 < 0)
      idum2 = idum2 + IM2;
   j = 1 + iy / NDIV;
   iy = iv[j] - idum2;
   iv[j] = seed;
   if ( iy < 1 )
      iy = iy + IMM1;

   return MIN ( AM * iy, RNMX);

};

//   C  (C) Copr. 1986-92 Numerical Recipes Software *L1^"i5..

//**********************************************************************
//*  Numerical Recipes software, used in the spirit of academic research.
//*  This function is a renamed version of expdev from Numerical Recipes.
//**********************************************************************

double RANDOM_NUMBER_GENERATOR :: ExpDev ( void ) {
   double dum;

   do {
      dum = GetRandomNumber ( );
   } while ( dum == 0.0 );

   return - log ( dum );

};

//   C  (C) Copr. 1986-92 Numerical Recipes Software *L1^"i5..

