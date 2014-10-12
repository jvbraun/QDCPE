
//
// File: DNASeqDt.cpp
//
// Purpose: This file defines the interface to a DNA sequence data file.
// This file must be in the FASTA format.  This object encapsulates the data
// file and provides a consistent interface for the analyzer and other
// interested components.
//
//	Author: Ron Braun
//
// History:
//		08/06/98	RKB	Created.
//

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "DNASeqDt.h"
#include "Globals.h"


/////////////////////////////
// class DNA_SEQUENCE_DATA //
/////////////////////////////

/////////////////
// constructors

// create sequence of length N.  sequence data is not populated.

DNA_SEQUENCE_DATA :: DNA_SEQUENCE_DATA ( long N ) {
   alphaLength = 4;
   length = N;
   data = new long[length+1];
};

////////////
// methods

//**********************************************************************
//*  Calculate the sample mean Y(l1,l2] and store in mu parameter.
//*  If l1 and l2 are not specified, the mean of the entire sequence
//*  is computed.
//**********************************************************************

void DNA_SEQUENCE_DATA :: ComputeMean ( double mu[5], long l1, long l2 ) {

   // if l2 is default 0 value, replace it with sequence length.  this gets
   // around the limitation of not being able to use an object attribute
   // in the parameter default.
   if ( l2 == 0 ) l2 = length;

   for ( int j1 = 1; j1 <= 4; j1++ )
      mu[j1] = 0.0;

   for ( int i = l1 + 1; i <= l2; i++ )
      mu[data[i]] = mu[data[i]] + 1;

   for ( int j = 1; j <= 4; j++ )
      mu[j] = mu[j] / ( l2 - l1 );

};

//**********************************************************************
//*  This function finds the quasi-deviance and the Pearson measure
//*  for the multinomial case.  Here the l1 and l2 are the change-points.
//*  If no change-points are specified, the function returns the
//*  quasi-deviance of the entire sequence.
//**********************************************************************

double DNA_SEQUENCE_DATA :: QuasiDeviance ( long l1, long l2 ) {
   // computed quasi-deviance of subsequence (l1, l2] of y
   double q;
   // mean of the subsequence (l1, l2] of Y
   double mu[5];

   // if default 0 for l2 was used, replace it with length.  this gets
   // around the limitation of not being able to use an object attribute
   // in the parameter default.
   if ( l2 == 0 ) l2 = length;

   // Get the mean.
   ComputeMean ( mu, l1, l2 );

   // Calculate the quasi-deviance.
   q = 0.0;

   for ( int i = 1; i <= l2 - l1; i++ ) {
      q = q - 2 * log ( mu[data[l1+i]] );

   };

   return q;

};

////////////////
// destructors

DNA_SEQUENCE_DATA :: ~DNA_SEQUENCE_DATA ( void ) {
      delete[] data;
};


////////////////////////////////////////////////////////
// subclass MULTINOMIAL_SEQUENCE of DNA_SEQUENCE_DATA //
////////////////////////////////////////////////////////

//////////////////////
// utility functions

/* old methodology

//**********************************************************************
//*  This function generates individual multinomial observations.
//*  Note that the p-vector is a probability mass function.  M is the
//*  number of multinomials.  We will be concerned with M=1.
//* Note: M parameter removed since it is dead code. RKB
//**********************************************************************

long MULTINOMIAL_SEQUENCE ::
        mlndev ( double p[5] ) {
   double unidev;

   unidev = randGen.GetRandomNumber ( );

   if ( unidev < p[1] ) {
      return 1;
   } else if ( unidev < p[1] + p[2] ) {
      return 2;
   } else if ( unidev < p[1] + p[2] + p[3] ) {
      return 3;
   } else {
      return 4;
   };

};

//*********************************************************************
//*  Generates change-points from a Poisson process with the given
//*  rate.  We need to truncate the process at the bottom or we will
//*  run into problems!
//*  If RATE=0 then we do not really have any change-points (ha ha);
//*  if RATE>0 then we use a Poisson process with mean number of
//*  change-points equal to RATE.
//*********************************************************************

void MULTINOMIAL_SEQUENCE ::
        CreateNewTau0 ( long RATE ) {
   double total, X;
   int repeatLoop;

   if ( RATE == 0 ) {
      rTrue = 0;
      tau0[0] = 0.0;
      tau0[1] = 1.0;
   } else {
      do {
         repeatLoop = FALSE;

         rTrue = 0;
         tau0[0] = 0.0;
         total = 0.0;

         do {
            X = randGen.ExpDev ( );
            total = total + X / RATE;

            if ( total <= 1 ) {
               rTrue = rTrue + 1;

               if ( rTrue > rMax ) break;

               tau0[rTrue] = total;
            };

         } while ( total <= 1 );

//*  Check to see if we have somehow come up with overwhelmingly many
//*  change-points.  Overwhelmingly >= rMax.
         if ( rTrue > rMax )
            repeatLoop = TRUE;

         else {
            tau0[rTrue+1] = 1.0;

//*  We need to make sure that the minimum width is not too small.
//*  In this case we mean that the change-points should occur more than
//*  10% of the unit interval from each other.
            for ( int r = 1; r <= rTrue + 1; r++ ) {
                if ( ( tau0[r] - tau0[r-1] ) <= .1 ) {
                   repeatLoop = TRUE;
                   break;
                };

            };

         };

      } while ( repeatLoop );

   };  // end if

};

//*********************************************************************
//*  Generate the means.  In this incarnation, the first mean is
//*  generated by uniform(epsilon,1-epsilon)^4, normalized, and
//*  the jumps are chosen to be uniform(-2,2) on the logistic scale.
//*********************************************************************

void MULTINOMIAL_SEQUENCE ::
        CreateNewMu0 ( void ) {

   double unidev, MUtot;

   for ( int j1 = 1; j1 <= alphaLength; j1++ ) {
      unidev = randGen.GetRandomNumber ( );
      (*mu0) ( 1, j1 ) = 0.9 * unidev + 0.05;
   };

   MUtot = (*mu0) ( 1, 1 ) + (*mu0) ( 1, 2 ) +
           (*mu0) ( 1, 3 ) + (*mu0) ( 1, 4 );

   for ( int j2 = 1; j2 <= alphaLength; j2++ ) {
      (*mu0) ( 1, j2 ) = (*mu0) ( 1, j2 ) / MUtot;
   };

   for ( int r = 2; r <= rTrue + 1; r++ ) {
      for ( int j = 1; j <= alphaLength; j++ ) {
         unidev = randGen.GetRandomNumber ( );
         (*mu0) ( r, j ) = log ( (*mu0) ( r-1, j ) /
                                 ( 1 - (*mu0) ( r-1, j ) ) );
         (*mu0) ( r, j ) = (*mu0) ( r-1, j ) + 4.0 * unidev - 2.0;
         (*mu0) ( r, j ) = exp ( (*mu0) ( r, j ) ) /
                           ( 1 + exp ( (*mu0) ( r, j ) ) );
      };

      MUtot = (*mu0) ( r, 1 ) + (*mu0) ( r, 2 ) +
              (*mu0) ( r, 3 ) + (*mu0) ( r, 4 );

      for ( int j3 = 1; j3 <= alphaLength; j3++ ) {
         (*mu0) ( r, j3 ) = (*mu0) ( r, j3 ) / MUtot;
      };

   };

};

*/

//**********************************************************************
//*  This function generates individual multinomial observations.
//*  Note that the p-vector is a probability mass function.
//*
//*  alphabet is the character string which gives the category
//*  labels.  For example, the DNA alphabet is "acgt".
//*
//*  The result of calling MULTINOMIAL_SEQUENCE is a string made
//*  up of letters from alphabet.
//**********************************************************************

long MULTINOMIAL_SEQUENCE :: mlndev ( double p[] ) {

   double unidev = randGen.GetRandomNumber ( );

   // We assume that the quantities of interest are in
   // p[1],...,p[alphabetLength] and alphabet[1],...,alphabet[alphabetLength].
   for ( int alphaPos = 1; alphaPos <= alphaLength - 1; alphaPos++) {
      if ( unidev <= p[alphaPos] ) {
         return alphaPos;

      } else {
         unidev = unidev - p[alphaPos];

	   };

   };

   return alphaLength;

};

//*********************************************************************
//*  Generates change-points from a uniform distribution on the unit
//*  interval.  We use a well-known simulation trick to get the order
//*  statistics of R draws from the uniform.
//*********************************************************************

void MULTINOMIAL_SEQUENCE :: CreateNewTau0 ( long rate,
                                             double minRatio ) {

   boolean repeatLoop;

   // save number of change points
   rTrue = rate;

   // Set the left-hand end-point.
   tau0[0] = 0.0;

   do {
      repeatLoop = FALSE;

      // The vector (V_1,...,V_R), where V_i is the
      // ith cumulative sum of R+1 exponential(1) random variables
      // divided by the total of all of them, is distributed
      // as the order statistics from a draw of R uniform deviates.
      //
      // By including the (R+1)th, we get tau0[RTrue+1]=1 automatically.

      for ( int cpNum = 1; cpNum <= rTrue + 1; cpNum++ ) {
         tau0[cpNum] = randGen.ExpDev ( ) + tau0[cpNum-1];

      };

      // Now check to make sure the change-points respect the minimum
      // segment length.

      for ( int cpNum2 = 1; cpNum2 <= rTrue + 1; cpNum2++ ) {
         tau0[cpNum2] = tau0[cpNum2] / tau0[rTrue+1];

         if ( ( tau0[cpNum2] - tau0[cpNum2-1] ) <= minRatio ) {
            repeatLoop = TRUE;

            break;

         };

      };

   } while ( repeatLoop );

};

//*********************************************************************
//*  Generate the true mean parameters.
//*
//*  We will use uniform [-.4,.4], [-.6,.6], and [-.8,.8] to represent
//*  small, medium, and large respectively.  jumpSize gives the parameter.
//*********************************************************************

void MULTINOMIAL_SEQUENCE :: CreateNewMu0 ( double jumpSize ) {

   double total, temp;

   // Generate the first segment mean parameter values from a set
   // of uniform order statistics.

   (*mu0) (1, 0) = 0;

   for ( int alphaPos = 1; alphaPos <= alphaLength; alphaPos++ ) {
     (*mu0) ( 1, alphaPos ) = (*mu0) ( 1, alphaPos-1 ) + randGen.ExpDev ( );

   };

   total = (*mu0) ( 1, alphaLength ) + randGen.ExpDev ( );

   for ( int alphaPos2 = 1; alphaPos2 <= alphaLength; alphaPos2++) {
	   (*mu0) ( 1, alphaPos2 ) = (*mu0) ( 1, alphaPos2 ) / total;

   };

   // Now generate the successive segment true mean parameter values.

   for ( int r = 2; r <= rTrue + 1; r++ ) {

	   // We allow a jump in each parameter, constrained to small,
	   // medium, or large.  Things are not perfect, since we must
	   // normalize the values, but it is the relative comparison which
	   // is important.

	   for ( int alphaPos3 = 1; alphaPos3 <= alphaLength; alphaPos3++ ) {

		   // First create a uniform [-jumpSize, jumpSize] variate.

		   temp = randGen.GetRandomNumber ( );
		   temp = 2 * temp * jumpSize - jumpSize;

		   // Now transform the previous value onto the real line,
		   // perturb it, and transform it back to the unit interval.

		   (*mu0) ( r, alphaPos3 ) = log ( (*mu0) ( r-1, alphaPos3 ) /
                                         ( 1 - (*mu0) ( r-1, alphaPos3 ) ) );
         (*mu0) ( r, alphaPos3 ) = (*mu0) ( r-1, alphaPos3 ) + temp;
         (*mu0) ( r, alphaPos3 ) = exp ( (*mu0) ( r, alphaPos3 ) ) /
                                   ( 1 + exp ( (*mu0) ( r, alphaPos3 ) ) );

      };

	   // Finally, normalize the values.

     total = 0;

	  for ( int alphaPos4 = 1; alphaPos4 <= alphaLength; alphaPos4++ ) {
		  total = total + (*mu0) ( r, alphaPos4 );

	  };

     for ( int alphaPos5 = 1; alphaPos5 <= alphaLength; alphaPos5++ ) {
        (*mu0) ( r, alphaPos5 ) = (*mu0) ( r, alphaPos5 ) / total;

     };

   };

};

//**********************************************************************
//*  This routine takes the mean step-function defined by TAU0 and MU0
//*  and generates a vector of multinomial variates following
//*  those means.  This vector is used to populate the sequence data.
//**********************************************************************

void MULTINOMIAL_SEQUENCE ::
      GenerateNewSequence ( long rate,
                            double jumpSize,
                            long minSegmentLength ) {

   double *muTemp = new double[alphaLength+1];

   // obtain a set of change-points and the mean function
   // which goes with them
   CreateNewTau0 ( rate, ( (double) minSegmentLength ) / length );
   CreateNewMu0 ( jumpSize );

   for ( int r = 1; r <= rTrue + 1; r++ ) {
      for ( int i = int ( length * tau0[r-1] ) + 1;
            i <= int ( length * tau0[r] );
            i++ ) {
         for ( int j = 1; j <= alphaLength; j++ ) {
            muTemp[j] = (*mu0) ( r, j );
         };

         data[i] = mlndev ( muTemp );

      };

   };

   delete[] muTemp;

};

/////////////////////////////////////////////////////
// subclass DATAFILE_SEQUENCE of DNA_SEQUENCE_DATA //
/////////////////////////////////////////////////////

/////////////////
// constructors

// create datafile sequence from a given datafile object
DATAFILE_SEQUENCE :: DATAFILE_SEQUENCE ( FASTA_SEQUENCE_FILE *seqFile ) :
   DNA_SEQUENCE_DATA ( seqFile->SequenceLength ( ) ) {
   // save pointer to sequence file
   sequenceFile = seqFile;

   // construct sequence data from file
   GenerateNewSequence ( );

};

////////////
// methods

// takes the genetics sequence datafile object and generates a
// sequence from it
void DATAFILE_SEQUENCE :: GenerateNewSequence ( void ) {
   sequenceFile->ResetSequenceFile ( );

   // add each element of file to sequence
   for ( int i = 1; i <= sequenceFile->SequenceLength ( ); i++ ) {
      switch ( sequenceFile->GetNextSymbol ( ) ) {
         case 'g' :
         case 'G' :
            data[i] = 1;

         break;

         case 'c' :
         case 'C' :
            data[i] = 2;

         break;

         case 'a' :
         case 'A' :
            data[i] = 3;

         break;

         case 't' :
         case 'T' :
            data[i] = 4;

         break;

         default:
            printf ( "Error!  Illegal character in sequence datafile. \n" );

            fflush ( stdin );
            getchar ( );
            exit ( 1 );

      };

   };

};

////////////////
// destructors


