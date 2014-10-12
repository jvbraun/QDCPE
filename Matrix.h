
//
// File: Matrix.h
//
// Purpose: This file defines a template for a MATRIX type; a dynamic
// two-dimensional array of some generic type.
//
//	Author: Ron Braun
//
// History:
//		09/21/98	RKB	Created.
//

#ifndef MATRIXH
#define MATRIXH

#include <malloc.h>

//////////////////
// class MATRIX //
//////////////////

// implements a dynamic two-dimensional array of some generic type

template < class TYPE >
class MATRIX {
   // pointer to the 2D array of data
   TYPE *matrixData;

public:
   // upper bound of first dimension
   long dimension1;
   // upper bound of second dimension
   long dimension2;

   /////////////////
   // constructors

   // construct matrix given two dimensions
   MATRIX ( long dim1, long dim2 ) {
      // initialize dimensions
      dimension1 = dim1;
      dimension2 = dim2;

      // allocate actual matrix of data
      matrixData = ( TYPE * ) malloc ( ( dimension1 + 1 ) *
                                       ( dimension2 + 1 ) *
                                       sizeof ( TYPE ) );

   };

   ////////////
   // methods

   // access a particular element x,y of matrix
   TYPE &Data ( long x, long y ) {
      return *( ( TYPE * ) matrixData + ( x * ( dimension2 + 1 ) + y ) );
   };

   TYPE &operator () ( long x, long y ) {
      return Data ( x, y );
   };

   ////////////////
   // destructors

   // free up memory allocated for matrix
   ~MATRIX ( void ) {
      free ( matrixData );
   };

};

#endif

