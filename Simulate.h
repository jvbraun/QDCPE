
//
// File: Simulate.h
//
// Purpose: This file defines the simulators used to investigate the
//    finite-sample properties of an analysis algorithm.  A simulator
//    generates multinomial test sequences and applies the analysis process
//    to the sequences.  The results are collected and summarized.
//
//	Author: Ron Braun
//
// History:
//		07/31/98	RKB	Created.
//

#ifndef SIMULATEH
#define SIMULATEH


/////////////////////
// class SIMULATOR //
/////////////////////

// this class organizes all of the simulation scenarios developed for
// various analyses

class SIMULATOR {
public:
   /////////////////
   // constructors
   //

   ////////////
   // methods
   //

   // dump simulation code here, elaborate structure later!
   virtual void Simulate ( void ) { };

   ////////////////
   // destructors
   //

};


/////////////////////////////////////
// subclass BENCHMARK of SIMULATOR //
/////////////////////////////////////

// this simulator corresponds to the original simulation code developed
// in fortran by J. Braun.  it is retained as a benchmark for use in code
// development, since it generates a known set of outputs.

class BENCHMARK : public SIMULATOR {
public:
   /////////////////
   // constructors
   //

   ////////////
   // methods
   //

   // dump simulation code here from original fortran program
   void Simulate ( void );

   ////////////////
   // destructors
   //

};


///////////////////////////////////
// subclass PHASE_I of SIMULATOR //
///////////////////////////////////

// this simulator corresponds to the first series of simulation described
// in section 4 or "Multiple Change-Point Fitting Via Quasi-Likelihood,
// With Application To DNA Sequence Seqmentation".

class PHASE_I : public SIMULATOR {
public:
   /////////////////
   // constructors
   //

   ////////////
   // methods
   //

   // perform phase I simulation
   void Simulate ( void );

   ////////////////
   // destructors
   //

};


////////////////////////////////////
// subclass PHASE_II of SIMULATOR //
////////////////////////////////////

// this simulator corresponds to the second series of simulation described
// in section 4 or "Multiple Change-Point Fitting Via Quasi-Likelihood,
// With Application To DNA Sequence Seqmentation".

class PHASE_II : public SIMULATOR {
public:
   /////////////////
   // constructors
   //

   ////////////
   // methods
   //

   // perform phase I simulation
   void Simulate ( void );

   ////////////////
   // destructors
   //

};

#endif

