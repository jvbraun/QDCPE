
//
// File: GenFile.h
//
// Purpose: This file provides the interface to a genetics sequence
//    datafile.
//
//	Author: Ron Braun
//
// History:
//		10/25/98 RKB   Created.
//

#ifndef GENFILEH
#define GENFILEH

#include "Globals.h"

///////////////////////////////
// class FASTA_SEQUENCE_FILE //
///////////////////////////////

// dna sequence datafile must be in FASTA format

class FASTA_SEQUENCE_FILE {
   char *_fileName;		// dna sequence datafile name
   FILE *_fp;           // dna sequence datafile in fasta format
   int _seqLen;         // length of sequence

public:
	// constructors
   FASTA_SEQUENCE_FILE ( char *fName );

   // methods
	int OpenDNASequenceFile ( void );
   	// opens dna sequence file, consumes the header information, and
      // sets up the object to read from the beginning of the sequence.
      // returns TRUE if open is successful and file type is recognized
      // as FASTA.  returns FALSE otherwise.
   boolean ResetSequenceFile ( void );
      // returns to the beginning of the sequence.
   int SequenceLength ( void );
      // returns the length of the sequence
   char GetNextSymbol ( void );
      // fetches the next symbol in the sequence.
   void CloseDNASequenceFile ( void );
   	// closes the dna sequence file

   // destructors
   ~FASTA_SEQUENCE_FILE ( void );

};


#endif
