
//
// File: GenFile.cpp
//
// Purpose: This file provides the interface to a genetics sequence
//    datafile.
//
//	Author: Ron Braun
//
// History:
//		10/25/98 RKB   Created.
//     1/30/99 RKB   Pulled out use of seqio.c and hardcoded
//                   requirement that file by FASTA format.
//

#include <string.h>
#include <stdio.h>

#include "GenFile.h"
#include "Globals.h"


///////////////////////////////
// class FASTA_SEQUENCE_FILE //
///////////////////////////////

/////////////////
// constructors
FASTA_SEQUENCE_FILE :: FASTA_SEQUENCE_FILE ( char *fName ) {
      _fileName = new char [ strlen ( fName ) + 1 ];
      strcpy ( _fileName, fName );

      _seqLen = 0;

};

////////////
// methods
int FASTA_SEQUENCE_FILE :: OpenDNASequenceFile ( void ) {
   if ( ( _fp = fopen ( _fileName, "r" ) ) == NULL ) {
      printf ( "Couldn't open file %s!\n", _fileName );
      return ( FALSE );
   };

   // consume comment line and make sure it starts with a '>'
   if (! ResetSequenceFile() ) {
      printf ( "File %s is not FASTA format!\n", _fileName);
      return ( FALSE );
   }

   // determine length
   _seqLen = 0;
   while (GetNextSymbol() != EOF)
      _seqLen++;

   if ( _seqLen == 0 ) {
      printf ( "No sequence located in file!\n" );
      return ( FALSE );
   };

   ResetSequenceFile();

   printf ( "Opened file %s of length %d\n", _fileName, _seqLen);

   return TRUE;

};

boolean FASTA_SEQUENCE_FILE :: ResetSequenceFile ( void ) {
   // go to beginning of file
   fseek(_fp, 0, SEEK_SET);

   // consume first (comment) line of FASTA file.
   // the comment line should be <= 80 chars as recommended for
   // FASTA format, but we'll give 'em 511 chars to play with.
   char commentLine[512];
   fgets(commentLine, 512, _fp);

   // make sure comment line starts with '>' and the line was less
   // than 512 characters
   if ((commentLine[0] != '>') ||
       (strlen(commentLine) == 512))
      return FALSE;

   return TRUE;

};

int FASTA_SEQUENCE_FILE :: SequenceLength ( void ) {
   return _seqLen;

};

char FASTA_SEQUENCE_FILE :: GetNextSymbol ( void ) {
   char nextChar;

   while ((nextChar = getc(_fp)) != EOF) {

      // make sure character is valid in DNA alphabet
      if (strchr(DNA_ALPHABET, nextChar) != NULL)
         return nextChar;

   }

   return EOF;

};

void FASTA_SEQUENCE_FILE :: CloseDNASequenceFile ( void ) {
    fclose ( _fp );

};

////////////////
// destructors
FASTA_SEQUENCE_FILE :: ~FASTA_SEQUENCE_FILE ( void ) {
   delete[] _fileName;

};
