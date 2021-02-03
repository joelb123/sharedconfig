/*
  
  \input cnoweb
  
  Program: weighbor                             \hfill\break
  File: io.c                                    \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 6-Feb-97                             \hfill\break
  
  \title{io.c}
  \job{Weighbor}
  
  \synopsis{File I/O Routines}
  
  
*/

/*__*
 *__*
 *__* WEIGHBOR 1.2 [21-Jun-01]
 *__*
 *__* Copyright (c) 1997-2001 
 *__*    Los Alamos National Laboratories
 *__*    The Rockefeller University
 *__*
 *__* This file is part of the WEIGHBOR program that can be used free of
 *__* charge in academic research and teaching. Any commercial use of
 *__* this software requires prior permission and possibly a license. For
 *__* commercial use contact: W. Bruno at billb@t10.lanl.gov Los Alamos
 *__* National Laboratories makes no representations about the
 *__* suitability of this software for any purpose. It is provided 
 *__* "as is" without express or implied warranty
 *__*	
 *__* For those using this program for research please use the following 
 *__* citation for this program:
 *__*
 *__*    Bruno, W. J., Socci, N. D. and Halpern, A. L., ``Weighted 
 *__*        Neighbor Joining: A Fast Approximation to Maximum-Likelihood
 *__*        Phylogeny Reconstruction,'' Molecular Biology and Evolution, 
 *__*        17(1): 189-197, (2000).
 *__*
 *__* 
 *__*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <float.h>
#include "weighbor.h"
#include "matrix.h"

static char version[] __attribute__ ((unused)) = "$Id: io.c,v 1.15 2001/06/22 19:05:45 nds Exp $";

/*
  
  \section{Introduction}
  
  These routines handle the basic file i/o for the weighbor
  program. The the expection of the tree printing routines which are
  in the \|tree| module.
  
  
*/

/*
  
  \section{DBL_EPSILON not defined}
  
  If you got
 
      "*** DBL_EPSILON is not defined in your version ..."

  when trying to compile, read the following:
  
  For setting various parameters in weigbor we need to know the 
  precision of your machine.  This is taken to mean the smallest 
  floating point number greater than 1.0, minus 1.0.  In ANSI
  C this has been defined by the preprocessor variable:
  DBL_EPSILON
  which is normally defined in float.h For some reason either 
  your compiler defines it in a different place, calls it something
  different or it is not defined at all.  If you can find the correct
  value or know what it is, use #define and set DBL_EPSILON below. 
  
  If you can not find it or have no idea what it should be you can 
  uncomment one of the defines below.  The 1.e-09 value will in 
  principle work with any ANSI C implementation; however, this may 
  not provide sufficient precision for weighbor to run on files 
  with many taxa, very large sequence length, or extremly large 
  distances.  The smaller number (9.e-15) will work on any "normal"
  computer, and is probably the better choice unless you're running
  on something strange. 
  
*/

/* #define DBL_EPSILON  1.e-09 */ /* should work on absolutely any machine */
/* #define DBL_EPSILON  9.e-15 */ /* should work better on almost any machine */

#ifndef DBL_EPSILON
#error *** DBL_EPSILON is not defined in your version of the
#error *** C compiler. Read the comments in file io.c 
#error *** for possible solutions
#endif



/*
  
  \subsection{File opening}
  
  Open a file for reading or writing and do a fp pointer
  check. 
  
*/

FILE *
openRead(char *filename)
{
  
  FILE *fp;
  
  fp = fopen(filename, "r");
  if(!fp)
    {
      fprintf(stderr, "%s::%d\n", __FILE__, __LINE__);
      fprintf(stderr, "%s could not be open for reading\n", filename);
      exit(1);
    }
  
  return(fp);
  
}

FILE *
openWrite(char *filename)
{
  
  FILE *fp;
  
  fp = fopen(filename, "w");
  if(!fp)
    {
      fprintf(stderr, "%s::%d\n", __FILE__, __LINE__);
      fprintf(stderr, "%s could not be open for writing\n", filename);
      perror("weighbor::io::openWrite");
      exit(1);
    }
  
  return(fp);
  
}

/*
  
  \subsection{Read Data File}
  
  Read in a Phylip data file. 
  
  The matrix D is reassigned by this routine so it should not be
  pointing to any allocated memory otherwise that memory will be
  lost.
  
  If there is an error in the file the function prints and error and
  exits the program
  
  \item\ fp - is a filepointer to the file to read from. Had better not
  be null
  
  \item\ N - Number of taxz (note this is a pointer to an int)
  
  \item\ D - The distance matrix $N\times N$; it is automatically
  symmetrized by this routine (again a pointer since we change it in
  this function)
  
  \item\ names - An array of strings holding the taxa
  names. Memory dynamically allocated be the routine
  
*/


/*
  print an error and die
*/

void
printError(char *s)
{
  perror("weigbor::io");
  fprintf(stderr, "%s\n", s);
  exit(1);
}


int
readPhylipFile(FILE *fp, int *N, MatrixT *D, char **names[])
{
  
  int i, j;
  double tmp, dMAX = 0;
#define NAMELEN 10
#define BUFSIZE 128
  char buf[BUFSIZE+1];
  
  BooleanT lowerTriag = False;
  BooleanT validFile = True;
  BooleanT inNum = False;
  int numWords = 0;
  
  /*
    Do make sure that the file we are reading is in standard 
    phylip format make sure that the first line of the 
    distance matrix is an integer by itself which will tell
    use the size of the matrix. If it is a float or 
    if there is anything else on the line then flag an invalid 
    format error and stop reading it.
  */
  
  do {
    
    /* This outer do loop eats any blank lines */
    
    if(!fgets(buf, BUFSIZE, fp)) {
      return(0);
    }
    
    i=0;
    while(buf[i] && validFile && numWords<=1) {
      if(isdigit(buf[i])) {
	if(!inNum) {
	  inNum = True;
	  ++numWords;
	} 
      } else if(isspace(buf[i])) {
	if(inNum) {
	  inNum=False;
	}
      } else 
	validFile = False;
      ++i;
    }
    
  } while(validFile && numWords==0);
  
  if(!validFile || numWords>1) {
    fprintf(stderr, "Invalid PHYLIP file format. First line\n");
    fprintf(stderr, "must be contain one integer. The current\n");
    fprintf(stderr, "line is\n\n");
    fprintf(stderr, "%s\n", buf);
    exit(1);
  }
  
  *N = atoi(buf);
  
  if((*N)<1) 
    return(0);
  
  *D = matrix(*N);
  
  /* allocate N string pointers, note we index from $0\ldots N-1$ */
  
  *names = (char **)malloc((*N)*sizeof(char *));
  if(!(*names)) printError("File format error-Out of memory");
  
  for(i=0;i<(*N);++i) {
    int s;
    int limit = lowerTriag?i+1:(*N);
    char c;
    /* 
       read the taxon name. Taxon names are assumed to be in phylip
       format which means the first 10 characters are the name
       including any embedded spaces. Embedded spaces are converted to
       underscores, trailing spaces are removed. 
       
       If there is a newline or tab before 10 characters then take
       whatever is up to the newline or tab as the name. Print a warning.
       
    */
    
    
    /* First read the trailing newline from the last line 
       and any trailing or leading white space*/
    
    while( isspace(c=fgetc(fp)) )
      ;
    
    buf[0] = c;
    for(s=1;s<NAMELEN;++s) {
      c = fgetc(fp);
      
      if(c==0 || c=='\n' || c=='\t') {
	buf[s] = 0;
	break;
      }
      buf[s] = c;
    }
    
    /*
      Have read NAMELEN (10) characters some of which might be 
      space. Now if there if the next char is not a whitespace
      then add it and and more char's up to the first whitespace
      for the name. This gives us a hybrid nameing rule.
      Space are allowed in the first 10 characters; if position 11 is a 
      space then spot; otherwise extend the name until the first 
      white space. So names longer than 10 characters are allowed 
      but we only allow spaces (ASCII 20) in the first 10 positions. 
      
      N.B. Tabs and Newlines always delimit the names so I tab with 
      in the first 10 positions marks the end of the name. 
    */
    
    if(s==NAMELEN) {
      while(!isspace(c=fgetc(fp))) {
	if(s<BUFSIZE) 
	  buf[s++] = c;
      }
    }
    buf[s] = 0;
    
#ifdef DEBUG
    fprintf(stderr, "<<%d>>%d\n[%s]\n", s, strlen(buf),buf);
#endif
    
    --s;
    if(s>=BUFSIZE) s=BUFSIZE-1;
    
    for(;s>=0 && isspace(buf[s]);--s)
      ;
    buf[s+1] = 0;
    s=0;
    while(buf[s]) {
      if(isspace(buf[s]))
	buf[s] = '_';
      ++s;
    }
    
    /* allocate enough memory and move string from tmp buf to name
       array */
    
    (*names)[i] = (char *)malloc((strlen(buf)+1)*sizeof(char));
    if(!((*names)[i])) printError("File format error-Out of memory");
    strcpy((*names)[i], buf);
    
    /* read in the $i^{\rm th}$ row */
    
    for(j=0;j<limit;++j)
      {
	
	if(fscanf(fp, "%lf", &tmp)!=1) {
	  
	  if(j!=i && j!=i+1) {
	    fprintf(stderr, "%s::%d\n", __FILE__, __LINE__);
	    printError("File format error-2");
	  } else {
	    lowerTriag = True;
	    break;
	  }
	} else {
	  if(tmp>dMAX) 
	    dMAX = tmp;
	  (*D)[i][j] = tmp;
	}
      }
    
  }
  
  /* Symmetrize Matrix */
  
  if(lowerTriag == False) {
    for(i=0;i<(*N);++i) 
      for(j=0;j<(*N);++j)
	if(i!=j) (*D)[i][j] = (*D)[j][i] = .5*((*D)[i][j] + (*D)[j][i]);
  } else {
    for(i=0;i<(*N)-1;++i) 
      for(j=i+1;j<(*N);++j)
	(*D)[i][j] = (*D)[j][i];
  }
  
#define LIMITLEN 1
#ifdef LIMITLEN 
  {  
    double maxV00 = 20.0*SQR((B-1)/B);
    double minVar0 = 0.02/SQR(L);
    double MAXD = ((B-1)/(2.0*B))*log((B-1.0)*L)
      +sqrt(maxV00*pow( ((*N)-1)*2*DBL_EPSILON*maxV00/minVar0, (-1.0/3.0)));
    
    if(printLevel>0) {
      fprintf(outfile, "Setting max length to MAXD=%g\n\n", 
	      MAXD);
    }
    if(sigma2tinv(DBL_MAX) < MAXD)
      MAXD = sigma2tinv(DBL_MAX);
    
    for(i=0;i<(*N)-1;++i) 
      for(j=i+1;j<(*N);++j)
	if((*D)[i][j] > MAXD) {
	  dMAX = MAXD;
	  (*D)[i][j] = MAXD;
	  (*D)[j][i] = MAXD;
	}
  }
#endif  
  
  {
    double maxVar0, minVar0;
    maxVar0 = DMIN(sigma2t(dMAX)-2*sigma2t(dMAX/2.0), 20.0*SQR((B-1)/B));
    minVar0 = 0.02/SQR(L);
    
    if(maxVar0 * 2.0 * DBL_EPSILON >= minVar0/((*N)-1)) {
      fprintf(stderr, "%s::%d\n", __FILE__, __LINE__);      
      fprintf(stderr, "Double precision arithmetic is insufficient for the\n");
      fprintf(stderr, "current values of N=%d, L=%g, and d_MAX=%g\n",
	      *N, L, dMAX);
      fprintf(stderr, "Maybe try the N^4 version of weigbor\n");
      exit(0);
    }
    
    minVar = minVar0*pow( ((*N)-1)*2*DBL_EPSILON*maxVar0/minVar0, (2.0/3.0) );
    maxVar = maxVar0*pow( ((*N)-1)*2*DBL_EPSILON*maxVar0/minVar0, (-1.0/3.0));
    
    if(printLevel>0)
      {
	fprintf(outfile,"setting variance bounds to %g ... %g\n",minVar,maxVar);
      }
  }
  
  return(1);
  
}

/* \endc */ 
