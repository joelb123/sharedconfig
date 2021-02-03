/*
  
  \input cnoweb
  
  Program: weighbor                             \hfill\break
  File: io.h                                    \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 6-Feb-97                             \hfill\break
  
  \title{io.h}
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

/* $Id: io.h,v 1.4 2001/06/22 03:47:41 nds Exp $ */

#ifndef __IO_H
#define __IO_H

#include "matrix.h"

FILE * openRead(char *filename);
FILE * openWrite(char *filename);
void printError(char *);
int readPhylipFile(FILE *fp, int *N, MatrixT *D, char **names[]);

#endif /* __IO_H */

/* \endc */
