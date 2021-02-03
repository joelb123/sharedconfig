/*
  
  \input cnoweb
  
  Program: weighbor                             \hfill\break
  File: matrix.h                                \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 4-Feb-97                             \hfill\break
  
  \title{matrix.h}
  \job{Weighbor}
  
  \synopsis{Include file for various linear algebra data types: matrix
  and vector} 
  
  
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

/* $Id: matrix.h,v 1.4 2001/06/22 03:47:41 nds Exp $ */

#ifndef __MATRIX_H
#define __MATRIX_H

#include <math.h>

/*
  
  \section{Introduction}
  
  This file contains a rudimentary vector and square matrix datatype.
  
*/

/*
  
  \section{Data Types}
  
  Here are the basic data types: \|MatrixT| which is a 2D dimensional
  array of doubles. 
  It will have the \"C" index convention; ie $i=0\ldots N-1$. 
  The data will be stored in C (row-first)
  ordering. \|VectorT| is a 1D array of doubles again with the index
  going from 0 to $N-1$.
  
  I also include the ever useful BooleanT here.
  
*/

typedef double **MatrixT;
typedef double *VectorT;

typedef enum {False=0, True=-1} BooleanT;


/*
  
  \section{Function Definitions}
  
  \subsection{Allocation and Deallocation}
  
  This functions allocate \& deallocate matrices and vectors. 
  Note all the matrices used by weighbor are square so matrix takes
  just one argument.
  
  
*/

MatrixT matrix(int);
VectorT vector(int);

void freeVector(VectorT);
void freeMatrix(MatrixT);

/*
  
  \subsection{Linear Algebra Routines}
  
  Here are the function def's for any linear algebra routines that
  weighbor will used. They will simply be wrappers around other
  standard library routines
*/

void setMM(int, MatrixT, MatrixT);
void delRowCol(int, MatrixT, int);


/* \section{Printing} */

void printMatrix(int, MatrixT);
void printVector(int, VectorT);

/* \section{Miscellaneous Functions} */

/* 
   
   Complimentary error function: in source file \|calerf.c| Using our
   own code since the clib of gcc only computes to single precision
   even though it returns a double.
   
*/

double derf(double x); 
double derfc(double x); 
double derfcx(double x); 

#endif /* __MATRIX_H */

/* \endc */ 




