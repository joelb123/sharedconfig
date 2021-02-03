/*
  
  \input cnoweb
  
  Program: weighbor                             \hfill\break
  Version: II
  File: weighbor.h                              \hfill\break
  Author: N. D. Socci                           \hfill\break
  
  \title{weighbor.h}
  \job{Weighbor}
  
  \synopsis{\"Weigh"ted neigh\"bor"-joining.
  This program will created a evolutionary tree based on distance
  between sequences. The tree built attempts to minimize the distances
  in it versus the actually distances inputted to the program. The
  output is the tree built along with the branch lengths.}
  
  
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

/* $Id: weighbor.h,v 1.15 2001/07/04 01:38:53 nds Exp $ */

#ifndef __WEIGHBOR_H
#define __WEIGHBOR_H

/*
  __attribute__ is only available in ANSI C compilers and
  since we use it just to suppress a warning get rid of it 
  if __STDC__ is not defined which likely means this is not
  an ANSI C compiler
*/


#ifndef __GNUC__
#define __attribute__(X) 
#endif

#include "matrix.h"
#include "tree.h"

#ifdef MENUIO
#define VERSION_INFO "Weighbor ver 1.2.1 [ 3-Jul-01] MenuIO\n"
#else
#define VERSION_INFO "Weighbor ver 1.2.1 [ 3-Jul-01]\n"
#endif

/*
  \section{Macros}
*/

/* 
   
   \subsection{Distance \|D(a,b)|}
   
   This macro accesses the distance array given indices into the node
   array.  
   
   Note it assumes that the distance matrix is called \|mD| and that the
   array of nodes is called \|nodes|
   
*/

#define INDX(a) (nodes[(a)]->ind)

#define D(a,b) (mD[nodes[(a)]->ind][nodes[(b)]->ind])


/*
  Same for the \|C| array $c_i$ (eq 0.39) 
*/

#define C(a) (vC[nodes[(a)]->ind])

/*
  This macros are for access to the \|mS|, \|mDelB| and \|mDel2B|
  arrays
*/

#define S(a,b)     (mS[nodes[(a)]->ind][nodes[(b)]->ind])
#define DelB(a,b)  (mDelB[nodes[(a)]->ind][nodes[(b)]->ind])
#define Del2B(a,b) (mDel2B[nodes[(a)]->ind][nodes[(b)]->ind])

/*
  
  \section{Global Varibles}
  
  The \|nodes| global varibles is the list of currently active nodes (taxa) 
  which still need to be joined.
  
  The matrix \|mD| is the distance matrix. The indices into this matrix 
  stored in each taxon's node varible. Ie to find the distance between 
  taxon $i$ and $j$ use: 
  
  \qquad\|mD[nodes[i]->ind][nodes[j]->ind]|
  
  The macro \|D(a,b)| is used for this purpose.
  
*/

/*
  
  We need a small positive number $\epsilon$ in cases where all the
  $\sigma$'s turn out to be zero. However, we can not make epsilon too
  small as the algorithm requires that the following be true: 
  $$
  (x+{1\over \epsilon})-{1\over\epsilon} = x
  $$
  when $x$ is of order 1. This means that $\epsilon$ must be greater
  then the machine precision and probably a lot greater. We set
  $\epsilon=10^{-9}$ which is about the square root of most double
  precision and also is roughly the minimal noise for sequence the
  size of the human genome
  
*/

/*
  
  These are the unrenormalized vaules for \|EPSILON| and \|MINB|; ie,
  before we renormalize the vaules for sequence length. \|MINB| is
  used in the calculation of the b's (eq 0.13-0.15 and 0.22-0.24)
  
*/

#define EPS_BARE (1e-9)
#define MINB_BARE (0.0)

#ifdef __WEIGHBOR_C   /* only declare globals once */

NodeT **nodes;    /* array of nodes left to be joined */
MatrixT mD;       /* Distance Matrix */
VectorT vC;       /* Renormalization vector */

/* 
   
   The following global matricies hold the old vaules of the $\Delta
   b$ and other related matrices. They are used in the new $N^3$
   updateing of this matrices by the function \|calcb|. They are
   stored in the format as the \|mD| matrix and are access via the
   \|S|, \|DelB| and \|Del2B| macros which use the node structures to
   get the right indicies.
   
*/

MatrixT mS;        /* Old value of $s_{ij}$ matrix (eq 0.9) */
MatrixT mDelB;     /* $\Delta b_{ij}$ (eq 0.8) */
MatrixT mDel2B;    /* $\Delta^2b_{ij}$ */

/*
  
  Make these varibles globals to make life easier so we do not
  constantly have to pass them around
  
  Also need them hanging around so we can calculate $R(q(i),i,j')$
  after we have decided which pair to join in \|build.c|
  
*/

MatrixT s;       /* $s_{ij}$ eq 0.9 */
MatrixT deltaB;  /* $\Delta b_{ij}$ eq 0.8 */
MatrixT delta2B; /* $\Delta^2 b_{ij}$ */
MatrixT oldDeltaB; /* Save orginal value of \|deltaB| */

/*
  
  These varible keep track of the last joined pair \|(ta,tb)| and the
  newest taxon \|tc|. They are needed in \|calcb| to update the
  $\Delta b$ values in constant time.
  
  N.B. these are really indicies into the \|nodes|. 
  
*/

int ta, tb, tc;

#define DEFAULT_LENGTH (500.0)
double L=DEFAULT_LENGTH;        /* Length of sequences */

/*
  
  Use a default sequence length of 500
  
*/

double EPSILON=(EPS_BARE)/500.0; 
double MINB=(MINB_BARE)/500.0; 

#define DEFAULT_B (4.0)  /* Default number of DNA base types used in the noise and
			    dissimiliarty functions */

double B=DEFAULT_B;
double sigBi; 
double sigBB; 
double sigBiBBLinv; 
double sigBBi;

double maxVar;    /* maximum variance */
double minVar;    /* minimum variance initalized in io.c */

FILE *outfile;    /* File for auxillary output */

int printLevel=0; /* Controls the amount of output */

BooleanT checkQQI=True; /* Turns on checking of $q[q[i]]=i$ */
BooleanT recalcB=False;  /* controls whether $\Delta b$ is recalculated with 
			    the improved $\sigma$'s */

BooleanT oldZflag=False; /* Determines which algorithm is used to calculate
			    $z(i,j)$ */

BooleanT useSigmaBar=False; /* Controls whether the barred version of 
			       $\sigma^2_{i,j;j'}$ is used in the calculation 
			       of R(i,j,j') */

BooleanT useBarValues=True; /* Controls whether the barred version of 
			       $d_{ij}$, $\Delta b$ and $s$ are used in
			       the calculation of R(i,j,jp) */

BooleanT extendedTourn=False; /* Activate the extended tournament */

BooleanT warnFlag=False; /* use to limit the printing of the non-deterministic
			    warning to once per data file */

BooleanT n_Flag=False; /* Change the normalization of $R(i,j,j')$ */
BooleanT w_Flag=False; /* Controls the printing of warnings */
BooleanT x_Flag=False; /* Change normalization of $R(i,j)$ */
BooleanT ZZ_Flag=False;/* Use the new PQ calculation */

#else  /* #ifdef __WEIGHBOR_C */

extern NodeT **nodes;    /* array of nodes left to be joined */
extern MatrixT mD;        /* Distance Matrix */
extern VectorT vC;

extern MatrixT mS;        /* Old value of $s_{ij}$ matrix (eq 0.9) */
extern MatrixT mDelB;     /* $\Delta b_{ij}$ (eq 0.8) */
extern MatrixT mDel2B;    /* $\Delta^2b_{ij}$ */

extern MatrixT s;       /* $s_{ij}$ eq 0.9 */
extern MatrixT deltaB;  /* $\Delta b_{ij}$ eq 0.8 */
extern MatrixT delta2B; /* $\Delta^2 b_{ij}$ */
extern MatrixT oldDeltaB; /* Save orginal value of \|deltaB| */

extern int ta, tb, tc;

extern double L;
extern double EPSILON;
extern double MINB;

extern double B;
extern double sigBi;   /* (1.0/B)     To optimize \|sigma| functions */
extern double sigBB;   /* ((B-1.0)/B) */
extern double sigBiBBLinv; /* ((B-1.0)/(B*B))/L */
extern double sigBBi;  /* (B/(B-1.0)) */

extern double maxVar;    /* maximum variance */
extern double minVar;    /* minimum variance initalized in io.c */

extern FILE *outfile;

extern int printLevel;

extern BooleanT checkQQI;
extern BooleanT recalcB;
extern BooleanT oldZflag; 
extern BooleanT useSigmaBar;
extern BooleanT useBarValues;
extern BooleanT extendedTourn;
extern BooleanT n_Flag;
extern BooleanT w_Flag;
extern BooleanT x_Flag;
extern BooleanT ZZ_Flag;

extern BooleanT warnFlag;

#endif /* #ifdef __WEIGHBOR_C */

/*
  
  \subsection{Compiler Directives}
  
  Here are the flags to control the compliation of various options
  
  \|CHECK_QQOFI| controls whether the code to check that $q(q(i))=i$
  is used in \"build.c". Set flag to turn this code on.
  
  \|LOGFILE| controls whether extra info is dumped to the weighbor.out
  file
  
*/


#undef CHECK_QQOFI
#define LOGFILE 1

/*
  
  \section{function defs}
  
*/

void calc_b(int N, MatrixT b, MatrixT deltaB, MatrixT delta2B, MatrixT s);
void recalc_b(int N, MatrixT b, MatrixT deltaB, MatrixT delta2B, MatrixT s);
void calc_q(int N, int* q, VectorT R, MatrixT b, 
	    int* q2, VectorT LLR, VectorT Zscore);
double calcR(int N, int i, int j, int jp);
double calcR2(int N, int i, int j, int k, MatrixT b);
void calc_z1(int N, VectorT z, int* q, int i);
void calc_z2(int N, VectorT z, int* q, int* q2, MatrixT b, int i);
double calcZ2(int N, int i, int q, int q2, MatrixT b);
double calcPhi(int N, int i, int ip);

double sigma2_3(int, int, int);
double sigma_na(double, double);
double sigma2_3p(int, int, int, double);
double sigma2t(double);
double sigma2tinv(double);
double sigord(int,int,int,int,double,double,double);

int maxVector(int, VectorT);

double SQR(double);
double DMAX(double, double);
double DMIN(double, double);
double RECT(double);

double expm1(double x);
double log1p(double x);

#endif

/* \endc */ 
