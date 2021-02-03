/*
  
  \input cnoweb
  
  Program: weighbor                             \hfill\break
  File: calcb.c                                 \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 6-Feb-97                             \hfill\break
  
  \title{calcb.c}
  \job{Weighbor}
  
  Calculate the $b_{i;j}$ matrix.
  
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



#include <math.h>
#include <stdlib.h>

#include "tree.h"
#include "matrix.h"
#include "weighbor.h"

#include <assert.h>

static char version[] __attribute__ ((unused)) = "$Id: calcb.c,v 1.7 2001/07/04 01:38:53 nds Exp $";


/*
  
  \section{calc_b}
  
  Calculate the $b_{i;j}$ matrix eq 0.7 and also calculate and return
  $\Delta b_{ij}$ (eq 0.8), $s_{ij}$ (eq 0.9), and the new terms;
  $\Delta^2b_{ij}$ (eq 0.14), and $\overline{\sigma^2_{ij}}$ (eq 0.15)
  
*/

void
calc_b(int N, MatrixT b, MatrixT deltaB,
       MatrixT delta2B, MatrixT s)
{
  
  int i, j, k;
  
  if(printLevel>3)
    fprintf(outfile, "Previous pair (%d,%d)-->%d\n", 
	    ta, tb, tc);
  
  for(i=0;i<N-1;++i)
    for(j=i+1;j<N;++j) {
      
      s[i][j]       = 0.0;
      deltaB[i][j]  = 0.0;
      delta2B[i][j] = 0.0;
      
      
      /*#ifndef HASH1	*/
      if(S(i,j)<0) 
	/*#endif*/
	{
	  
	  for(k=0;k<N;++k)
	    if(k!=i && k!=j)
	      {
		
		/*
		  epsilon added to handle the cases when the two 
		  $\sigma^2_{ik;j}$ are zero (ie \|sigfac=0|
		*/
		
		double weight=1.0/(sigma2_3(i,k,j)
				   +sigma2_3(j,k,i)+EPSILON);
		
		s[i][j]       += weight;
		deltaB[i][j]  += weight*(D(i,k)-D(j,k));
		delta2B[i][j] += weight*(D(i,k)-D(j,k))*(D(i,k)-D(j,k));
		
		if(s[i][j]<0) {
		  printf("i %d j %d k %d weight %f s_ij %f\n", 
			 i, j, k, weight, s[i][j]);
		  exit(0);
		  assert(s[i][j]>0.0);
		  
		}
	      }
	}
      
      /* #ifndef HASH1 */
      
      else
	{
	  
	  /* 
	     
	     We have already calculated $s_{ij}$ for this $(i,j)$
	     pair.  So now just update by subtracting out the
	     previously joined pair (\|ta|,\|tb|) and adding in the
	     new taxon \|tc|.
	     
	  */
	  
	  double weighta=1.0/(sigma2_3(i,ta,j)
			      +sigma2_3(j,ta,i)+EPSILON);
	  double weightb=1.0/(sigma2_3(i,tb,j)
			      +sigma2_3(j,tb,i)+EPSILON);
	  double weightc=1.0/(sigma2_3(i,tc,j)
			      +sigma2_3(j,tc,i)+EPSILON);
	  
	  
	  s[i][j] = S(i,j) - weighta - weightb + weightc;
	  
	  deltaB[i][j] 
	    = DelB(i,j)
	    - weighta*(D(i,ta)-D(j,ta)) 
	    - weightb*(D(i,tb)-D(j,tb))
	    + weightc*(D(i,tc)-D(j,tc));
	  
	  delta2B[i][j] 
	    = Del2B(i,j)
	    - weighta*(D(i,ta)-D(j,ta))*(D(i,ta)-D(j,ta)) 
	    - weightb*(D(i,tb)-D(j,tb))*(D(i,tb)-D(j,tb))
	    + weightc*(D(i,tc)-D(j,tc))*(D(i,tc)-D(j,tc));
	  
	  
	} 
      /* #endif */
      
      /*
	assert(s[i][j]>0)
	assert( s[i][j]*delta2B[i][j] > -1)
      */
      
      if(s[i][j] <= 0 || delta2B[i][j] <= -1) {
	fprintf(stderr, "%s::%d\n", __FILE__, __LINE__);
	fprintf(stderr, "Round off error in the calculation\n");
	fprintf(stderr, "of s_ij and/or delta^2B_ij. Can\n");
	fprintf(stderr, "not use the N^3 version of\n");
	fprintf(stderr, "weighbor for this tree\n");
	exit(1);
      }
      
    }
  
  for(i=0;i<N-1;++i)
    for(j=i+1;j<N;++j) 
      {
	S(i,j) = S(j,i) = s[j][i] = s[i][j];
	
	DelB(i,j) = deltaB[i][j];
	DelB(j,i) = -deltaB[i][j];
	
	deltaB[i][j]  /= s[i][j];
	deltaB[j][i] = -deltaB[i][j];
	
	Del2B(i,j) = delta2B[i][j];
	Del2B(j,i) = delta2B[i][j];
	
	delta2B[i][j] /= s[i][j];
	delta2B[j][i] = delta2B[i][j];
      }
  
  if(printLevel>2) {
    fprintf(outfile,"Delta b_{ij}=\n");
    for(i=0;i<N;++i) {
      for(j=0;j<N;++j) 
	fprintf(outfile,"%9.4g ", deltaB[i][j]);
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
    
    fprintf(outfile,"Delta^2 b_{ij}=\n");
    for(i=0;i<N;++i) {
      for(j=0;j<N;++j) 
	fprintf(outfile,"%9.4g ", delta2B[i][j]);
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
    
    fprintf(outfile,"s_{ij}=\n");
    for(i=0;i<N;++i) {
      for(j=0;j<N;++j) 
	fprintf(outfile,"%9.4g ", s[i][j]);
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
  }
  
  /*
    Finally, calculate $b_{i;j}$ (eq 0.7)
  */
  
  for(i=0;i<N;++i)
    for(j=0;j<N;++j) 
      {
	
	if( D(i,j) >= fabs(deltaB[i][j]) )
	  
	  b[i][j] = (deltaB[i][j] + D(i,j))/2.0;
	
	else if( D(i,j) < -deltaB[i][j] )
	  
	  b[i][j] = 0.0;
	
	else
	  
	  b[i][j] = D(i,j);
	
      }
  
}

/*
  
  \section{\|recalc_b|}
  
  Redo the calculation of $b_{i;j}$ (eq 0.7) $\Delta b_{ij}$ (eq 0.8),
  $s_{ij}$ (eq 0.9), $\Delta^2b_{ij}$ (eq 0.14) using the new
  $\sigma^{2'}_{ij;k}$ and the old values of $\Delta b_{ij}$.
  
*/

void
recalc_b(int N, MatrixT b, MatrixT deltaB,
	 MatrixT delta2B, MatrixT s)
{
  
  int i, j, k;
  
  MatrixT old_deltaB = matrix(N); /* Tmp matrix to hold old values of 
				     $\delta b_{ij}$ */
  
  setMM(N, deltaB, old_deltaB);
  
  
  for(i=0;i<N;++i)
    for(j=0;j<N;++j) 
      {
	
	s[i][j]       = 0.0;
	deltaB[i][j]  = 0.0;
	delta2B[i][j] = 0.0;
	
	if(i!=j) 
	  {
	    for(k=0;k<N;++k)
	      if(k!=i && k!=j)
		{
		  
		  /*
		    epsilon added to handle the cases when the two 
		    $\sigma^2_{ik;j}$ are zero (ie \|sigfac=0| */
		  
		  /*
		    N.B. We are using the new \|sigma2_3p| function
		  */
		  
		  
		  double weight=1.0/(sigma2_3p(i,k,j,old_deltaB[i][j])
				     +sigma2_3p(j,k,i,old_deltaB[j][i])
				     +EPSILON);
		  
		  s[i][j]       += weight;
		  deltaB[i][j]  += weight*(D(i,k)-D(j,k));
		  delta2B[i][j] += weight*(D(i,k)-D(j,k))*(D(i,k)-D(j,k));
		  
		  if(s[i][j]<0) {
		    printf("Recalc N %d i %d j %d k %d weight=%g sij=%g " 
			   "sigma %g %g [%g %g] deltaB %g %g\n", 
			   N, i, j, k, weight, s[i][j], 
			   sigma2_3p(i,k,j,old_deltaB[i][j]),
			   sigma2_3p(j,k,i,old_deltaB[j][i]),
			   sigma2_3p(k,i,j,old_deltaB[k][j]),
			   sigma2_3p(k,j,i,old_deltaB[k][i]),
			   old_deltaB[i][j], old_deltaB[j][i]);
		    exit(0);
		    assert(s[i][j]>0.0);
		  }
		}
	    
	    deltaB[i][j]  /= s[i][j];
	    delta2B[i][j] /= s[i][j];
	    
	    
	  }
      }
  
  
  if(printLevel>2)
    {
      fprintf(outfile,"Delta b_{ij}=\n");
      for(i=0;i<N;++i) {
	for(j=0;j<N;++j) 
	  fprintf(outfile,"%9.4g ", deltaB[i][j]);
	fprintf(outfile,"\n");
      }
      fprintf(outfile,"\n");
      
      fprintf(outfile,"Delta^2 b_{ij}=\n");
      for(i=0;i<N;++i) {
	for(j=0;j<N;++j) 
	  fprintf(outfile,"%9.4g ", delta2B[i][j]);
	fprintf(outfile,"\n");
      }
      fprintf(outfile,"\n");
      
      fprintf(outfile,"s_{ij}=\n");
      for(i=0;i<N;++i) {
	for(j=0;j<N;++j) 
	  fprintf(outfile,"%9.4g ", s[i][j]);
	fprintf(outfile,"\n");
      }
      fprintf(outfile,"\n");
    }
  
  /*
    Finally, calculate $b_{i;j}$ (eq 0.7)
  */
  
  for(i=0;i<N;++i)
    for(j=0;j<N;++j) 
      {
	
	if( D(i,j) >= fabs(deltaB[i][j]) )
	  
	  b[i][j] = (deltaB[i][j] + D(i,j))/2.0;
	
	else if( D(i,j) < -deltaB[i][j] )
	  
	  b[i][j] = 0.0;
	
	else
	  
	  b[i][j] = D(i,j);
	
      }
  
  
  freeMatrix(old_deltaB);
  
}

/* \endc */ 
