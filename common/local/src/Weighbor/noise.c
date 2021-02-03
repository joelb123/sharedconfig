/*
  
  \input cnoweb
  
  Program: weighbor                             \hfill\break
  File: noise.c                                 \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 6-Feb-97                             \hfill\break
  
  \title{noise.c}
  \job{Weighbor}
  
  
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


#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "weighbor.h"

static char version[] __attribute__ ((unused)) = "$Id: noise.c,v 1.10 2001/06/22 19:23:35 nds Exp $";

/*
  
  \section{Dissimilarity---\|D(d)|}
  
  Equation 0.1
  
*/

#define BM1OB ((B-1.0)/B)  /* $(b-1)/b$ */

#define DISS(d) (-BM1OB*(expm1(-(d)/BM1OB)))

/*
  
  \section{Total Variance---\|sigma2t(d)|}
  
  Equation 0.3
  
  $$
  \sigma^2_t(d) = n*e^{2d b/(b-1)}D(d)(1-D(d))
  $$
  
  and its inverse (eq 0.29)
  $$
  \sigma^2_{\rm inv}(x) = {b-1 \over b} \ln\left( 2[x b^2 L + (b-1)^2] \over
  b\sqrt{4 x (b-1) L+(b-1)^2}+ (b-1)(b-2) \right)
  $$
  
  
  
*/

double
sigma2t(double d)
{
  
  /*  return( exp(2.0*d/BM1OB)*DISS(d)*(1-DISS(d))/((double)L) ); */
  
  double sub1; 
  
  sub1 = expm1(d*sigBBi);
  
  /* printf("#%16.12f %17.10g %04X\n", d, sigBB*(sub1*sub1*sigBi + sub1)/(L), 
     CalcHash((unsigned char *)(&d))); */
  
  /* return( sigBB*(sub1*sub1*sigBi + sub1)/(L)); */
  
  return( sigBiBBLinv*sub1*(sub1 + B));
  
  
}


double
sigma2tinv(double x)
{
  double ya;
  double y =   sigBB*log(2.0*(x*B*B*L+(B-1.0)*(B-1.0))/
			 (B*sqrt(4.0*x*(B-1.0)*L+(B-1.0)*(B-1.0))
			  +(B-1.0)*(B-2.0)));
  if(x <= 0.0) return(  0.0);
  if(sqrt(x) > (B-1.0)/(sqrt(L*(B-1.0))*DBL_EPSILON))
    ya = sigBB*(log(x)+log(L*B*B/(B-1.0)))/2.0;
  else if(4*x > sqrt(DBL_EPSILON)*(B-1.0)/L)
    {
      ya =   
	sigBB*log1p((2.0*x*B*B*L+B*((B-1.0)-sqrt((B-1.0)*(B-1.0)+4.0*x*(B-1.0)*L)))/
		    (B*sqrt(4.0*x*(B-1.0)*L+(B-1.0)*(B-1.0))
		     +(B-1.0)*(B-2.0)));
    }
  else
    {
      ya = sigBB*log1p((2.0*x*L*B*(B-1.0))/
		       (B*sqrt(4.0*x*(B-1.0)*L+(B-1.0)*(B-1.0))
			+(B-1.0)*(B-2.0)));
    }
  
  return(ya);
  
}


/*
  
  \section{\|sigma2_3|---$\sigma^2_{ij;k}$}
  
  Non-additive noise for the joining of taxa $i$ and $j$ about taxon $k$.
  See equation 0.4
  
  $$
  \sigma^2_{ij;k} = \sigma^2(d_{ij}) -  \sigma^2(d_{il}) -  \sigma^2(d_{lj})
  $$
  
  The two boolean varibles \|bareI| and \|bareJ| are used to control the
  recursion and prevent the calculation of the renormalized sigma when we
  want the bare sigma for a specific index.
  
*/ 

double
sigma2_3(int i, int j, int k)
{
  
  double d_il, d_lj;
  /*  double args; */
  double sigma;
  
  
  if(i==j || i==k || k==j)
    return(0.0);
  
  if(nodes[i]->ind > nodes[j]->ind) {
    int t = j;
    j = i;
    i = t;
  }
  
#ifdef DEBUG
  
  fprintf(stdout, ">>sigma[%d,%d,%d]\n", i, j, k);
  fprintf(stdout, "indices %d=%d %d=%d %d=%d\n",
	  i, nodes[i]->ind,
	  j, nodes[j]->ind,
	  k, nodes[k]->ind);
  
  fprintf(stdout, "D(%d,%d)=%lg\n", i, j, D(i,j));
  fprintf(stdout, "D(%d,%d)=%lg\n", i, k, D(i,k));
  fprintf(stdout, "D(%d,%d)=%lg\n", j, k, D(j,k));
  
#endif
  
  
  /*
    Equation 0.5
  */
  
  
  if( D(j,k) > (D(i,j)+D(i,k)) ) 
    {
      d_il = 0.0;
      d_lj = D(i,j) - d_il;
      sigma = 0.0;
    }
  else if( D(i,j) > (D(i,k)+D(j,k))) 
    {
      if(D(i,k)<D(j,k))
        d_il = D(i,k);
      else
        d_il = D(i,j)-D(j,k);
      
      d_lj = D(i,j) - d_il;
      sigma = sigma2t(D(i,j)+C(i)+C(j))- sigma2t(d_il+C(i)) - sigma2t(d_lj+C(j));
    }
  else if( D(i,k) > (D(i,j)+D(j,k)) )
    {
      d_il = D(i,j);
      d_lj = D(i,j) - d_il;
      sigma = 0.0;
    }
  else
    {
      d_il = (D(i,j) + D(i,k) - D (j,k))/2;
      d_lj = D(i,j) - d_il;
      sigma = sigma2t(D(i,j)+C(i)+C(j))- sigma2t(d_il+C(i)) - sigma2t(d_lj+C(j));
    }
  
  
  if(sigma < minVar) sigma = minVar;
  if(sigma > maxVar) sigma = maxVar;
  
  /*
    args = (D(i,j)+C(i)+C(j)) - (d_il+C(i)) - (d_lj+C(j));
    
    if(sigma<0.0) {
    if(-sigma<1e-10 || fabs(args)<1e-10){
    sigma=0.0;
    }
    else
    {
    printf("ERROR: sigma2_3(%d,%d,%d)=%g < 0.0\n",
    i,j,k, sigma);
    printf("%g - %g - %g %g\n", D(i,j)+C(i)+C(j),
    d_il+C(i), d_lj+C(j), 
    D(i,j)+C(i)+C(j)-d_il-C(i)-d_lj-C(j));
    sigma=0.0;
    }
    }
  */
  
#ifdef DEBUG
  
  fprintf(stdout, "d_ij=%6lg  %.16lg\n", D(i, j), sigma2t(D(i,j)) );
  fprintf(stdout, "d_il=%6lg  %.16lg\n", d_il, sigma2t(d_il) ); 
  fprintf(stdout, "d_lj=%6lg  %.16lg\n", d_lj, sigma2t(d_lj) );
  fprintf(stdout, "sigma2_3(%d,%d;%d)=%g %016lx\n",
	  i, j, k, sigma, sigma);
  
#endif
  
#ifdef DEBUGOFF
  fprintf(stdout, "sigma2_3(%d,%d;%d)=%g\n",
	  i, j, k, sigma);
#endif
  
#ifdef DEBUGOFF
  printf("%3d %3d %3d | %12.6g %016lx\n", 
	 INDX(i), INDX(j), INDX(k), sigma, sigma);
#endif
  
#ifdef DEBUG9
  {
    unsigned int *dp = (unsigned int *)(&sigma);
    printf("(%3d %3d %3d) [%3d %3d %3d] | %12.6g %08x--%08x | (%8.4f %8.4f %8.4f) [%8.4f %8.4f]\n", 
	   INDX(i), INDX(j), INDX(k), i, j, k, sigma, dp[0], dp[1],
	   D(i,j), D(i,k), D(j,k), C(i), C(j));
  }
#endif
  
  return(sigma);
  
}


double
sigma_na(double X, double Y)
{
  
  double sigma = sigma2t(X+Y)-sigma2t(X)-sigma2t(Y);
  
  if(sigma<0.0) {
    
    /*
     * Some machines seem (notably the Intel ones) to be having
     * roundoff problems so if sigma is small and negative or relatively
     * small (compared to terms) then just set it to zero
     */
    
    if(-sigma<1e-4 || -sigma/(sigma2t(X+Y)+sigma2t(X)+sigma2t(Y)) < 1e-8)
      sigma=0.0;
    else
      {
	printf("ERROR: sigma_na(%g,%g)=%g < 0.0 | X+Y=%g\n",
	       X, Y, sigma, X+Y);
	printf("%12g - %12g - %12g\n", 
	       sigma2t(X+Y),sigma2t(X),sigma2t(Y));
	exit(0);
      }
  }
  
  return(sigma);
  
}

double
sigma2_3p(int i, int k, int j, double deltaBij)
{
  
  double bi, bj, bk;
  
  bi = DMAX(0.5*(D(i,j)+deltaBij),MINB);
  bj = DMAX(0.5*(D(i,j)-deltaBij),MINB);
  
  bk = (
	DMAX(D(i,k)-bi,MINB)/(sigma2_3(i,k,j)+EPSILON)
	+
	DMAX(D(j,k)-bj,MINB)/(sigma2_3(j,k,i)+EPSILON)
	)
    /( 1.0/(sigma2_3(i,k,j)+EPSILON) + 1.0/(sigma2_3(j,k,i)+EPSILON) );
  
  return(sigma_na(bi+C(i),bk+C(k)));
  
}

double
sig4(int i, int k, int j, int l, double dPQ, double diP, double djP)
{
  
  double dQk;
  double sig;
  double sigtot = sigma2_3(i,k,j)+sigma2_3(j,k,i);
  
  dQk = -dPQ;
  
  if(sigtot>0.0) {
    dQk += ( sigma2_3(j,k,i)*(D(i,k)-diP) + sigma2_3(i,k,j)*(D(j,k)-djP) )/sigtot;
  }
  
  /*
    if(dQk<0.0) {
    fprintf(stderr, "dQk==%g sig23[%g %g] ik(%d%d)=%g jk=%g\n", dQk,
    sigma2_3(j,k,i), sigma2_3(i,k,j), i,k,D(i,k), D(j,k));
    dQk = 0.0;
    }
  */
  
  sig=sigma2t(D(i,k)+C(i)+C(k))-sigma2t(diP+C(i))-sigma2t(dPQ)-sigma2t(dQk+C(k));
  
  if(sig<0.0 && outfile) {
    fprintf(outfile, "sig4=%g (%d,%d,%d,%d) %g-%g-%g-%g\n"
	    "PQ=%g iP=%g jP=%g Qk=%g\n d(ik)=%g d(jk)=%g\n", 
	    sig, i,j,k,l,
	    sigma2t(D(i,k)+C(i)+C(k)),
	    sigma2t(diP+C(i)),
	    sigma2t(dPQ),
	    sigma2t(dQk+C(k)),
	    dPQ, diP, djP, dQk,
	    D(i,k), D(j,k));
  }
  
  
  return(sig);
}

double
sigord(int i, int j, int k, int l, double dPQ, double diP, double djP)
{
  
  double sig;
  
  if(dPQ<0.0) return(EPSILON);
  
  sig = 
    sigma2_3(i,k,j) + sigma2_3(i,k,l) - sig4(i,k,j,l,dPQ,diP,djP)
    + sigma2_3(j,l,i) + sigma2_3(j,l,k) - sig4(j,l,i,k,dPQ,djP,diP);
  
  if(sig<0.0 && outfile) {
    fprintf(outfile, "sigord=%g (%d,%d,%d,%d) %g+%g-%g+%g+%g-%g PQ=%g iP=%g jP=%g\n", 
	    sig, i,j,k,l,
	    sigma2_3(i,k,j), sigma2_3(i,k,l), sig4(i,k,j,l,dPQ,diP,djP),
	    sigma2_3(j,l,i), sigma2_3(j,l,k), sig4(j,l,i,k,dPQ,djP,diP),
	    dPQ, diP, djP);
  }
  
  return(sig);
  
}



/* \endc */ 
