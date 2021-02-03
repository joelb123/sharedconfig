/*
  
  \input cnoweb
  
  Program: weighbor                             \hfill\break
  File: calcq.c                                 \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 6-Feb-97                             \hfill\break
  
  \title{calcq.c}
  \job{Weighbor}
  
  Calculate the $q(i)$ array (eq 0.12). Additionally the $b_{i;j}$ 
  matrix is also calculated.
  
  Note this is the new version which calculates the pairs by
  minimizing the residual of the weighted least squares.
  
  
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

static char version[] __attribute__ ((unused)) = "$Id: calcq.c,v 1.6 2001/06/22 19:05:45 nds Exp $";

#ifndef M_SQRT1_2
#define M_SQRT1_2   0.70710678118654752440      /* 1/sqrt(2) */
#endif


/*
  
  \section{$R(i,j,j')$}
  
  Calculate the cheaper $R$ function to determine $q(i)$ by a
  tournament among all $j$'s
  
*/

double
calcR(int N, int i, int j, int jp)
{
  
  /* modifications by billb for 0.2.8.5 */
  double lnerfcj, lnerfcjp;
  double argj, argjp;
  double Dbpp_jjp, Dbpp_jpj;
  double Dbpp_ijp, Dbpp_ij;
  double Dbpp_jjpB, Dbpp_jpjB; /*B at end means "bar" or $\overline x$*/
  double spp_jjp, spp_jpj;
  double spp_ijp,spp_ij;
  double sijpjBB, sijjpBB;
  double sijpjB, sijjpB, sjjpiB;
  double norm;
  double normj, normjp;
  double dijpB,dijB;
  
  double bi, bj, bjp;
  
#define SIGMA23P(x,y,z) (sigma2_3p(x,y,z,oldDeltaB[x][z]))
  
  if(recalcB)
    {
      norm = 1.0/(SIGMA23P(j,i,jp)+SIGMA23P(jp,i,j)+EPSILON);
      normj = 1.0/(SIGMA23P(i,j,jp)+SIGMA23P(jp,j,i)+EPSILON);
      normjp = 1.0/(SIGMA23P(i,jp,j)+SIGMA23P(j,jp,i)+EPSILON);
    }
  else
    {
      norm = 1.0/(sigma2_3(j,i,jp)+sigma2_3(jp,i,j)+EPSILON);
      normj = 1.0/(sigma2_3(i,j,jp)+sigma2_3(jp,j,i)+EPSILON);
      normjp = 1.0/(sigma2_3(i,jp,j)+sigma2_3(j,jp,i)+EPSILON);
    }      
  
  if(norm<0.0) {
    fprintf(stderr, "%s::%d\n", __FILE__, __LINE__);
    if(recalcB) 
      fprintf(stderr, "Norm < 0 i=%d j=%d jp=%d norm=%g sigs %g %g\n",
	      i,j,jp,norm,
	      SIGMA23P(i,j,jp),SIGMA23P(i,jp,j));
    else
      fprintf(stderr, "Norm < 0 i=%d j=%d jp=%d norm=%g sigs %g %g\n",
	      i,j,jp,norm,
	      sigma2_3(i,j,jp),sigma2_3(i,jp,j));
    exit(1);
  }
  
#ifdef DEBUG
  if(recalcB)
    fprintf(stdout, "i=%d j=%d jp=%d normj=%g sigs %g %g\n",
	    i,j,jp,normj,
	    SIGMA23P(j,i,jp),SIGMA23P(j,jp,i));
  else
    fprintf(stdout, "i=%d j=%d jp=%d normj=%g sigs %g %g\n",
	    i,j,jp,normj,
	    sigma2_3(j,i,jp),sigma2_3(j,jp,i));
#endif
  
  spp_jjp = (s[j][jp] - norm) + EPSILON;
  spp_jpj = (s[jp][j] - norm) + EPSILON;
  spp_ijp = (s[i][jp] - normj) + EPSILON;
  spp_ij = (s[i][j] - normjp) + EPSILON;
  
#ifdef DEBUG
  fprintf(stdout, "spp_jjp=%g, spp_jpj=%g, spp_ijp=%g, spp_ij=%g\n",
	  spp_jjp, spp_jpj, spp_ijp, spp_ij); 
  fprintf(stdout, "    spp_ijp=%g, s[i][jp]=%g, normj=%g\n",
	  spp_ijp, s[i][jp], normj);
#endif
  
  Dbpp_jjp = 
    (s[j][jp]*deltaB[j][jp] - (D(i,j)-D(i,jp))*norm)/(spp_jjp/*+EPSILON*/);
  
  Dbpp_jpj = 
    (s[jp][j]*deltaB[jp][j] - (D(i,jp)-D(i,j))*norm)/(spp_jpj/*+EPSILON*/);
  
  Dbpp_ijp = 
    (s[i][jp]*deltaB[i][jp] - (D(i,j)-D(jp,j))*normj)/(spp_ijp/*+EPSILON*/);
  Dbpp_ij = 
    (s[i][j]*deltaB[i][j] - (D(i,jp)-D(j,jp))*normjp)/(spp_ij/*+EPSILON*/);
  
  
  if(useSigmaBar)
    {
      sijjpB = sigma2_3(i,j,jp);
      sijpjB = sigma2_3(i,jp,j);
      sjjpiB = sigma2_3(j,jp,i);
    }
  else    
    {
      bj  = DMAX( (D(j, jp)+Dbpp_jjp)/2.0, MINB );
      bjp = DMAX( (D(j, jp)-Dbpp_jjp)/2.0, MINB );
      
      bi  = 
	( 
	 DMAX(D(i,j)-bj, MINB)/(sigma2_3(i,j,jp)+EPSILON)
	 +
	 DMAX(D(i,jp)-bjp, MINB)/(sigma2_3(i,jp,j)+EPSILON)
	 )/(1.0/(sigma2_3(i,j,jp)+EPSILON)+1.0/(sigma2_3(i,jp,j)+EPSILON));
      
      sijjpB = sigma_na(bi+C(i),bj+C(j));
      sijpjB = sigma_na(bi+C(i),bjp+C(jp));
      sjjpiB = sigma_na(bj+C(j),bjp+C(jp));
    }
  
  if(useBarValues)
    {
      sijpjBB = 1.0/(1.0/(sijpjB+1.0/spp_jjp)+1.0/(sjjpiB+1.0/spp_ijp));
      sijjpBB = 1.0/(1.0/(sijjpB+1.0/spp_jjp)+1.0/(sjjpiB+1.0/spp_ij));
      
      dijpB = (D(i,jp)/(sijpjB+1.0/spp_jjp)+D(j,jp)/(sjjpiB+1.0/spp_ijp))*
	sijpjBB;
      dijB = (D(i,j)/(sijjpB+1.0/spp_jpj)+D(j,jp)/(sjjpiB+1.0/spp_ij))*
	sijjpBB;
      
      Dbpp_jjpB = (Dbpp_jjp/(sijpjB+1.0/spp_jjp)+Dbpp_ijp/(sjjpiB+1.0/spp_ijp))*
	sijpjBB;
      Dbpp_jpjB = (Dbpp_jpj/(sijjpB+1.0/spp_jpj)+Dbpp_ij/(sjjpiB+1.0/spp_ij))*
	sijjpBB;
    }
  else
    {
      
      sijpjBB = 1.0/spp_jjp + sijpjB;
      sijjpBB = 1.0/spp_jpj + sijjpB;
      
      dijpB = D(i,jp);
      dijB = D(i,j);
      
      Dbpp_jjpB = Dbpp_jjp;
      Dbpp_jpjB = Dbpp_jpj;
      
    }
  
  argj  = (D(i,j)-Dbpp_jjpB-dijpB)
    /sqrt(2.0*(sijpjBB+sijjpB));
  
  argjp = (D(i,jp)-Dbpp_jpjB-dijB)
    /sqrt(2.0*(sijjpBB+sijpjB));
  
  
#ifdef DEBUG
  fprintf(stdout, "ijjp=%d%d%d, bi=%g, bj=%g, bjp=%g,\n\t"
	  "C(i)=%g, C(j)=%g, C(jp)=%g\n",
	  i,j,jp,bi, bj, bjp, C(i), C(j), C(jp));
  fprintf(stdout, "sijjpB=%g, sijpjB=%g, sjjpiB=%g\n",sijjpB, sijpjB, sjjpiB);
  fprintf(stdout, "sijpjBB=%g sijjpB=%g, sijjpBB=%g sijpjB=%g\n\n",
	  sijpjBB,sijjpB,sijjpBB,sijpjB);
#endif
  
#ifdef OFF
  if(argjp!=-argj) {
    fprintf(stderr, "%s::%d\n", __FILE__, __LINE__);
    
    fprintf(stderr, "FATAL ERROR argjp!=-argj\n");
    fprintf(stderr,"%d%d%d %.16g %.16g    %.16g %.16g  sig %g %g"
	    "bi %g bj %g bjp %g norm %g spp %g %g diff %g\n",
	    i,j,jp,argj, argjp,Dbpp_jjp,Dbpp_jpj ,
	    sigma2_3p(i,j,jp), sigma2_3p(i,jp,j), bi, bj, bjp, norm,
	    spp_jjp, spp_jpj, argjp+argj);
    exit(0);
  }
#endif
  
  
  if(argj>0.0)
    lnerfcj = log(derfcx(argj)) - argj*argj;
  else
    lnerfcj = log(derfc(argj));
  
  if(argjp>0.0)
    lnerfcjp = log(derfcx(argjp)) - argjp*argjp;
  else
    lnerfcjp = log(derfc(argjp));
  
  
#ifdef DEBUG
  fprintf(stdout,">>>%lf, %lf || %lf %lf  %g %g\n"
	  "===========================\n\n",
	  s[i][j]*(delta2B[i][j]-SQR(deltaB[i][j]))
	  -s[i][jp]*(delta2B[i][jp]-SQR(deltaB[i][jp])),
	  -2.0*(lnerfcj-lnerfcjp), lnerfcj, lnerfcjp, argj, argjp);
#endif
  
  
  if(!n_Flag)
    return( (s[i][j]*(delta2B[i][j]-SQR(deltaB[i][j]))
	     -s[i][jp]*(delta2B[i][jp]-SQR(deltaB[i][jp])))/(((double)N)-3.0)
	    -2.0*(lnerfcj-lnerfcjp));
  else {
    
    double Y;
    
    Y = (sigma2_3(i,jp,j) + sigma2_3(j,jp,i) + 1.0/spp_ij)
      /(sigma2_3(i,jp,j) + sigma2_3(j,jp,i)+(((double)N)-3.0)/spp_ij);
    
    return(Y*(s[i][j]*(delta2B[i][j]-SQR(deltaB[i][j]))
	      -s[i][jp]*(delta2B[i][jp]-SQR(deltaB[i][jp])))
	   -2.0*(lnerfcj-lnerfcjp));
  }
}

/*
  
  \section{Calculate $q(i)$}
  
  Calculate the $q(i)$ using equation 0.11
  
  This function calls \|calc_b| to calculate the $b_{i;j}$, $\Delta
  b_{ij}$, $\Delta^2 b_{ij}$ and $s_{ij}$ matrices
  
  We also return the residule $R(i,q(i))$ and matrix $b_{i;j}$.
  
*/

void
calc_q(int N, int *q, VectorT R, MatrixT b, 
       int* q2, VectorT LLR, VectorT Zscore)
{
  
  int i, j, jwin;
  double Rmin = -1.0; /* best R */ 
  
  /* for best speed I think we could omit the following (test it first!)
     2 variables.  Without them different second best values are
     printed with -v, but who cares? */
  double R2ndmin = -1.0; /* 2ndbest R */
  double R3rdmin = -1.0; /* 3rdbest R ; note first 2 might be the same */
  
  VectorT z       = vector(N); /* $z(i,q(i))$ */
  
  /*
    First compute $b_{i;j}$, $\Delta b_{ij}$, $\Delta^2b_{ij}$,
    $s_{ij}$ */
  
  calc_b(N, b, deltaB, delta2B, s);
  
  if(recalcB)
    {
      /* 
	 
	 This is the original $\Delta b_{ij}$. The one before the
	 recalculation using the new sigma's. We need to save this for the
	 proper calculation of $s''_{ij}$ which also uses the new sigma's.
	 
      */
      
      /*
	Save the original value of \|deltaB| for use in calculating $s''$
	and $\Delta b''$
      */
      
      setMM(N,deltaB,oldDeltaB);
      
      /* 
	 Now recompute $b_{i;j}$, $\Delta b_{ij}$, $\Delta^2b_{ij}$,
	 $s_{ij}$ using the previously calculted values of $\Delta b_{ij}$
	 and the new noise function.
      */
      
      recalc_b(N, b, deltaB, delta2B, s);
    }
  
  /* 
     
     First calculate the $q(i)$ array using equation 0.11 to run a
     single elmination tournament among the $j$'s 
     
  */
  
  
  if(printLevel>3)
    fprintf(outfile,"args:\n");
  
  for(i=0;i<N;++i)
    {
      
      double minR=0, tmpR;
      int j2win=-1;
      
      if(i!=0)
	jwin = 0;
      else
	jwin = 1;
      for(j=0;j<N;++j)
	if(j!=i && j!=jwin)
	  {
	    if(calcR(N,i,j,jwin)<0.0)
	      jwin = j;
	  }
      
      if(printLevel>3) {
	fprintf(outfile,"basics:\n");
	fprintf(outfile,"%f %f \n",calcR(N,0,2,1),calcR(N,0,2,3));
      }
      
      /*
	Now make sure the \|jwin| is really the winner
      */
      
      
      for(j=0;j<N;++j) 
	if(j!=i && j!=jwin && (((tmpR=calcR(N,i,j,jwin))<minR) || j2win==-1))
	  {
	    minR = tmpR;
	    j2win = j;
	  }
      
      if(minR<0)
	{
	  if(printLevel>1) {
	    fprintf(outfile, "WARNING: non-transitive R(i,j,j'); tree could depend on input order\n");	    
	    fprintf(outfile, "New winner %d (orig winner was %d)\n\n", 
		    j2win, jwin);
	  }
	  
	  if(!warnFlag) {
	    fprintf(stderr, "WARNING: non-transitive R(i,j,j'); tree could depend on input order\n");
	    warnFlag=True;
	  }
	  
	  /* New winner swap, jwin and j2win */
	  
	  q[i] = j2win;
	  q2[i] = jwin;
	}
      else
	{
	  
	  q[i] = jwin;
	  q2[i] = j2win;
	  
	  if(minR==0.0) { 
	    if(printLevel>1) {
	      fprintf(outfile, "WARNING: Tie in R(i,j,j'); tree likely to depend on input order\n");
	      fprintf(outfile, "Tie %d==%d\n\n", j2win, jwin);
	    }
	    
	    if(!warnFlag) {
	      fprintf(stderr, "WARNING: Tie in R(i,j,j'); tree likely to depend on input order\n");
	      warnFlag=True;
	    }
	  }
	  
	  
	}
      
      if(printLevel>2)
	fprintf(outfile, "i=%d q[i]=%d q2[i]=%d minR=%f\n",
		i, q[i], q2[i], minR);
      
      
      /*
	
	The \|LLR| array holds the 3 index R values for the winning
	pair and the second best taxon.
	
      */
      
      
      LLR[i] = -0.5*calcR(N,i,q[i],q2[i]);
      
    }
  
  
  /* 
     
     now calculate $z(i,q(i))$ and then the full residule (eq 0.10)
     
     Two different methods for calculating $z(i,j)$ choose the old one 
     if the \|oldZflag| is set else use the newer ($N^2$) method.
     
  */
  
  
  for(i=0;i<N;++i)
    {
      double arg;
      double lnerfc;
      
      if(!x_Flag) {
	R[i] = 
	  s[i][q[i]]*(delta2B[i][q[i]]-SQR(deltaB[i][q[i]]))/(((double)N)-3.0) ;
	
      } else {
	R[i] = 
	  2.0*s[i][q[i]]*(delta2B[i][q[i]]-SQR(deltaB[i][q[i]]))
	  /(((double)N)-2.0) ;
      }
      
      /* only do z if R looks ok so far */
      if(R[i] < R3rdmin || R3rdmin==-1.0)
	{
	  if(!oldZflag)
	    calc_z2(N, z, q, q2, b, i);
	  else
	    calc_z1(N, z, q, i);
	  /*
	    
	    Finally calculate $R(i,q(i))$ eq 0.10
	    
	  */
	  
	  arg = -z[i]*M_SQRT1_2;
	  
	  
	  /*
	    Save \|Zscore| for use in \|build.c|
	  */
	  Zscore[i] = z[i];
	  
	  if(arg>0.0)
	    lnerfc = log(0.5*derfcx(arg)) - arg*arg;
	  else
	    lnerfc = log(0.5*derfc(arg));
	  
	  R[i] -= 
	    2.0*lnerfc;
	  
	  
	  if(R[i] < R3rdmin || R3rdmin==-1.0) R3rdmin = R[i];
	  if(R[i] < R2ndmin || R2ndmin==-1.0) {R3rdmin = R2ndmin; R2ndmin = R[i];}
	  if(R[i] < Rmin || Rmin==-1.0) {R2ndmin = Rmin; Rmin=R[i];}
	  
	  if(printLevel>2)
	    fprintf(outfile, "z(%d,%d)=%g\n", i, q[i], z[i]); 
	}	  
      
      
    }
  
  if(printLevel>0)
    fprintf(outfile, "\n");
  
  freeVector(z);
  
}

double calcR2(int N, int i, int q, int q2, MatrixT b)
{
  
  double arg = -calcZ2(N, i, q, q2, b)*M_SQRT1_2;
  double lnerfc;
  
  if(arg>0.0)
    lnerfc = log(0.5*derfcx(arg)) - arg*arg;
  else
    lnerfc = log(0.5*derfc(arg));
  
  return(s[i][q]*(delta2B[i][q]-SQR(deltaB[i][q]))/(((double)N)-3.0) - 2.0*lnerfc);
  
}


/* \endc */ 
