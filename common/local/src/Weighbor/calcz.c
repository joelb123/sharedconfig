/*
  
  \input cnoweb
  
  Program: weighbor                             \hfill\break
  File: calcz.c                                 \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 6-Feb-97                             \hfill\break
  
  \title{calcz.c}
  \job{Weighbor}
  
  Calculate the $z(i,q(i))$ array equation 0.23
  
  
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
#include <assert.h>

#include "tree.h"
#include "matrix.h"
#include "weighbor.h"

static char version[] __attribute__ ((unused)) = "$Id: calcz.c,v 1.5 2001/06/22 19:05:45 nds Exp $";

/*
  
  \section{\|calc_z1|}
  
  Calculate the $Z$-score of each $(i, j)$ pair. This is the older
  $N^3$ method and will be called if the \|oldZflag| is set.
  
*/

void
calc_z1(int N, VectorT z, int *q, int i)
{
  
  double PQ(int, int, int, int);
  double sigma2PQ(int, int, int, int);
  VectorT klpair=vector(N);
  
  
  int k, j=q[i];
  double PQ_ij=0.0;           /* $PQ(ij)$ (eq 0.16)            */
  double s_PQ_ij=0.0;         /* $s_{PQ}(ij)$                  */
  double sigma_PQ_ij=0.0;     /* $\sigma^2_{PQ}(ij)$ (eq 0.25) */
  
  
  /*
    reset pair matrix
  */
  
  for(k=0;k<N;++k)
    klpair[k] = -1;
  
  /*
    
    \subsection{Equation 0.16}
    
    $$
    PQ(ij) \equiv {1\over s_{PQ}(ij)}
    \sum_{k} {1\over \sigma^2_{PQ}(ij)}
    \min_l[PQ(ij,kl)]
    $$
    
    \'N.B.' $j\equiv q(i)$, ie we want $PQ(i\,q(i))$
    
    Note we calculate 
    $$
    \sum_{k} {1\over \sigma^2_{PQ}(ij)}
    \min_l[PQ(ij,kl)]
    $$
    and
    $$
    \sum_{k} {1\over \sigma^2_{PQ}(ij)} \equiv s_{PQ}(ij)
    $$
    the normalization factor and also
    $$
    \sigma^2_{PQ}(ij)\equiv{1\over s_{PQ}(ij)}
    +{1\over N-2}\sum_{k\ne i,j}\sigma^2_{ij;k}
    $$
    
    
  */
  
  for(k=0;k<N;++k)
    {
      int l;
      double sigPQ;
      double minPQ;
      int minl;
      
      /*
	$k \ne i$ or $j$
      */
      
      if(k==i || k==j)
	continue;
      
      /*
	get an inital $l\ne i$ or $j=q(i)$ or $k$
      */
      
      l=0;
      while(l==i || l==j || l==k)
	++l;
      
      assert(l<N); /* there had better be at least one l left */
      
      minPQ=PQ(i,j,k,l);
      minl = l;
      
      /*
	
	\subsection{$\min_l[PQ(ij,kl)]$}
	
	Find the minimum over $l$ of $PQ(ij,kl)$
	
      */
      
      for(++l;l<N;++l)
	if(l!=i && l!=j && l!=k && PQ(i,j,k,l)<minPQ)
	  {
	    minPQ=PQ(i,j,k,l);
	    minl=l;
	  }
      
      /*
	Accumulate the second term of eq. 0.25
	$$
	{1\over N-2}\sum_{k\ne i,j}\sigma^2_{ij;k}
	$$
	where $j=q(i)$
      */
      
      if(k!=i && k!=j)
	{
	  sigma_PQ_ij += sigma2_3(i,j,k);
	}
      
      /*
	Keep track of the (k,l) pairs to avoid duplicates; i.e. pair 
	(l,k)
      */
      
      klpair[k] = minl;
      
      /*
	Now check for duplicates
      */
      
      if(minl<k && klpair[minl] == k)
	{
	  continue; /* skip rest of loop; do not add in this term */
	}
      
      /*
	
	\|sigPQ| = $\sigma^2_{PQ}(ij,kl)$. Add $\epsilon$ to
	handle cases when it equals zero
	
      */
      
      sigPQ = sigma2PQ(i,j,k,minl) + EPSILON;
      
      assert(sigPQ>0.0);
      
      PQ_ij += minPQ/sigPQ;
      s_PQ_ij += 1.0/sigPQ;
      
    } /* \|for(k=0;k<N;++k)| */
  
  
  /*
    
    Finish calculation of 
    $$
    \sigma^2_{PQ}(ij) \equiv{1\over s_{PQ}(ij)}
    +{1\over N-2}\sum_{k\ne i,j}\sigma^2_{ij;k}
    $$
    
  */
  
  sigma_PQ_ij = 1.0/s_PQ_ij + (0.25)*(1.0/((double)N-2.0))*sigma_PQ_ij;
  
  /*
    
    And now the $Z$-score (eq 0.23)
    $$
    z(i,q(i)) \equiv PQ(i\,q(i))/\sqrt{\sigma^2_{PQ}(ij)}
    $$
    
    \'N.B.' Remember \|PQ_ij|=$s_{PQ}(ij) PQ(ij)$ so need to divide
    by normalization \|s_PQ_ij|
    
  */
  
  calcPhi(N, i, j);
  
  z[i] = PQ_ij/(s_PQ_ij*sqrt(sigma_PQ_ij));
  
  
  
  freeVector(klpair);
  
}

/*
  
  \section{\|calc_z2|}
  
  This is the newer method that runs in $N^2$ time and which does a
  heuristic min over both $(k,l)$.
  
  The \|b| matrix is passed for used in calculating the new $PQ$ (eq
  38-42). 
  
*/

void
calc_z2(int N, VectorT z, int* q, int* q2, MatrixT b, int i)
{
  
  int j, k, l;
  double PQ(int, int, int, int);
  double sigma2PQ(int, int, int, int);
  
  int minl, mink;
  double mArgl=0, mArgk=0;
  
  j = q[i];
  k = q2[i];
  
  /*
    
    min over $l$
    $$
    \min_l PQ(ij,kl)/\sigma^2_{PQ}(ij,kl)
    $$
    
  */
  
  minl=-1;
  for(l=0;l<N;++l)
    if(l!=i && l!=j && l!=k)
      {
	double argl = PQ(i,j,k,l)
	  /(sqrt(
		 sigma2PQ(i,j,k,l)
		 +(sigma2_3(i,j,k)+sigma2_3(i,j,l))/8.0
		 +EPSILON));
	
	if(minl==-1 || argl<mArgl)
	  {
	    minl = l;
	    mArgl = argl;
	  }
	
      }
  
  /*
    
    min over $k$
    $$
    \min_k PQ(ij,kl)/\sigma^2_{PQ}(ij,kl)
    $$
    
  */
  
  mink=-1;
  for(k=0;k<N;++k)
    if(k!=i && k!=j && k!=minl)
      {
	double argk 
	  = PQ(i,j,k,minl)
	  /(sqrt(
		 sigma2PQ(i,j,k,minl)
		 +(sigma2_3(i,j,k)+sigma2_3(i,j,minl))/8.0
		 +EPSILON));
	
	if(mink==-1 || argk<mArgk)
	  {
	    mink = k;
	    mArgk = argk;
	  }
	
      }
  
  /*
    if \|-Z| flag used then do the new $PQ$ calculation using the 
    already calculated PQ
  */
  
  if(!ZZ_Flag) {
    z[i] = mArgk;
  } else {
    
    double dPQ;
    double PQnew, sigord_ijkl, sigord_ijlk, sigNew;
    
    k = mink;
    l = minl;
    dPQ = PQ(i,j,k,l);
    /* if(dPQ<0.0) dPQ = 0.0; */
    
    sigord_ijkl = sigord(i,j,k,l,dPQ,b[i][j],b[j][i]);
    sigord_ijlk = sigord(i,j,l,k,dPQ,b[i][j],b[j][i]);
    
    if(sigord_ijkl+sigord_ijlk >0.0) { 
      PQnew = 0.5*(
		   (sigord_ijlk*(D(i,k)+D(j,l)) + sigord_ijkl*(D(i,l)+D(j,k)))
		   /(sigord_ijlk+sigord_ijkl) 
		   - D(i,j) - D(k,l) );
      
      sigNew = 0.25*(
		     (sigord_ijlk*sigord_ijkl)/(sigord_ijlk+sigord_ijkl)
		     +0.5*(sigma2_3(k,l,i)+sigma2_3(k,l,j)));
    } else {
      PQnew = 0.5*(-D(i,j)- D(k,l));
      sigNew = 0.125*(sigma2_3(k,l,i)+sigma2_3(k,l,j));
    }
    
    z[i] = PQnew/sqrt(sigNew +0.125*(sigma2_3(i,j,k)+sigma2_3(i,j,l))+EPSILON);
    
  }
  
}

/*
  
  \subsection{\|calcZ2|}
  
  This is the newer method that runs in $N^2$ time and which does a
  heuristic min over both $(k,l)$. This function just calculates Z for
  just one (i,q(i)) pair.
  
*/

double
calcZ2(int N, int i, int q, int q2, MatrixT b)
{
  
  int j, k, l;
  double PQ(int, int, int, int);
  double sigma2PQ(int, int, int, int);
  int minl, mink;
  double mArgl=0, mArgk=0;
  
  j = q;
  k = q2;
  
  /*
    
    min over $l$
    $$
    \min_l PQ(ij,kl)/\sigma^2_{PQ}(ij,kl)
    $$
    
  */
  
  minl=-1;
  for(l=0;l<N;++l)
    if(l!=i && l!=j && l!=k)
      {
	double argl 
	  = PQ(i,j,k,l)
	  /(sqrt(
		 sigma2PQ(i,j,k,l)
		 +(sigma2_3(i,j,k)+sigma2_3(i,j,l))/8.0
		 +EPSILON));
	
	if(minl==-1 || argl<mArgl)
	  {
	    minl = l;
	    mArgl = argl;
	  }
	
      }
  
  /*
    
    min over $k$
    $$
    \min_k PQ(ij,kl)/\sigma^2_{PQ}(ij,kl)
    $$
    
  */
  
  mink=-1;
  for(k=0;k<N;++k)
    if(k!=i && k!=j && k!=minl)
      {
	double argk 
	  = PQ(i,j,k,minl)
	  /(sqrt(
		 sigma2PQ(i,j,k,minl)
		 +(sigma2_3(i,j,k)+sigma2_3(i,j,minl))/8.0
		 +EPSILON));
	
	if(mink==-1 || argk<mArgk)
	  {
	    mink = k;
	    mArgk = argk;
	  }
	
      }
  
  /*
    if \|-Z| flag used then do the new $PQ$ calculation using the 
    already calculated PQ
  */
  
  if(!ZZ_Flag) {
    return(mArgk);
  } else {
    
    double dPQ;
    double PQnew, sigord_ijkl, sigord_ijlk, sigNew;
    
    k = mink;
    l = minl;
    dPQ = PQ(i,j,k,l);
    /*	if(dPQ<0.0) dPQ = 0.0; */
    
    sigord_ijkl = sigord(i,j,k,l,dPQ,b[i][j],b[j][i]);
    sigord_ijlk = sigord(i,j,l,k,dPQ,b[i][j],b[j][i]);
    
    if(sigord_ijkl+sigord_ijlk >0.0) { 
      PQnew = 0.5*(
		   (sigord_ijlk*(D(i,k)+D(j,l)) + sigord_ijkl*(D(i,l)+D(j,k)))
		   /(sigord_ijlk+sigord_ijkl) 
		   - D(i,j) - D(k,l) );
      
      sigNew = 0.25*(
		     (sigord_ijlk*sigord_ijkl)/(sigord_ijlk+sigord_ijkl)
		     +0.5*(sigma2_3(k,l,i)+sigma2_3(k,l,j)));
    } else {
      PQnew = 0.5*(-D(i,j)- D(k,l));
      sigNew = 0.125*(sigma2_3(k,l,i)+sigma2_3(k,l,j));
    }
    
    return(PQnew/sqrt(sigNew +0.125*(sigma2_3(i,j,k)+sigma2_3(i,j,l))+EPSILON));
  }
}




/*
  
  \section{\|PQ(ij;kl)|}
  
  Equation 0.22. Note this equation can be written has the following form
  $$
  {1\over 2}\left[{ {{\alpha\over z_1}+{\beta\over z_2}}
  \over{{1\over z_1}+{1\over z_2}}}-d_{ij}-d_{kl}\right]
  $$
  where
  $\alpha=(d_{ik}+d_{jl})$, $\beta=(d_{il}+d_{jk})$, 
  $z_1=\min(\sigma^2_{ik;j},\sigma^2_{ik;l})
  +\min(\sigma^2_{jl;i},\sigma^2_{jl;k})$ and 
  $z_2=\min(\sigma^2_{il;j},\sigma^2_{il;k})
  +\min(\sigma^2_{jk;i},\sigma^2_{jk;l})$
  
  However if this expression is unstable if $z_1$ or $z_2$ is zero even 
  though the expression is still well defined.
  
  The expression is evaluated in the following form
  $$
  {1\over 2}\left[{{\alpha z_2 + \beta z_1}\over{z_2+z_1}}-d_{ij}-d_{kl}\right]
  $$
  which is finite as long as either is greater than zero. If both are zero 
  then the expression has the following limit
  $$
  {1\over 2}\left[{{\alpha + \beta}\over{2}}-d_{ij}-d_{kl}\right]
  $$
  
*/


double
PQ(int i, int j, int k, int l)
{
  
  double alpha, beta, z1, z2;
  
  alpha = D(i,k)+D(j,l);
  beta  = D(i,l)+D(j,k);
  
  z1 = DMIN(sigma2_3(i,k,j),sigma2_3(i,k,l))
    + DMIN(sigma2_3(j,l,i),sigma2_3(j,l,k));
  
  z2 = DMIN(sigma2_3(i,l,j),sigma2_3(i,l,k))
    + DMIN(sigma2_3(j,k,i),sigma2_3(j,k,l));
  
  if( (z1+z2) == 0.0 )
    return( (((alpha+beta)/2.0)-D(i,j)-D(k,l))/2.0 );
  else
    return( (((alpha*z2+beta*z1)/(z2+z1)) - D(i,j) - D(k,l))/2.0 );
  
}


/*
  
  \section{\|sigma2PQ(ij,kl)|}
  
  As in the case of $PQ(ij, kl)$, $\sigma^2_{PQ}(ij,kl)$ (equation 0.16)
  can also be rewritten in the following form:
  $$
  {1\over 4}\left[ {{z_1 z_2}\over{z_1+z_2}} 
  + (\sigma^2_{kl;i}+\sigma^2_{kl;j})/2 \right]
  $$
  where $z_1$ and $z_2$ are as defined in the comment for the $PQ(ij, kl)$
  function. Again this will be finite as long as either is non-zero. If
  both are then the limit is just zero.
  
*/

double 
sigma2PQ(int i, int j, int k, int l)
{
  
  double z1, z2;
  
  z1 = DMIN(sigma2_3(i,k,j),sigma2_3(i,k,l))
    + DMIN(sigma2_3(j,l,i),sigma2_3(j,l,k));
  
  z2 = DMIN(sigma2_3(i,l,j),sigma2_3(i,l,k))
    + DMIN(sigma2_3(j,k,i),sigma2_3(j,k,l));
  
  if( (z1+z2) == 0.0 )
    return( (sigma2_3(k,l,i)+sigma2_3(k,l,j))/8.0 );
  else
    return( ((z1*z2)/(z1+z2)
	     + (sigma2_3(k,l,i)+sigma2_3(k,l,j))/2.0)/4.0 );
  
}

/*
  \|calcPhi| is used in the renormalization phase and is called from 
  \|build.c|
*/

double
calcPhi(int N, int i, int ip)
{
  
  double PQ(int, int, int, int);
  double sigma2PQ(int, int, int, int);
  VectorT klpair=vector(N);
  int k;
  double PQ_iip=0.0;           /* $PQ(ii')$ (eq 0.16)      */
  double s_PQ_iip=0.0;         /* $s_{PQ}(ii')$            */
  double sigma_iip=0.0;        /* \sigma^2_{ii'} (eq 0.27) */
  
  /*
    reset pair matrix
  */
  
  
  for(k=0;k<N;++k)
    klpair[k] = -1;
  
  
  for(k=0;k<N;++k)
    {
      int l;
      double sigPQ;
      double minPQ;
      int minl;
      
      /*
	$k \ne i$ or $i'$
      */
      
      if(k==i || k==ip)
	continue;
      
      /*
	get an inital $l\ne i$ or $i'=q(i)$ or $k$
      */
      
      l=0;
      while(l==i || l==ip || l==k)
	++l;
      
      assert(l<N); /* there had better be at least one l left */
      
      minPQ=PQ(i,ip,k,l);
      minl = l;
      
      /*
	
	\subsection{$\min_l[PQ(ii',kl)]$}
	
	Find the minimum over $l$ of $PQ(ii',kl)$
	
      */
      
      for(++l;l<N;++l)
	if(l!=i && l!=ip && l!=k && PQ(i,ip,k,l)<minPQ)
	  {
	    minPQ=PQ(i,ip,k,l);
	    minl=l;
	  }
      
      
      /*
	Keep track of the (k,l) pairs to avoid duplicates; i.e. pair 
	(l,k)
      */
      
      klpair[k] = minl;
      
      /*
	Now check for duplicates
      */
      
      if(minl<k && klpair[minl] == k)
	{
	  continue; /* skip rest of loop; do not add in this term */
	}
      
      /*
	
	\|sigPQ| = $\sigma^2_{PQ}(ii',kl)$. Add $\epsilon$ to
	handle cases when it equals zero
	
      */
      
      
      sigPQ = sigma2PQ(i,ip,k,minl) + EPSILON;
      
      assert(sigPQ>0.0);
      
      PQ_iip += minPQ/sigPQ;
      s_PQ_iip += 1.0/sigPQ;
      
      
    } /* \|for(k=0;k<N;++k)| */
  
  
  /*
    
    Calculate $\sigma^2_{ii'}$ equation (0.27)
    
  */
  
  sigma_iip=0.0;
  for(k=0;k<N;++k)
    if(k!=i && k!=ip)
      sigma_iip += sigma2_3(i,ip,k);
  
  sigma_iip /= ((double)N)-2.0;
  
  /*
    
    And now $\phi$ eq (0.26)
    
  */
  
  
  freeVector(klpair);
  
  return( (-PQ_iip/2.0)
	  /( (1.0/(sigma_iip+EPSILON)) + s_PQ_iip/4.0 + EPSILON ) );
  
}



/* \endc */
