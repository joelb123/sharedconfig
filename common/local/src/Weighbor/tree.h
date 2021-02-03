/*
  
  \input cnoweb
  \input epsf
  
  Program: weighbor                             \hfill\break
  File: tree.h                                  \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 4-Feb-97                             \hfill\break
  
  \title{tree.h}
  \job{Weighbor}
  
  \synopsis{This file contains the various data structures used by the
  weighbor program. In particular those to build the best distance tree}
  
  
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

/* $Id: tree.h,v 1.5 2001/06/22 19:05:45 nds Exp $ */

#ifndef __TREE_H
#define __TREE_H

#include <stdio.h>
#include "matrix.h"

/* \section{NodeT}
   
   \|NodeT| is the basic binary tree node. It will hold not only the
   pointer to its childern but the length corresponding to the branch
   from the node to its immediate ancestor (father node). 
   
   The terminal nodes will hold the name of the taxon's they represent
   (read in from the input file). The non-terminal nodes will have
   \|NULL| pointers here for the name string, this will be used to
   distinguish the terminals from the non-terminal nodes. 
   
   Each non-terminal node will also hold the renormalized
   $\sigma^2_\infty(i\bar{i})$ and $\sigma^2_\infty(i'\bar{i})$. This
   value will be used in the recursive calculation of the renormalized
   weights.
   
*/

typedef double LengthT;   /* Type representing the length of a branch */

typedef struct _NodeT {
  
  char *name;            /* The name of the taxon, \|NULL| for
			    non-terminal */
  
  int ind;      /* This is used to hold the nodes position in the 
		   distance matrix. */
  
  LengthT rho;  /* The length of the branch leading back to the nodes
                   father */
  
  struct _NodeT *child_r, *child_l;  /* Pointers to the childern */
  
  int cind_r, cind_l; /* index numbers in to the \|nodes| array for the
			 two childern */
  
} NodeT;


/*
  
  \section{RootNodeT}
  
  The evolutionary tree is not quite a binary tree. It really is an
  un-rooted tree (graph) of tri-nodes; ie each vertex has 3
  edges. This can be represented as a rooted binary tree but the only
  problem is the root node must be a tri-node. 
  
  \vskip1cm\centerline{\epsfxsize=4in\epsfbox{tree.eps}}
  
  Note: any non-terminal node in the unrooted tree can serve as the
  root. The root node will be determined by how the joining algorithm
  progresses and will be created by the last three unjoined taxa.
  
  \|RootNodeT| is the root node data structure and will have 3 child
  pointers. Since it is a non-terminal node it will not have any name
  field. And since it has no ancestor it does not have a length field.
  
  And once we have the root we do not have to worry about
  renormalizing anymore weights so do not need the \|sigma2inf| values
  either.  
  
*/

typedef struct {
  
  NodeT *child_r, *child_m, *child_l;
  
} RootNodeT;


/*
  
  \section{Function Defs}
  
  Here are the various function definitions for processing the tree.
  
*/

NodeT *createNode();
RootNodeT *createRootNode();
void printTree(FILE*, RootNodeT *);
void deleteTree(RootNodeT *);
void deleteNode(NodeT *);

/*
  
  \section{Tree macros}
  
  \subsection{Operations on node pointers}
  
  These macros take as arguments pointers to the \|NodeT| struct.
  
*/

#define LEAF(n) (!((n)->child_r))   /* test if node \|n| is a leaf */
#define RHO(n)  ( (n)->rho )

#define RIND(n) ( (n)->child_r->ind ) /* return the left and right  */
#define LIND(n) ( (n)->child_l->ind ) /* child node indicies        */ 

#define RIGHT(n) ( (n)->child_r )
#define LEFT(n) ( (n)->child_l )

/*
  
  \subsection{Operations on node index numbers}
  
  These macros take as arguments index numbers into the global
  \|nodes| array.
  
*/


#define SIGMA2INF_R(nidx) ( nodes[(nidx)]->sigma2inf_r )
#define SIGMA2INF_L(nidx) ( nodes[(nidx)]->sigma2inf_l )

#define IND(nidx)  ( nodes[(nidx)]->ind )


/*
  
  This macro access the $c_{i;j}$ array stored in each node
  $j$=\|nidx| must be a index into the \|nodes| array
  
  
*/

#endif /* __TREE_H */

/* \endc */





