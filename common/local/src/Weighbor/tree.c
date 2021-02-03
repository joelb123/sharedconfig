/*
  
  \input cnoweb
  
  Program: weighbor                             \hfill\break
  File: tree.c                                  \hfill\break
  Author: N. D. Socci                           \hfill\break
  Created: 5-Feb-97                             \hfill\break
  
  \title{tree.c}
  \job{Weighbor}
  
  \synopsis{Functions to manipulate the trees built weigbor
  program.}
  
  
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

#include <string.h>
#include <stdlib.h>
#include "tree.h"
#include "matrix.h"
#include "weighbor.h"

static char version[] __attribute__ ((unused)) = "$Id: tree.c,v 1.5 2001/06/22 19:05:45 nds Exp $";

/*
  \section{Memory Managment}
  
  These functions are just wrappers around \|malloc| and make sure that
  no errors have occured and that the data structures are initialized
  to known default values.
  
  Note there are no destructor (free's) as the tree does not need to
  be destroyed.
  
*/

/*
  \subsection{\|createNode|}
*/

NodeT *createNode()
{
  
  NodeT *np;
  
  np = (NodeT *)malloc(sizeof(NodeT));
  if(!np)
    {
      perror("weighbor::tree::createNode");
      exit(1);
    }
  
  np->name=0;
  np->rho = 0.0;
  np->ind = 0;
  np->child_l = NULL;
  np->child_r = NULL;
  np->cind_r = -1;
  np->cind_l = -1;
  
  return(np);
  
}


/*
  \subsection{\|createRootNode|}
*/

RootNodeT *createRootNode()
{
  
  RootNodeT *rp;
  rp = (RootNodeT *)malloc(sizeof(RootNodeT));
  if(!rp)
    {
      perror("weighbor::tree::createRootNode");
      exit(1);
    }
  
  rp->child_l = NULL;
  rp->child_m = NULL;
  rp->child_r = NULL;
  
  return(rp);
  
}

/*
  
  \subsection{Cleanup}
  
  \subsubsection{\|deleteTree|}
  
  Delete all nodes in the tree. Do so recursively from the leaves 
  first.
  
*/

void
deleteTree(RootNodeT *t)
{
  
  /* 
     
     First make sure each sub-tree of each child is empty and then
     delete the child
     
  */
  
  if(t->child_r->child_l || t->child_r->child_r)
    deleteNode(t->child_r);
  free(t->child_r);
  
  if(t->child_m->child_l || t->child_m->child_r)
    deleteNode(t->child_m);
  free(t->child_m);
  
  if(t->child_l->child_l || t->child_l->child_r)
    deleteNode(t->child_l);
  free(t->child_l);
  
  /* 
     Now delete the root
  */
  
  free(t);
  
}

/*
  
  \subsubsection{\|deleteNode|}
  
  Delete the childern of a node by first recusively 
  deleting all decendent nodes.
  
*/

void 
deleteNode(NodeT* n)
{
  
  if(n->child_r->child_r || n->child_r->child_l)
    deleteNode(n->child_r);
  free(n->child_r);
  
  if(n->child_l->child_r || n->child_l->child_l)
    deleteNode(n->child_l);
  free(n->child_l);
  
}


/*
  
  \section{Printing}
  
  \subsection{wrapPrint}
  
  \|wrapPrint| takes care of the actuall printing. I buffers the
  output until it exceeds the line width. It then finds the best line
  break it can and prints the line.
  
  When done wrapPrint must be called with a string argument of NULL to
  flush the remaining buffer.
  
*/

/*
  
  moveMem for stdc libraries that do not have movmem. Note this is
  not a general memmove it will fail if the blocks overlap and
  src<dest!!  This should never happen in this program but check
  anyway.
  
  Will not work if pointers are in different segments. Hopefully it
  will not be compiled on a segmented machine.
  
*/


void
moveMem(char *dest, char *src, int n)
{
  int i;
  
  if(src<dest) {
    fprintf(stderr, "ALGORITHM ERROR\n"
	    "Can not call moveMem with overlapping blocks"
	    "where src is less than destination\n"
	    "weighbor:tree:moveMem\n");
    exit(1);
  }
  
  for(i=0;i<n;++i)
    *(dest++) = *(src++);
  
}

void
wrapPrint(FILE *fp, char *s)
{
  
#define LINEWIDTH 80
  
  static char printBuf[2*LINEWIDTH];
  
  if(!s) { /* s==NULL is the last call, dump whatever is remaining */
    
    fprintf(fp, "%s\n", printBuf);
    printBuf[0]=0;
    return;
    
  } else {
    
    if(strlen(printBuf)+strlen(s)>(2*LINEWIDTH-1))
      {
	fprintf(stderr, "FATAL ERROR--Out of static space\n"
		"weighbor:tree:wrapPrint:: length=%d\n",
		strlen(printBuf)+strlen(s));
	exit(1);
      }
    
    strcat(printBuf, s);
    
    if(strlen(printBuf) >= LINEWIDTH) {
      
      char *breakpoint;
      
      /*
	First try to find a common to break the line
      */
      
      breakpoint = strrchr(printBuf, ',');
      
      if(breakpoint) {
	
	*breakpoint=0;
	fprintf(fp, "%s,\n", printBuf);
	++breakpoint;
	moveMem(printBuf, breakpoint, strlen(breakpoint)+1);
	
      } else {
	
	/*
	  try a ')'
	*/
	
	breakpoint = strrchr(printBuf, ')');
	
	if(breakpoint) {
	  
	  *breakpoint=0;
	  fprintf(fp, "%s)\n", printBuf);
	  ++breakpoint;
	  moveMem(printBuf, breakpoint, strlen(breakpoint)+1);
	  
	} else {
	  
	  /* Break at line break */
	  
	  char tmpc = printBuf[LINEWIDTH];
	  printBuf[LINEWIDTH] = 0;
	  
	  fprintf(fp, "%s\n", printBuf);
	  
	  printBuf[LINEWIDTH]=tmpc;
	  moveMem(printBuf, &(printBuf[LINEWIDTH]), 
		  strlen(&(printBuf[LINEWIDTH]))+1);
	  
	}
      }
      
    }
  }
}

/*
  
  \subsection{printBinaryTree}
  
  This is a local function called by printTree. It does the actual
  traversal of the tree.
  
*/


void printBinaryTree(FILE *fp, NodeT *np)
{
  
  char buf[128];
  
  if(LEAF(np)) {
    sprintf(buf, "%s:%f", np->name, np->rho);
    wrapPrint(fp, buf);
  }
  else {
    wrapPrint(fp, "(");
    printBinaryTree(fp, np->child_l);
    wrapPrint(fp, ",");
    printBinaryTree(fp, np->child_r);
    sprintf(buf, "):%f", np->rho);
    wrapPrint(fp, buf);
  }
  
}

/*
  
  \subsection{printTree}
  
  Function to print a tree. Calls \|printBinaryTree| on all childern of
  the root node.
  
*/

void printTree(FILE *fp, RootNodeT *np)
{
  
  wrapPrint(fp, "(");
  printBinaryTree(fp, np->child_l);
  wrapPrint(fp, ",");
  printBinaryTree(fp, np->child_m);
  wrapPrint(fp, ",");
  printBinaryTree(fp, np->child_r);
  wrapPrint(fp, ");\n");
  wrapPrint(fp, NULL);
}


/* \endc */ 
