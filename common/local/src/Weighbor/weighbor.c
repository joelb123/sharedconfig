/*
  
  \input cnoweb
  
  Program: weighbor                             \hfill\break
  Version: II
  File: weighbor.c                              \hfill\break
  Author: N. D. Socci                           \hfill\break
  
  \title{weighbor.c}
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


/*
  
  \section{Introduction}
  
  The default input file is: \|infile|
  
  The default output file is \|treefile|
  
  The algorithm used is a weighted neighbor-joining one from W. Bruno
  (\'refs'). The weights are computed by the function \|weight(d1,d2)|.
  
  Important note. Array indices are from $0\ldots N-1$ in contrast to
  the first version which went from $1\ldots N$
  
*/

#define __WEIGHBOR_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tree.h"
#include "matrix.h"
#include "io.h"
#include "weighbor.h"

static char version[] __attribute__ ((unused)) = "$Id: weighbor.c,v 1.7 2001/06/22 19:05:45 nds Exp $";

BooleanT expertMode=False;

#ifdef OFF
void
getInput(char *buf, int SIZE)
{
  
  int i=0;
  char ch;
  
  do {
    
    ch=getc(stdin);
    if(ch=='\n'||ch==EOF)
      buf[i++]=(char)NULL;
    else
      buf[i++] = ch;
    
  }while(ch!='\n' && ch!=EOF && i<(SIZE-1));
  
  buf[SIZE-1] = 0;
  
}
#endif 

void
usage()
{
  fprintf(stderr, VERSION_INFO);
  fprintf(stderr, "usage: weighbor [-L num] [-b num] [-i name] [-o name] [-v[v[v]]] [-V]\n");
  fprintf(stderr, "         -L num   set the sequence length to num (default %g)\n", DEFAULT_LENGTH);
  fprintf(stderr, "         -b num   set number of bases (default %g)\n", DEFAULT_B);
  fprintf(stderr, "         -i name  read input from file\n");
  fprintf(stderr, "         -o name  write output to file\n");
  fprintf(stderr, "         -v       turn on logfile and dump level 1 info\n");
  fprintf(stderr, "         -vv      turn on logfile and dump level 2 info\n");
  fprintf(stderr, "         -vvv     turn on logfile and dump level 3 info\n");
  fprintf(stderr, "         -V       print version information\n");
  fprintf(stderr, "\nFor those using this program for research please use the following\n");
  fprintf(stderr, "citation for this program:\n\n");
  fprintf(stderr, "    Bruno, W. J., Socci, N. D. and Halpern, A. L., ``Weighted\n"); 
  fprintf(stderr, "       Neighbor Joining: A Fast Approximation to Maximum-Likelihood\n");
  fprintf(stderr, "       Phylogeny Reconstruction,'' Molecular Biology and Evolution,\n" );
  fprintf(stderr, "       17(1): 189-197, (2000).\n");
  
  if(expertMode) {
    fprintf(stderr, "\n  Advance (expert) options [-XqrzSTeZ]\n");
    fprintf(stderr, "         -X       turn on expert mode (allows use of these options)\n");
    fprintf(stderr, "           n      change normalization of R(i,j,j')\n");
    fprintf(stderr, "           q      turn off q[q[i]]==i checking\n");
    fprintf(stderr, "           r      turn on recalculation of Delta B\n");
    fprintf(stderr, "           z      use the old (N^3) method to calculate z(i,j)\n");
    fprintf(stderr, "           S      use barred version of sigma for R(i,j,j')\n");
    fprintf(stderr, "           T      use unbarred versions of d, Delta B and S for R(i,j,j')\n");
    fprintf(stderr, "           e      extended tournament option\n");
    fprintf(stderr, "           w      turn on non-deterministic warnings\n");
    fprintf(stderr, "           x      change normaliztion of R(i,j)\n");
    fprintf(stderr, "           Z      Calculate PQNew and sigmaNew for z-score\n");
  }
  exit(0);
}

int
main(int argc, char *argv[])
{
  
  int i;
  FILE *fp=stdout, *inputfp=stdin;
  
  int N;            /* Number of taxa  */
  MatrixT D;        /* Distance Matrix */
  char **taxon;     /* names of taxa   */
  BooleanT lengthFlag = False; /* set if the -L option was given */
  
  RootNodeT *tree;
  RootNodeT *buildTree(int, MatrixT, char **);
  
#ifdef MENUIO
  char inputbuf[256];
#endif
  
#ifndef MENUIO
  
  for(i=1;i<argc;++i)
    {
      char option;
      int optidx;
      
      if(argv[i][0]=='-')
	option = argv[i][1];
      else
	usage();
      
      optidx=1;
      while(optidx && (option=argv[i][optidx++])) 
	{
	  
	  switch(option)
	    {
	      
	    case 'V':
	      fprintf(stderr, VERSION_INFO);
	      break;
	      
	    case 'i':
	      ++i;
	      optidx=0;
	      if(i<argc) {
		
		inputfp = fopen(argv[i], "r");
		if(!inputfp)
		  {
		    perror(argv[i]);
		    exit(1);
		  }
	      }
	      else
		usage();
	      
	      break;
	      
	    case 'o':
	      ++i;
	      optidx=0;
	      if(i<argc) {
		
		fp = fopen(argv[i], "w");
		if(!fp)
		  {
		    perror(argv[i]);
		    exit(1);
		  }
	      }
	      else
		usage();
	      
	      break;
	      
	    case 'L':
	      ++i;
	      optidx=0;
	      if(i<argc) {
		
		L = atof(argv[i]);
		if(L<=0.0)
		  {
		    fprintf(stderr, "Must specify a length eg: -L 1000\n");
		    fprintf(stderr, "Note the space after the `L'\n\n");
		    usage();
		  }
		EPSILON = (EPS_BARE)/L;
		MINB    = (MINB_BARE)/L;
		lengthFlag = True;
	      }
	      else 
		usage();
	      
	      break;
	      
	    case 'b':
	      ++i;
	      optidx=0;
	      if(i<argc) {
		
		B = atof(argv[i]);
		if(B<=0.0)
		  {
		    fprintf(stderr, 
			    "Must specify a number greater than zero eg: -b 2.3\n");
		    fprintf(stderr, "Note the space after the `b'\n\n");
		    usage();
		  }
	      }
	      else 
		usage();
	      
	      break;
	      
	    case 'v':
	      ++printLevel;
	      break;
	      
	    case 'X':
	      expertMode = True;
	      optidx=2;
	      while((option=argv[i][optidx++]))
		{
		  printf("Expert Option [%c]\n", option);
		  switch(option)
		    {
		      
		    case 'e':
		      extendedTourn=True;
		      break;
		      
		    case 'S':
		      useSigmaBar=True;
		      break;
		      
		    case 'T':
		      useBarValues=False;
		      break;
		      
		    case 'n':
		      n_Flag = True;
		      break;
		      
		    case 'q':
		      checkQQI = False;
		      break;
		      
		    case 'r':
		      recalcB = True;
		      break;
		      
		    case 'w':
		      w_Flag = True;
		      break;
		      
		    case 'x':
		      x_Flag = True;
		      break;
		      
		    case 'z':
		      oldZflag = True;
		      break;
		      
		    case 'Z':		      
		      oldZflag = False;
		      ZZ_Flag  = True;
		      break;
		      
		    default:
		      usage();
		    }
		}
	      optidx=0;
	      break;
	      
	    default:
	      usage();
	    }
	}
    }
#else
  
  /*
    This code is only compiled if \|MENUIO| is defined. 
    It gets user input from stdin for the various program 
    parameters. 
  */
  
  printf(VERSION_INFO);
  printf("\nNOTE: to accept the default vaule just press ENTER (RETURN) by itself\n\n");
  
  printf("Enter input file (default is infile) : ");
  scanf("%s", inputbuf);
  if(inputbuf[0]!='=') 
    inputfp = fopen(inputbuf, "r");
  else
    inputfp = fopen("infile", "r");
  if(!inputfp) {
    perror("weighbor:Invalid filename:");
    exit(0);
  }
  
#define DEFOUTFILE "treefile"
  
  printf("\nEnter output file (default is %s) : ", DEFOUTFILE);
  scanf("%s", inputbuf);
  if(inputbuf[0]!='=') 
    fp = fopen(inputbuf, "w");
  else
    fp = fopen(DEFOUTFILE, "w");
  if(!fp) {
    perror("weighbor:Invalid filename:");
    exit(0);
  }
  
  printf("\nEnter L, the length of the sequences used to create the\n");
  printf("      distance matrix (default is 500) : ");
  do {
    scanf("%s", inputbuf);
    if(inputbuf[0]!='=')
      L = atof(inputbuf);
    else 
      L = 500;
    if(L<=0.0)
      {
	fprintf(stderr, "Must specify a length greater than zero\n");
      }
    else 
      {
	EPSILON = (EPS_BARE)/L;
	MINB    = (MINB_BARE)/L;
      }
    
  }while(L<=0.0);
  
  printf("\nEnter b, the number of allowed characters in the alphabet\n");
  printf("      (Jukes-Cantor model of DNA would be 4.0, Jukes-Cantor\n");
  printf("      model of proteins would be 20.0; default is 4.0) : ");
  do {
    scanf("%s", inputbuf);
    if(inputbuf[0]!='=')
      B = atof(inputbuf);
    else
      B = 4.0;
    if(B<=1.0)
      {
	fprintf(stderr, 
		"Must specify a number greater than one\n");
      }
  } while(B<=1.0);
  
#endif
  
  
#ifndef MENUIO
  if(!lengthFlag)
    fprintf(stderr, "No length was set. Using the default L=%g\n\n", L);
#endif
  
  if(recalcB && n_Flag) {
    fprintf(stderr, "\nWARNING: Not recommended that you use the -Xr and -Xn flags together\n");
  }
  
  if(B<1.01) {
    
    if(B<1.0) {
      fprintf(stderr, "ERROR: B must be greater than 1.0 now set to %.15g\n", B);
      exit(1);
    }
    
    fprintf(stderr, 
	    "WARNING: weighbor is not expected to work for values of b close to 1.0\n");
    fprintf(stderr,  
	    "You have set b to %.15g\n\n", B);
    
  }
  
  /*
    initialize runtime constants for optimization
  */
  
  sigBi  = (1.0/B);     
  sigBB  = ((B-1.0)/B); 
  sigBBi = (B/(B-1.0));
  sigBiBBLinv = ((B-1.0)/(B*B*L)); 
  
  /*
    
    Read in the Phylip data file. 
    
  */
  
  /* open the auxillary output file */
  
  if(printLevel>0)
    {    
      outfile = fopen("weighbor.out", "w"); 
      if(!outfile)
	printError("weighbor::main::open(outfile)");
    }
  else 
    outfile = NULL;
  
  while(readPhylipFile(inputfp, &N, &D, &taxon))
    {
      
      int i;
      
      if(printLevel>0) {
	fprintf(outfile, "\n============================================");
	if(printLevel>1) {
	  fprintf(outfile, "\nL = %g, b = %g\n", L, B);
	  fprintf(outfile,"\nInput Distance Matrix = \n");
	  for(i=0;i<N;++i) 
	    {
	      int j;
	      for(j=0;j<N;++j)
		fprintf(outfile, "%7.4g ", D[i][j]);
	      fprintf(outfile, "\n");
	    }
	  fprintf(outfile, "\n");
	}
      }
      
      /*
	Reset any global switch and variables
      */
      
      warnFlag=!(w_Flag);
      
      switch(N) {
	
	/* Handle trival cases separately */
	
      case 1: /* 1 taxon; just print the name*/
	
	fprintf(fp, "(%s)\n", taxon[0]);
	if(printLevel>0)
	  fprintf(outfile, "\n(%s)\n\n", taxon[0]);
	
	break;
	
      case 2: /* 2 taxon; trivial tree with each branch $=D_{0,1}/2$ */
	
	fprintf(fp, "(%s:%f,%s:%f)\n", 
		taxon[0], .5*D[0][1], taxon[1], .5*D[0][1]);
	if(printLevel>0)
	  fprintf(outfile, "\n(%s:%f,%s:%f)\n\n", 
		  taxon[0], .5*D[0][1], taxon[1], .5*D[0][1]);
	break;
	
	
      default: /* 3 or more need to build the tree */
	
	tree = buildTree(N, D, taxon);
	printTree(fp, tree);
	if(printLevel>0) {
	  fprintf(outfile, "\n");
	  printTree(outfile, tree);
	}
	deleteTree(tree);
	
	
	/*
	  Now free the memory used by the name strings
	*/
	
	for(i=0;i<N;++i)
	  free(taxon[i]);
	
	break;
	
      }
    }
  
  
  if(printLevel>0)
    fclose(outfile);
  fclose(fp);
  
#ifdef MENUIO
  printf("done\n");
#endif
  
  return(0);
  
}


/* \endc */ 
