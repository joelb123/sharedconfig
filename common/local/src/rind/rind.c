/* rind.c
* This program is copyright 1996 by the Regents of the University of California
* All rights reserved.
* This version is intended for non-commercial use only.
* This program written by Bill Bruno.*/ 
/* "Release 4" added virtual counts to deal with gamma dependence */
/* "Release 5" applies EM to branch lengths as well as probs */
/* "Release 6" tweeked for quicker convergence and uses 2 treescripts*/
/* "Release 6b" changed treescript choice and fixed 2 sequence bug*/
/* "Release 6c" allow local tree topology changes*/
/* "Release 6.4" make file i/o more robust, other user-friendly features */
/* "Release 6.4.1" fixed problem with names of length 10*/
/* "Release 6.4.2" add -v option*/
/* "Release 6.4.3" changed default gap - from J to X; detect spaces in names;
                   made PSCOUNT work; recommend value .65 */
/* "Release 6.4.4" change criterion for topology change attempt       */
/* "Release 6.4.5" allow -P option       */
/* "Release 6.5" put branchlengths in natural units; added exponential
                   prior on distances; change default to weighbor, -s
		   for script */
/* "Release 6.9" corrected Selfchange for all-gap positions */
/* "Release 6.9.1" minor memory improvements  */
/* "release6.9.3 bug fixes for alpha 64 bit processors */
/* "release6.9.4 bug fix for other */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ctype.h>
#include <unistd.h>
#include <sysexits.h>

#define MaxNoSeq 3001     /* maximum number of sequences CHECK THIS VALUE*/
#define MaxLengthSeq 33001  /* maximum length of alignment CHECK THIS VALUE*/
/* It's ok if these are bigger than you need; but can't have
spaces in data, only returns between sequences*/

double Pscount= 0.0;  /* pseudo-count used in oldprobupdate*/

double dist[MaxNoSeq][MaxNoSeq];
char  seq[MaxNoSeq][MaxLengthSeq]; /* store the sequences here*/

/* 0.0 gives maximum likelihood; 1.0 is uniform prior */
/* note: likelihood values ignore the posterior! */
/* .65 worked in Halpern & Bruno, 1998 */

#define DEBUG 0 
#define FINALTOL .05 /* convergence criterion for counts (main loop)*/
#define COUNTTOL .000000005 /* convergence criterion for counts (modelupdate)*/
#define MAINITMAX 20 /*max iterations in converge (main loop)*/
#define MAXUPIT 100 /*max iterations in modelupdate and improvetree*/
#define MINTREEIMPIT 4 /*min iterations in improvetree before giving up*/

#define TREEMAX 6 /* >= 3*/
/*maximum number of trees saved*/
#define TREEMIN 3  /* don't throw away trees unless this many are left */
/*I recommend 3 or more; smaller values may cause trouble */

#define MINDIST 0.0001 /* minimum allowed distance in initial tree */
/* 0 for MINDIST might be dangerous; use something < 1/(20*LengthSeq) */
/* if the tree has a branch length shorter than MINDIST when read 
from treefile, it is increased to MINDIST */
/* 6c : only branches shorter than MINDIST are candidates for 
local rearrangements */
#define MaxSlowSeqs 16 /* max seqs for which to use slowtreescript */

#define LTOL .05 /* likelihood convergence criterion for oldprobupdate*/
#define LIKECUTOFF 7.0 /* exp(-LIKECUTOFF) is not going to effect average*/
#define REALIMPROVEMENT 0.5 /* loglike improved enough to try new tree*/

#define BIGNO 1E+37 /* approximate infinity */
/* the following used only by bbisect and findzero */
#define EPS .000000000000003 /* machine precision (?) 10x smaller gives underflow */
#define TOL .000001 /* convergence criterion for bbisect*/
#define ITMAX 50 /*max iterations in bbisect*/

/* P. Fenimore, April 2009: Increased NAMELENGTH from 10 to 100. */
#define NAMELENGTH      100   /* number of characters max. in species name */
/*extern char *malloc();*/

void early_exit();


typedef char naym[NAMELENGTH];

typedef struct tstruct {
  double value;
  double error;
} Double_with_error;

typedef struct tdata { /*renormdata:*/
  double pcharvect[26]; /*22 of these always needed*/
  double *pnvect[26];  /* may only need one of these*/
  double data[2*MaxNoSeq]; /* storage for pnvects, doubled for B and Z */
  int maxn[26];  /* and one of these*/
  double virtcts[26]; /*newly added for virtual counts; (caused by gamma)*/
  double virtctsp[26]; /*derivative of virtcts*/
} Renormdata;


typedef struct tnode *Treeptr;

typedef struct tnode { /* node in tree */
  int seqno;     /* name of node, if any */
  double length;        /* length from mother to node */
  double changes;       /* expected changes along branch for EM */
  Treeptr left;        /* left daughter address, if any */
  Treeptr right;       /* right daughter address, if any */
  Treeptr mother;      /* mother address, if any */
  char inc;           /* 'y' or 'n' for 'included' */ 
  char mcalc[26];           /* 'y' or 'n' for 'calculated' */
  char rcalc[26];           /* 'y' or 'n' for 'calculated' */
  char lcalc[26];           /* 'y' or 'n' for 'calculated' */
  Renormdata mrdata;
  Renormdata rrdata;
  Renormdata lrdata;
  char side;           /* whether node is l or r or m daughter */
} Treenode;


naym names[MaxNoSeq];
int LengthSeq, NoSeq;  /* length and number of sequences read*/
int Iter; /* no of iterations*/
int Treeit; /* which tree we're on */
int NoTrees; /* how many trees right now */
int Convgd; /* 1 means on last iteration */
char zerodist;  /* whether two seqs have zero distance */

double oldprob[MaxLengthSeq][26], rateequil[MaxLengthSeq],Rateequil,Selfchange,Pscount;
double Eff_Length;  /* effective length of sequence (ignore invar & all gaps)*/
double Posloglike, Virtcts, Virtctsp;

Double_with_error counts[MaxLengthSeq][26];
Double_with_error treecounts[TREEMAX][MaxLengthSeq][26];
Double_with_error probs[MaxLengthSeq][26];
double prevcounts[MaxLengthSeq][26];  /*counts in previous iteration*/
double oldcounts[MaxLengthSeq][26]; 
double virttreecounts[TREEMAX][MaxLengthSeq];
double virttreecountsp[TREEMAX][MaxLengthSeq];
double virtcounts[MaxLengthSeq];
double virtcountsp[MaxLengthSeq];
double entropy[MaxLengthSeq];
double loglikevect[TREEMAX];
double changenorms[2*MaxNoSeq-2];
Double_with_error answer;

int Pos;
char Char;
int npzeros; /*for finzero*/
int Col;
FILE *treefile;
int Verbose,UserTree,TreeScript; /* option to specify using fixed tree */

Treeptr treeloc[TREEMAX][MaxNoSeq],Top[TREEMAX];

FILE *Fopen(name,mode)
char *name, *mode;
{
FILE *tmp;

tmp = fopen(name,mode);

if(tmp == NULL)
  {
    fprintf(stderr, "could not open file ");
    perror(name);
    fprintf(stderr, "exiting \n");
    early_exit(132);
  }
return(tmp);
}

int comparechars(char1,char2)
char char1,char2;
{
if(char1==char2) return(1);

if(char1==' ' && char2=='_') return(1);
if(char1=='_' && char2==' ') return(1);

if('a'<= char1 && char1 <= 'z' && 'A'<=char2 && char2 <='Z' &&
   (char1-char2) == 'a'-'A')
  return(1);
if('a'<= char2 && char2 <= 'z' && 'A'<=char1 && char1 <='Z' &&
   (char2-char1) == 'a'-'A')
  return(1);

return(0);
}

int comparenames(name1,name2)
char *name1, *name2;
{
int i1, i2;
i1 = 0;
i2 = 0;
while((name1[i1] == ' ' || name1[i1]=='_') && i1 < NAMELENGTH-1)
  i1 ++;
while((name2[i2] == ' ' || name2[i2]=='_') && i2 < NAMELENGTH-1)
  i2 ++;
while(comparechars(name1[i1],name2[i2])==1 
      && i1 < NAMELENGTH-1 && i2 < NAMELENGTH-1 && name1[i1] != '\0'
      && name2[i2] != '\0')
  {
    i1 ++;
    i2 ++;
  }

if(name1[i1] == '\0' || (i1 == NAMELENGTH-1 && 
			 (name1[i1] == ' ' || name1[i1] == '_')))
  {
    while(i2 < NAMELENGTH-1 && (name2[i2]== ' ' || name2[i2]== '_'))
      i2 ++;
    if((i2 == NAMELENGTH-1 && (name2[i2]== ' ' || name2[i2] == '_' ))
       || name2[i2]== '\0')
      return(1);
  }
if(name2[i2] == '\0' || (i2 == NAMELENGTH-1 && 
			 (name2[i2] == ' ' || name2[i2] == '_')))
  {
    while(i1 < NAMELENGTH-1 && (name1[i1]== ' ' || name1[i1]== '_'))
      i1 ++;
    if((i1 == NAMELENGTH-1 && (name1[i1]== ' ' || name1[i1] == '_' ))
       || name1[i1]== '\0')
      return(1);
  }

if((i1 == NAMELENGTH-1 || i2 == NAMELENGTH-1 ) && 
   comparechars(name1[i1],name2[i2])==1)
  {
    return(1);
  }
else
  {
    return(0);
  }
}

/* 
   P. Fenimore, April 2009.
   Modified to read table format and retain backward compatibility 
*/
void readdata(char *fn)
{
FILE *fp = NULL;
const char eol[] = {'\n'}, sep[] = {'\t'};
char *line = NULL, *sp1, *sp2, *seq_p;
size_t len;
char tchar;
int tpos,tseq, scanret, error=0,j, ttseq, whitespace, whitespaceerror;
int seqlen[MaxNoSeq], startpos;

/* read sequences with no white space except carriage returns
after each sequence */
/* NAMELENGTH character name followed by white space, then sequences */

if((fp = fopen(fn,"r")) == NULL) {
  //printf("could not open file %s\n", fn);
  //error = 1;
  fprintf(stderr, "Unable to open alignment file.\n");
  perror(fn);
  exit(EX_NOINPUT);
}
for(tseq=0;tseq< MaxNoSeq;tseq++)  
     seqlen[tseq]=0;
	
LengthSeq=0;
tseq=-1;
scanret = 0;
tchar = 'e';
while (((line = fgetln(fp,&len)) != NULL) && (tseq < (MaxNoSeq))) {
  (void) strsep(&line, eol); /* This null-terminates the line of data. */
  tseq++;
  sp1 = NULL;
  sp2 = NULL;
  if ((sp1 = memchr(line, '\t', len)) == NULL) {
    if ((sp2 = memchr(line, ' ', len)) == NULL) {
      fprintf(stderr, "Invalid format in sequence %d\n", tseq);
      exit(EX_DATAERR);
    }
  }
  if (sp1 != NULL) {
    seq_p = strsep(&line, sep);
  }
}
while(scanret != EOF && tseq < MaxNoSeq)
   {
	tseq++;  
	whitespace=0;
	whitespaceerror=0;
	for (j = 0; j < NAMELENGTH && scanret != EOF; j++) 
	  {
	    tchar = getc(fp);
	    if(j==0)
	      {
		while(tchar == '\n')
		  tchar = getc(fp);
	      }
	    if(whitespace==1 && 
	       !(tchar == ' ' || tchar == '\t' || tchar=='\n') )
	      {
		whitespaceerror=1;
		error=1;
	      }
	    if(tchar == ' ' || tchar == '\t' || tchar=='\n') 
	      whitespace=1;


	    if(tchar != EOF)
	      names[tseq][j] = tchar;
	    else
	      scanret = tchar;
	  }
	tpos = -1;
	if(whitespaceerror==1)
	  {
	    fprintf(stderr,"Bad Name: contains internal whitespace\n");
	    fprintf(stderr,"sequence %d, %.10s\n",tseq,names[tseq]);
	    whitespaceerror=0;
	  }

	if(scanret != EOF)
	  {
	    ttseq = 0;
	    while(ttseq< tseq &&comparenames(names[tseq],names[ttseq])==0)
	      ttseq++;
	    if(ttseq < tseq && ttseq > 0 && (seqlen[ttseq-1] <= seqlen[ttseq]))
	      {
		printf("bad data file: either inadvertent name repeat in lines %d and %d\n",tseq+1,ttseq+1);
		printf("or change in order of sequences in different blocks\n");
		early_exit(258);
	      }
	    tseq=ttseq;
	    tpos = seqlen[tseq] -1;
	    startpos = tpos+1;
	  }
if(Verbose>=1)	printf("reading seq %d from pos %d\n",tseq,startpos); 
	while(tpos<MaxLengthSeq && tchar != '\n' && scanret != EOF)
		{
		tpos++;
		scanret = fscanf(fp,"%c",&tchar);
		if(tpos==startpos)
		  {
		    while((tchar == ' ' || tchar == '\t' ||
			  tchar=='\n') && scanret!=EOF)
		      scanret = fscanf(fp,"%c",&tchar);
		  }
			  
		if(tchar !='\n' && scanret != EOF)
		  {
		if(tchar=='-') tchar = 'X';
/* will ignore seq's with gaps */
/* J now means a gap */
		if(tchar=='x') tchar = 'X';
		if(tchar=='$') tchar = 'X';
		if(tchar=='*') tchar = 'X';

		    if((tchar-'A' >25) || (tchar < 'A'))
		      {
			if((tchar-'A' >25) || (tchar < 'A'))
			  {  
			    printf("illegal character! %c should be A-Z\n",tchar);
			    printf("tpos= %d, tseq=%d\n",tpos,tseq);
			    error = 1;
			  }
		      }
		  }
/*		if(tchar=='B' || tchar=='Z')
		  {
		    printf("illegal B or Z! %c not allowed\n",tchar);
		    error = 1;
		  }
*/		/* B and Z have special meanings */
/* B here means N or D; Z here means E or Q.
However, in the probability vectors I store
in B the probability of broken chain of 
inheritance, and in Z the total probability
of all characters */

		if(tchar=='O'|| tchar=='U')
		  {
		    printf("illegal O or U! %c never allowed\n",tchar);
		    error = 1;
		  }
/* J now means a gap */
		if((tpos == MaxLengthSeq) && (tchar != '\n'))
		  {
		    printf("sequence too long! MaxLengthSeq=%d\n",
			     MaxLengthSeq);
		    printf("tseq= %d\n",tseq);

		    while(tchar != '\n' && scanret!=EOF)
		      {
			scanret = fscanf(fp,"%c",&tchar);
			tpos++;
		      }
		    printf("this seq is %d positions long\n",tpos);

		    error = 1;
		  }
		if(tchar != '\n' && scanret != EOF)
		  {
		    seq[tseq][tpos]=tchar;
		    seqlen[tseq]++;
		  }
	      }       /* end while */
	if(LengthSeq == 0) 
	  LengthSeq = tpos;
	if(tpos > LengthSeq && tseq==NoSeq)
	  {
	    fprintf(stderr,"warning: seq's not all same length\n");
	    NoSeq++;
	  }
	if(tpos < LengthSeq && tpos > 0)
	  {
	    fprintf(stderr,"warning: seq's not all same length\n");
	    if(tseq==NoSeq)
	      NoSeq ++;
	  }
	if(tpos == LengthSeq) NoSeq ++;

	if((tseq == MaxNoSeq) && (scanret!= EOF)) 
	  {
	    printf("too many sequences! MaxNoSeq= %d\n",MaxNoSeq);
	    error = 1;
	  }
	if(tchar == '\n') tchar = 'r';
	
      }
LengthSeq=seqlen[0];
for(tseq=0;tseq<NoSeq;tseq++)
{
  if(seqlen[tseq]>LengthSeq)
    LengthSeq=seqlen[tseq];
}
for(tseq=0;tseq<NoSeq;tseq++)
{
  if(seqlen[tseq]!=LengthSeq)
    {
      printf("warning, seq's not same length: LengthSeq=%d, tseq=%d, seqlen[tseq]=%d NoSeq=%d\n",
	     LengthSeq,tseq,seqlen[tseq],NoSeq);
      for(tpos=seqlen[tseq];tpos<LengthSeq;tpos++)
	seq[tseq][tpos]='X';

    }
}

printf("read %d sequences of length %d\n", NoSeq, LengthSeq);
printf("here are first and last sequences\n");
tseq=0;  /*print out first and last seqs*/
for(tpos=0;tpos<LengthSeq;tpos++)
  printf("%c",seq[tseq][tpos]);
printf("\n");
tseq=NoSeq-1;  /*print out first and last seqs*/
for(tpos=0;tpos<LengthSeq;tpos++)
  printf("%c",seq[tseq][tpos]);
printf("\n");

/* treebuilding programs require at least 3 sequences, so: */

if(NoSeq <3) 
{
  for(tseq=NoSeq;tseq<3;tseq++)
    {
      sprintf(names[tseq],"DummySeq%d ",tseq);
      for(tpos=0;tpos<LengthSeq;tpos++)
	{
	  seq[tseq][tpos]=seq[0][tpos];
	}
    }
  NoSeq = 3;
}


fclose(fp);
if(error ==1)
{
  printf("exiting\n");
  early_exit(390);
}

}

double pck(prob,loc)
double prob;       /*check probability is in bounds (allow some error) */
int loc;
{
double retval;
int err=0;

retval = prob;

if(prob > 1.0) 
  {
    if(prob > 1.000001)
      {
	printf("probcheck: at %d prob= %f > 1\n",loc,prob);
	err = 1;
      }
    retval = 1.0;
  }
if(prob < 0.0) 
  {
    printf("probcheck: at %d prob= %f < 0\n",loc,prob);
    err = 1;
    retval = 0.0;
  }
if(prob <= 1.000001 && prob >= 0.0 && err == 0)
  return retval;
else
  {
    printf("probcheck: at %d prob = %f\n",loc,prob);
    err = 1;
    early_exit(425);
    exit(1);
  }
}

double logck(arg,loc)
double arg;       /*check arg is in bounds (allow some error) */
int loc;
{
double retval;
int err=0;

retval = arg;

if(arg < 0.0) 
  {
    printf("logcheck: error at %d arg= %f < 0\n",loc,arg);
    err = 1;
  }

if(arg == 0.0)
  {
    printf("logcheck: warning at %d arg = 0\n",loc);
    retval = 1.0e-100;
  }
if(arg >= 0.0 && err == 0)
  return retval;
else
  {
    printf("logcheck: at %d arg = %f\n",loc,arg);
    early_exit(454);
    exit(1);
  }
}

double nonnegck(num,loc) /* check nonnegative */
double num;
{
if(num>=0.0)
  return(num);
else
    {
      printf("%g should be >= 0! at %d\n",num,loc);
      early_exit(466);
      exit(1);
    }
}

double numck(x, loc)
double x;
int loc;
{
  if(x>= 0.0||x< 0.0)
    return(x);
  else
    {
      printf("non-number at loc %d!\n", loc);
      return(x);
    }
}

double square(x)
double x;
{
return(x*x);
}

void prevcountsinit()
{
int pos1, char1;

for(pos1=0;pos1<LengthSeq;pos1++)
  {
    for(char1=0;char1<26;char1++)
      {
	prevcounts[pos1][char1] = 0.0;
      }
  }
}

void oldcountsinit()
{
int pos1, char1;

for(pos1=0;pos1<LengthSeq;pos1++)
  {
    for(char1=0;char1<26;char1++)
      {
	oldcounts[pos1][char1] = 0.0;
      }
  }
}


void selfchangeupdate()
{
double pisq;
int tpos, tchar;
double norm=0.0;

pisq=0.0;
for(tpos=0;tpos<LengthSeq;tpos++)
  {
    if(rateequil[tpos] > 1.0)
      {
	norm += rateequil[tpos];
	pisq+= rateequil[tpos]*oldprob[tpos][0]*oldprob[tpos][0];
	for(tchar=2;tchar<26;tchar++)  /*skip 1; it equals 1.0 (not an a.a.)*/
	  {
	      pisq+= rateequil[tpos]*oldprob[tpos][tchar]*oldprob[tpos][tchar];
	  }
      }
  }
Selfchange = (pisq+EPS)/(Rateequil*Eff_Length+EPS);
 if(Selfchange>=.999999){Selfchange=0.999999;}
}

void rateequilinit()
{
double maxp;
int tpos, tchar;

Rateequil=0.0;
Eff_Length = EPS;
for(tpos=0;tpos<LengthSeq;tpos++)
  {
    maxp=pck(oldprob[tpos][0],6601);
    for(tchar=2;tchar<26;tchar++)  /*skip 1; it equals 1.0 (not an a.a.)*/
      {
	  if(oldprob[tpos][tchar]>maxp) {
	    maxp = oldprob[tpos][tchar];
	  }
      }

    if(pck(maxp,6602) <EPS) 
      {
/*	printf("no positive probs at position %d; all x's?",tpos);
*/
	maxp = 1.0;
      }

    if(maxp < 1.0)
      {
	Eff_Length+=1.0;
	rateequil[tpos]=nonnegck(1.0/pck(maxp,6605),6705);
	Rateequil += rateequil[tpos];
      }
    else
      {
      	rateequil[tpos]=1.0;
      }
  }
Rateequil /= (Eff_Length + EPS);
selfchangeupdate();
}

void oldprobinit1()
{
int pos1, k;
char char1;
double norm;

for(pos1=0;pos1<LengthSeq;pos1++) /*initialize oldprob*/
  {
    norm = 20.0;
    for(char1='A';char1<='Z';char1++)
      {
	oldprob[pos1][char1-'A'] = 1.0;
      }
    oldprob[pos1]['J'-'A'] = 0.0;
    for(char1='A';char1<='Z';char1++)
      {
	for(k=0;k<NoSeq;k++)
	  {
	    if(seq[k][pos1] ==  char1)
	      {
		if(char1 != 'B' && char1 != 'Z' && char1 != 'X' &&
		   oldprob[pos1][char1-'A'] < 1.5)
		  {
		    norm += 1.0;
		    oldprob[pos1][char1-'A'] += 1.0;
		  }
		if(char1 == 'B' && oldprob[pos1]['N'-'A'] < 1.5)
		  {
		    norm += 1.0;
		    oldprob[pos1]['N'-'A'] += 1.0;
		  }
		if(char1 == 'B' && oldprob[pos1]['D'-'A'] < 1.5)
		  {
		    norm += 1.0;
		    oldprob[pos1]['D'-'A'] += 1.0;
		  }
		if(char1 == 'Z' && oldprob[pos1]['Q'-'A'] < 1.5)
		  {
		    norm += 1.0;
		    oldprob[pos1]['Q'-'A'] += 1.0;
		  }
		if(char1 == 'Z' && oldprob[pos1]['E'-'A'] < 1.5)
		  {
		    norm += 1.0;
		    oldprob[pos1]['E'-'A'] += 1.0;
		  }
	      }
	  }
      }

    for(char1='A';char1<='Z';char1++) 
      {
	oldprob[pos1][char1-'A'] /= norm;
      }
    oldprob[pos1]['B'-'A'] = 1.0; /* for use in renorm--uninheritance has no gain*/
    oldprob[pos1]['X'-'A'] = 0.0;
    oldprob[pos1]['Z'-'A'] = 0.0;
    oldprob[pos1]['O'-'A'] = 0.0;
    oldprob[pos1]['U'-'A'] = 0.0;
  }

rateequilinit();

}

void oldprobinit2()
{
int pos1, k;
char char1;
double norm;

for(pos1=0;pos1<LengthSeq;pos1++) /*initialize oldprob*/
  {
    norm = 0.0;
    for(char1='A';char1<='Z';char1++)
      {
	oldprob[pos1][char1-'A'] = 0.0;
      }
    oldprob[pos1]['J'-'A'] = 0.0;
    for(char1='A';char1<='Z';char1++)
      {
	for(k=0;k<NoSeq;k++)
	  {
	    if(seq[k][pos1] == char1)
	      {
		norm += 1.0;
		if(char1 != 'B' && char1 != 'Z'&& char1 != 'X')
		  oldprob[pos1][char1-'A'] += 1.0;
		if(char1 == 'B')
		  {
		    oldprob[pos1]['N'-'A'] += 0.5;
		    oldprob[pos1]['D'-'A'] += 0.5;
		  }
		if(char1 == 'Z')
		  {
		    oldprob[pos1]['Q'-'A'] += 0.5;
		    oldprob[pos1]['E'-'A'] += 0.5;
		  }
	      }
	  }
      }

    if(norm==0.0) printf("position with all x's? pos=%d\n",pos1);
    for(char1='A';char1<='Z';char1++) 
      {
	oldprob[pos1][char1-'A'] /= norm;
      }
    oldprob[pos1]['B'-'A'] = 1.0; /* for use in renorm--uninheritance has no gain*/
    oldprob[pos1]['X'-'A'] = 0.0;
    oldprob[pos1]['Z'-'A'] = 0.0;
    oldprob[pos1]['O'-'A'] = 0.0;
    oldprob[pos1]['U'-'A'] = 0.0;
  }

rateequilinit();

}

void oldprobiniteq()
{
int pos1, k;
char char1;
double norm;

for(pos1=0;pos1<LengthSeq;pos1++) /*initialize oldprob*/
  {
    norm = 10.0;
    for(char1='A';char1<='Z';char1++)
      {
	oldprob[pos1][char1-'A'] = 0.50;
      }

    for(char1='A';char1<='Z';char1++)
      {
	for(k=0;k<NoSeq;k++)
	  {
	    if(seq[k][pos1] ==  char1)
	      {
		if(char1 != 'B' && char1 != 'Z' && char1 != 'X')
		  {
		    norm += 1.0/oldprob[pos1][char1-'A'];
		    oldprob[pos1][char1-'A'] += 1.0/oldprob[pos1][char1-'A'];
		  }
		if(char1 == 'B')
		  {
		    norm += 0.5/oldprob[pos1]['N'-'A'] +
		      0.5/oldprob[pos1]['D'-'A'];
		    oldprob[pos1]['N'-'A'] += 0.5/oldprob[pos1]['N'-'A'];
		    oldprob[pos1]['D'-'A'] += 0.5/oldprob[pos1]['D'-'A'];
		    
		  }
		if(char1 == 'Z')
		  {
		    norm += 0.5/oldprob[pos1]['Q'-'A'] +
		      0.5/oldprob[pos1]['E'-'A'];
		    oldprob[pos1]['Q'-'A'] += 0.5/oldprob[pos1]['Q'-'A'];
		    oldprob[pos1]['E'-'A'] += 0.5/oldprob[pos1]['E'-'A'];
		  }
	      }
	  }
      }

    for(char1='A';char1<='Z';char1++) 
      {
	oldprob[pos1][char1-'A'] /= norm;
      }
    oldprob[pos1]['B'-'A'] = 1.0; /* for use in renorm--uninheritance has no gain*/
    oldprob[pos1]['X'-'A'] = 0.0;
    oldprob[pos1]['Z'-'A'] = 0.0;
    oldprob[pos1]['O'-'A'] = 0.0;
    oldprob[pos1]['U'-'A'] = 0.0;
  }

rateequilinit();

}




double dlike(p,seq1,seq2)  /* with thanks to aaron halpern */
double p;
int seq1,seq2;
{
int tpos,tchar,tchar2;
double sumval=0.0,tprob, tprobm,gamma, pi;

for(tpos=0;tpos<LengthSeq;tpos++)
{
if(seq[seq1][tpos]!='X' && seq[seq2][tpos]!='X')
  {
    tchar = seq[seq1][tpos];
    tchar2 = seq[seq2][tpos];
    gamma = rateequil[tpos];          /* rate of equilibration= # of codons*/
    if(tchar == 'B')                 /* pi = prob of char at saturation */
      pi = oldprob[tpos]['N' - 'A'] + oldprob[tpos]['D' - 'A']; 
    else if(tchar == 'Z')
      pi = oldprob[tpos]['Q' - 'A'] + oldprob[tpos]['E' - 'A']; 
    else 
      pi = oldprob[tpos][tchar - 'A'];
    tprobm = pow(p,gamma-1.);             /* p^(gamma-1) */
    if(tprobm >= 0.0 && tprobm <= 1.0)
      {
	tprobm = tprobm;
      }
    else
      {
	tprobm = 1.0;
      }
    tprob = tprobm*p;
    if(tchar==seq[seq2][tpos])
      {
	sumval += (1.-pi)*gamma*tprobm/(tprob+(1-tprob)*pi);
      }
    else if((tchar == 'B' && (tchar2 == 'D' || tchar2 == 'N')) ||
	    (tchar == 'Z' && (tchar2 == 'E' || tchar2 == 'Q')) )
      sumval += (1.-pi)*gamma*tprobm/(tprob+(1-tprob)*pi);
    else if((tchar2 == 'B' && (tchar == 'D' || tchar == 'N')) ||
	    (tchar2 == 'Z' && (tchar == 'E' || tchar == 'Q')) )
      {
	if(tchar2 == 'B')                 /* pi = prob of char at saturation */
	  pi = oldprob[tpos]['N' - 'A'] + oldprob[tpos]['D' - 'A']; 
	if(tchar2 == 'Z')
	  pi = oldprob[tpos]['Q' - 'A'] + oldprob[tpos]['E' - 'A']; 
	sumval += (1.-pi)*gamma*tprobm/(tprob+(1-tprob)*pi);
      }
    else	    /* sequences disagree */
      {
	sumval -= gamma*tprobm/(1.-tprob);
      }
  }
}
return(sumval + 1/p); /*multiplying by p before differentiating the log
		      is like using an exponential prior on d */
}

double d2like(p,seq1,seq2)  /* second derivaative */
double p;
int seq1,seq2;
{
int tpos,tchar, tchar2;
double sumval=0.0,tprob, tprobm,tprobm2,gamma, pi;

for(tpos=0;tpos<LengthSeq;tpos++)
{
if(seq[seq1][tpos]!='X' && seq[seq2][tpos]!='X')
  {
    tchar = seq[seq1][tpos];
    tchar2 = seq[seq2][tpos];
    gamma = rateequil[tpos];          /* rate of equilibration= # of codons*/
    if(tchar == 'B' || tchar2 == 'B')    /* pi = prob of char at saturation */
      pi = oldprob[tpos]['N' - 'A'] + oldprob[tpos]['D' - 'A']; 
    else if(tchar == 'Z' || tchar2 =='Z')
      pi = oldprob[tpos]['Q' - 'A'] + oldprob[tpos]['E' - 'A']; 
    else 
      pi = oldprob[tpos][tchar - 'A'];
    tprobm2 = pow(p,gamma-2.0);             /* p^(gamma-1) */
    tprobm = tprobm2*p;
    tprob = tprobm*p;
    if(tchar==tchar2 || 
       ((tchar == 'B' && (tchar2 == 'D' || tchar2 == 'N')) ||
	(tchar == 'Z' && (tchar2 == 'E' || tchar2 == 'Q')) ) ||
       ((tchar2 == 'B' && (tchar == 'D' || tchar == 'N')) ||
	    (tchar2 == 'Z' && (tchar == 'E' || tchar == 'Q')) ))  /*agree */
      {
	sumval += (1.-pi)*gamma*
	  ((tprob*(1.-pi)+pi)*(gamma-1.)*tprobm2-tprobm*gamma*tprobm*(1.-pi))
	    /((tprob+(1-tprob)*pi)*(tprob+(1-tprob)*pi));
      }
    else  /* disagree */
      {
	sumval -= gamma*
	  ((1.-tprob)*(gamma-1.)*tprobm2+tprobm*gamma*tprobm)
	    /((1.-tprob)*(1.-tprob));
      }
  }
}
return(sumval);
}

double pmoodlike(p,seq1,seq2)
double p;
int seq1, seq2;
{
return(-1.0/dlike(p,seq1,seq2)-p);
}

double bbisect(x1,x2,tol,seq1,seq2,func) /* find root of func using bisection*/
double x1,x2,tol;                      /* func must be decreasing */
int seq1,seq2;
double (*func)();
{
int i=0;
double cx1,cx2,f1,f2,fnew,xnew;

cx1= x1;
cx2= x2;

f1= func(cx1,seq1,seq2);
f2= func(cx2,seq1,seq2);
if(f1 == f2) return((cx1+cx2)/2.);
if(f1 <= 0.0) return(cx1);
if(f2 >= 0.0) return(cx2);

 while(i<ITMAX && cx2-cx1 > TOL){
   xnew = (cx1+cx2)*0.5;
   fnew = func(xnew,seq1,seq2);
   if(fnew==0.0) return(xnew);
   else if(fnew > 0.0) cx1 = xnew;
   else cx2= xnew;
   i++;
 }
if(i==ITMAX) printf("too many its?: %d",i);

return(0.5*(cx1+cx2));
}

double findzero(tol,seq1,seq2,weightptr)
double tol, *weightptr;
int seq1,seq2;
{
  double p=EPS, tlike, cross=0.0, ncross= EPS;

  tlike = dlike(p, seq1, seq2);

  while((p < .98) && tlike <= 0.0 && cross == 0.0) 
    {
      p+= .02;
      tlike = dlike(p, seq1, seq2);
      if( -1.0 > tlike*p) 
	{
	  cross= p;
	}
      else ncross = p;
    }
  if(cross == 0.0)
    {
      p = bbisect(p,1.0-.5*tol,tol,seq1,seq2,dlike);
/*      printf("max like = %f\n",p);
*/
/*      *weightptr = -d2like(p,seq1,seq2); removed for speed */
      return(p);
    }
  else
    { /* infinite m.l. distance; must make finite estimate */
/*      printf("pmood(ncross) = %lg pmood(cross)= %lg cross=%f\n",
	     pmoodlike(ncross,seq1,seq2),pmoodlike(cross,seq1,seq2),cross);
*/
      p = bbisect(ncross,cross,tol,seq1,seq2,pmoodlike);
/*      printf("max like = 0; using %lg\n",0.5*p); */
      npzeros++;
      *weightptr = 4./(p*p);
/*      return(0.5*p); */ /* try something more "rigorous" */
/* if we integrate sqrt(2/pi*s^2) -log(p) e^(-p^2/2s^2) dp from 0 to infinity,
   (s is sigma, = root of pmoodlike above) and take exp(- result) we get
   s e^(-gamma/2)/sqrt(2) = .5298393548 s */
      return(0.5298393548*p);
    }
}

void early_exit(loc)
int loc;
{
printf("program NOT converged, exiting prematurely from line %d!\n",loc);
fflush(stdout);
fflush(stderr);
exit(1);
/* savecounts();
calcprobs();
printprobs();
printmodel();
calcentropy();
printent();
if(system("mv besttree oldbesttree")!= 0)
    printf("could not mv besttree!\n");
printbesttree();
exit(1);
*/
}


int writedistances() /* calculate and write distances*/
{
int seq1, seq2, j;
double p,d,weight;
FILE *fp;
/* FILE *fpe; 
int k;
*/

zerodist = 0;

if((fp=Fopen("infile","w")) == NULL) 
  {
    printf("can't open infile for writing");
    return(1);
  }

/* weightfile removed for speed
if((fpe=Fopen("weightfile","w")) == NULL) 
  {
    printf("can't open weightfile for writing");
    return(1);
  }
*/

fprintf(fp," %d\n",NoSeq);

npzeros = 0;
for(seq1=0;seq1<NoSeq;seq1++)
  {

    for(j=0;j<NAMELENGTH;j++)
	fprintf(fp,"%c",names[seq1][j]);
    fprintf(fp," ");

    for(seq2=0;seq2<NoSeq;seq2++)
      {
	if(seq2>=seq1)
	  {
	    p=pck(findzero(.0000001,seq1,seq2,&weight),5);
	    d=log(logck(1.0/p,5));
	    dist[seq1][seq2] = d;
	    if(d==0.0 && seq2>seq1) zerodist = 1;
	  }
	else
	  d= dist[seq2][seq1];
	fprintf(fp,"%f ",d*Rateequil*(1.0-Selfchange));
/*	fprintf(fpe,"%f ",weight*p*p);  removed for speed*/
      }
    fprintf(fp,"\n");
  }
if(npzeros > 0) printf("%d prs w/ inf. dist,",npzeros);
fclose(fp);
/*fclose(fpe);*/
 return(0);
}
/* -------- ReadTree ------------------------
 *
 */
Treeptr allocatetree()
{
Treenode *top;

if(NoSeq > 3) 
  {
    if((top = (Treenode *) malloc((2*NoSeq-2)*sizeof(Treenode)))==NULL)
      {
	printf("allocatetree: can't allocate tree memory for tree %d\n",
	       NoTrees+1);
	printf("will use %d trees\n",NoTrees);
      }
  }
else
  {
    if((top = (Treenode *) malloc((2*3-2)*sizeof(Treenode)))==NULL)
      {
	printf("allocatetree: can't allocate tree memory for tree %d\n",
	       NoTrees+1);
	printf("will use %d trees\n",NoTrees);
      }
  }

return(top);
}

int ReadTree( top, in_fp )
     Treeptr top;
     FILE *in_fp;
{    	
  Treenode *here;
  char a, punct[2],name[NAMELENGTH+1];
  int seqno, count, memcount, numneg=0, j, tseqno, found[MaxNoSeq];
  double worstneg=0.0;

  for(j=0;j<NoSeq;j++)
    found[j]=0;
  count = -1;
  memcount=0;
  here = &(top[memcount]);
  (*here).seqno = count;
  memcount ++;
  (*here).side = 't';
  count --;

  while ( top != NULL && (fscanf(in_fp,"%c",&a) == 1 )) {

    switch (a) {
    case '(':
      (*here).left = &(top[memcount]);
      memcount ++;
      (*(*here).left).side = 'l';
      (*(*here).left).mother = here;
      (*(*here).left).left = (Treenode *) NULL;
      (*(*here).left).right = (Treenode *) NULL;
      (*(*here).left).seqno = count;
      count--;
      (*here).right = &(top[memcount]);
      memcount ++;
      (*(*here).right).side = 'r';
      (*(*here).right).mother = here;
      (*(*here).right).left = (Treenode *) NULL;
      (*(*here).right).right = (Treenode *) NULL;
      (*(*here).right).seqno = count;
      count--;

      here = (*here).left;
      break;
    case ')':
      here = (*here).mother;
      break;
    case ',':
      if( (*(*here).mother).left == here)
	here = (*(*here).mother).right;
      else if( (*(*here).mother).right == here)
	{
	  (*(*here).mother).mother = &(top[memcount]);
	  memcount ++;
	  (*(*(*here).mother).mother).side = 'm';
	  (*(*here).mother).side = 'm';
	  (*(*(*here).mother).mother).mother = (*here).mother;
	  (*(*(*here).mother).mother).left = (Treenode *) NULL;
	  (*(*(*here).mother).mother).right = (Treenode *) NULL;
	  (*(*(*here).mother).mother).seqno = count;
	  count--;

	  here = (*(*here).mother).mother;
	}
      else
	{ top = NULL; printf("four way branch not allowed!\n");}
      break;
    case ':':
      if (fscanf(in_fp,"%lf", &((*here).length)) != 1) top = NULL;
      (*here).length /= ((Rateequil+EPS)*(1.0+EPS-Selfchange));
      if((*here).length < MINDIST) 
	{
	  if((*here).length < worstneg) worstneg = (*here).length;
	  (*here).length = MINDIST;
	  numneg ++;
	}
      break;
    case '\n':
      break;
    case ';':
      break;
    default:
      name[0] = a; /* read rest of name: */
      j=1;		  
      while ( fscanf(in_fp,"%c",&a) == 1  && 
	     a != ':' && a !=',') 
	{
	  if(j<NAMELENGTH) 
	    name[j] = a;
	  else
	    {
	      name[j] = '\0';
	      printf("name in treefile too long: %s%c\n",name,a);
	      early_exit(1222);
	    }
	  j++;
	}
      name[j] = '\0';
      punct[0]= a;
      punct[1] = '\0';


      seqno=-1;
      for(tseqno=0;tseqno< NoSeq;tseqno++)
	{
	  if(comparenames(name,names[tseqno])==1)
	    {
	      if(seqno != -1)
		{
		  printf("name matches 2 names!\n");
		  printf("seq %d:",tseqno);
		  for(j=0;j<NAMELENGTH;j++)
		    printf("%c",names[tseqno][j]);
		  printf("\n");
		  printf("seq %d:",seqno);
		  for(j=0;j<NAMELENGTH;j++)
		    printf("%c",names[seqno][j]);
		  printf("\n");
		  fprintf(stderr,"exiting...\n");
		  early_exit(1197);
		}
	      if(found[tseqno]!=0)
		{
		  printf("same sequence in tree twice!\n");
		  printf("seq %d:",tseqno);
		  for(j=0;j<NAMELENGTH;j++)
		    printf("%c",names[tseqno][j]);
		  printf("\n");
		  fprintf(stderr,"exiting...\n");
		  early_exit(1207);
		}
	      seqno = tseqno;
	      found[seqno]=1;
	    }
	}
      if(seqno<0)
	{
	  printf("did not match sequence\n");
		  for(j=0;j<NAMELENGTH;j++)
		    printf("%c",name[j]);
		  printf("\n");
	  early_exit(1216);
	}
      (*here).seqno = seqno;

/*    count ++; doesn't work--play it safe*/

      treeloc[Treeit][seqno] = here;
      if(seqno>= NoSeq)
	(*treeloc[Treeit][seqno]).inc = 'n';
      if (punct[0] == ',') {
	here = (*(*here).mother).right;
	break;
      }
      else {
	if (punct[0] == ':') {
	  if (fscanf(in_fp,"%lf", &((*here).length)) != 1) top = NULL;
	  nonnegck((*here).length,1261);
	  (*here).length /= ((Rateequil+EPS)*(1.0+EPS-Selfchange));
	  if((*here).length < MINDIST) 
	    {
	      if((*here).length < worstneg) worstneg = (*here).length;
	      (*here).length = MINDIST;
	      numneg ++;
	    }
	  nonnegck((*here).length,1267);
	  break;
	}
	else {
	  printf("after name=%s:\n",name);
	  printf("problem with punctuation--punct[0]=%c:\n",punct[0]);
	  printf("treefile ends prematurely?");
	  top = NULL;
	  break;
	}
      }
    }
  }
  if(NoSeq >= 3 && memcount != 2*NoSeq-2)
    {
      printf("did not fill tree!");
      printf("seqno = %d, memcount = %d",seqno,memcount);
      top=NULL;
    }
  if(NoSeq < 3 && memcount != 4)
    {
      printf("did not fill tree!");
      printf("seqno = %d, memcount = %d",seqno,memcount);
      top=NULL;
    }
  if (top != here) {
    printf("error: tree did not parse\n");
    if(top == NULL) printf("top == Null\n");
    printf("top = %d, here = %d\n",(*top).seqno,(*here).seqno);
    top = NULL;
  }
  else  if((*here).side != 'm')
    {
      printf("error: tree did not parse--should be unrooted\n");
    }
  (*here).length = (*(*here).mother).length;

  if(numneg > 0) printf(" %d lengths<%g, worst: %6g,",numneg,
			MINDIST,worstneg);

  for(j=0;j<NoSeq;j++)
    {
      if(found[j]!=1) fprintf(stderr,"didn't find taxon %d\n!",j);
    }

  if (top == NULL) return(1);
  else return(0);
}


int imin(a,b)
int a,b;
{
  if(a<=b) return(a);
  else return(b);
}

int vrenorm(length1,length2,ndata1,ndata2,ndata,calc)
double length1, length2; 
Renormdata *ndata1, *ndata2, *ndata;
char calc;
{
double p1, p2, omp1, omp2,   sum,tmp,tmpa,tmpb,gamma;
int char1,i,j,retval=0;
double *pcharvect1, *pcharvect2, *pcharvectret, *pnvect1, *pnvect2, *pnvectret;
int maxn1, maxn2, maxn;
double *virtcts1, *virtcts2, *virtctsret;
double *virtctsp1, *virtctsp2, *virtctspret;
double negvirtct1, negvirtct2, negvirtctp1, negvirtctp2;

pcharvect1 = (*ndata1).pcharvect;
pcharvect2 = (*ndata2).pcharvect;
pcharvectret = (*ndata).pcharvect;

virtcts1 = (*ndata1).virtcts;
virtcts2 = (*ndata2).virtcts;
virtctsret = (*ndata).virtcts;

virtctsp1 = (*ndata1).virtctsp;
virtctsp2 = (*ndata2).virtctsp;
virtctspret = (*ndata).virtctsp;


maxn1 = (*ndata1).maxn[Char-'A'];
maxn2 = (*ndata2).maxn[Char-'A'];
maxn = maxn1 + maxn2;

if((*ndata).maxn['Z'-'A'] +maxn> 2*NoSeq) 
{
  printf("vrenorm: index out of bounds\n");
  printf("maxn=%d, NoSeq=%d",(*ndata).maxn['Z'-'A'],NoSeq);
  early_exit(1319);
}
(*ndata).maxn[Char-'A'] = maxn;
if(maxn1>0)
     pnvect1 = ((*ndata1).pnvect)[Char-'A'];
if(maxn2>0)
     pnvect2 = ((*ndata2).pnvect)[Char-'A'];
if(maxn>0)
pnvectret = (*ndata).data + (*ndata).maxn['Z'-'A'];
if(maxn>0)
((*ndata).pnvect)[Char-'A']=pnvectret;

(*ndata).maxn['Z'-'A'] += maxn;

gamma = rateequil[Pos];          /* rate of equilibration*/

 if(DEBUG){
   printf("gam %lg l %lg\n",gamma,length1);
 }

p1 = pck(exp(-gamma*length1),-1);
p2 = pck(exp(-gamma*length2),-2);
omp1 = -expm1(-gamma*length1); /* 1-p1 */
omp2 = -expm1(-gamma*length2); /* 1-p2 */
    if(omp1 > 0)
       {
	 negvirtct1 = -gamma*length1*p1/omp1;
	 negvirtctp1 = (gamma*length1*p1*(omp1-gamma*length1)/omp1)/omp1;
       }
    else
       {
	 negvirtct1 = - 1.0;
	 negvirtctp1 = - 0.0;
       }

    if(omp2 > 0)
       {
	 negvirtct2 = -gamma*length2*p2/omp2;
	 negvirtctp2 = (gamma*length2*p2*(omp2-gamma*length2)/omp2)/omp2;
       }
    else
       {
	 negvirtct2 = - 1.0;
	 negvirtctp2 = - 0.0;
       }


/* these have been normalized to 1.0 so we will get rid of them :
sum1 = 1.0;
sum2 = 1.0;
*/

if(calc == 'n') /*not yet calculated pcharvect*/
{
  for(char1=0;char1<25;char1++)
  {
    tmpa = 0.0;
    tmpb = 0.0;
    if(char1==1) /* this is where we keep uninherited characters*/
      {
	tmpa = (omp1)*(omp2);
	tmpb =  -2.0*p1*p2*pcharvect1[1]*pcharvect2[1]; 
      }
    if(pcharvect1[char1] > 0.0 || pcharvect2[char1] > 0.0 || tmpa>0.0)
      {
	pcharvectret[char1] =  nonnegck(tmpa + tmpb +
	p1*pcharvect1[char1]*p2*pcharvect2[char1]/oldprob[Pos][char1] +
	  p1*pcharvect1[char1]*(omp2+p2*pcharvect2[1]) + 
	    p2*pcharvect2[char1]*(omp1+p1*pcharvect1[1]),101) ;
	virtctsret[char1] = 
	  tmpa*
	    (virtcts1[25]+virtcts2[25]+negvirtct1+negvirtct2)
	      +
	  tmpb*
	    (virtcts1[1]+virtcts2[1]+gamma*(length1+length2))               +
	  p1*pcharvect1[char1]*p2*pcharvect2[char1]/oldprob[Pos][char1]*
	    (virtcts1[char1]+virtcts2[char1]+gamma*(length1+length2))       +
          p1*pcharvect1[char1]*(omp2*
            (virtcts1[char1]+virtcts2[25]+gamma*length1+negvirtct2)+
				p2*pcharvect2[1]*
	    (virtcts1[char1]+virtcts2[1] +gamma*(length1+length2)) )        +
          p2*pcharvect2[char1]*(omp1*
            (virtcts2[char1]+virtcts1[25]+gamma*length2+negvirtct1)+
				p1*pcharvect1[1]*
	    (virtcts2[char1]+virtcts1[1] +gamma*(length1+length2)) )        ;  

	virtctspret[char1] = 
  tmpa*
    (square(virtcts1[25]+virtcts2[25]+negvirtct1+negvirtct2)
     +
     (virtctsp1[25]+virtctsp2[25]+negvirtctp1+negvirtctp2))
	      +
  tmpb*
    (square(virtcts1[1]+virtcts2[1]+gamma*(length1+length2)) +
     (virtctsp1[1]+virtctsp2[1]-gamma*(length1+length2)))
      +
  p1*pcharvect1[char1]*p2*pcharvect2[char1]/oldprob[Pos][char1]*
    (square(virtcts1[char1]+virtcts2[char1]+gamma*(length1+length2)) +
     (virtctsp1[char1]+virtctsp2[char1]-gamma*(length1+length2)))
     +
  p1*pcharvect1[char1]*omp2*
    (square(virtcts1[char1]+virtcts2[25]+gamma*length1+negvirtct2) +
     (virtctsp1[char1]+virtctsp2[25]-gamma*length1+negvirtctp2))
     +
  p1*pcharvect1[char1]*p2*pcharvect2[1]*
    (square(virtcts1[char1]+virtcts2[1] +gamma*(length1+length2)) +
     (virtctsp1[char1]+virtctsp2[1]-gamma*(length1+length2)))
     +
  p2*pcharvect2[char1]*omp1*
    (square(virtcts2[char1]+virtcts1[25]+gamma*length2+negvirtct1)+
     (virtctsp2[char1]+virtctsp1[25]-gamma*length2+negvirtctp1))
     +
  p2*pcharvect2[char1]*p1*pcharvect1[1]*
    (square(virtcts2[char1]+virtcts1[1] +gamma*(length1+length2))+
     (virtctsp2[char1]+virtctsp1[1]-gamma*(length1+length2)));  
      }
    else
      {
	pcharvectret[char1] = 0.0;
	virtctsret[char1] = 0.0;
	virtctspret[char1] = 0.0;
      }
  }

  sum = 0.0;
  virtctsret[25] = 0.0;
  virtctspret[25] = 0.0;
  for(char1=0;char1<25;char1++) 
    {
      sum += pcharvectret[char1];
      virtctsret[25] += virtctsret[char1];
      virtctspret[25] += virtctspret[char1];
    }
  if(sum ==0.0) 
    {
      retval=-1;
      printf("vrenorm: impossible tree!\n");
    }
  for(char1=0;char1<25;char1++) 
    {
      if(pcharvectret[char1] > 0.0)
	{
	  virtctsret[char1] /= pcharvectret[char1];
	  virtctspret[char1] /= pcharvectret[char1];
	  virtctspret[char1] -= square(virtctsret[char1]);

	}
      else
	{
	  virtctsret[char1] = 0.0;
	  virtctspret[char1] = 0.0;
	}
      pcharvectret[char1] /= sum;
    }
  pcharvectret[25] = log(logck(sum,251))+pcharvect1[25]+pcharvect2[25];
  virtctsret[25] /= sum;
  virtctspret[25] /= sum;
  virtctspret[25] -= square(virtctsret[25]);

  if(pcharvectret[25] > 0.000001) printf("warning: likelihood > 1 (log=%g)\n",
				    pcharvectret[25]);
}


for(i=0;i<maxn;i++)
{
  tmp = 0.0;
  for(j=i-imin(maxn1,i);j<imin(maxn2,i);j++)
    {
      if(i<0 || (i-j)-1 < 0) printf("indexing screwed up!");
      tmp += p1*pnvect1[(i-j)-1]*p2*pnvect2[j]/oldprob[Pos][Char-'A'];
    }
  pnvectret[i] = tmp;
}

for(i=0;i<maxn1;i++)
     pnvectret[i] += p1*pnvect1[i]*(omp2+p2*pcharvect2[1]);

for(i=0;i<maxn2;i++)
     pnvectret[i] += p2*pnvect2[i]*(omp1+p1*pcharvect1[1]);

/* some debugging below; still have to normalize pnvectret*/
tmp = 0.0;

for(i=0;i<maxn;i++)
     tmp += pnvectret[i];
     
if(tmp < 0.0) printf("warning: negative tmp %g",tmp);
if(maxn > 0 && tmp > 0.0)
{
  tmp = pcharvectret[Char-'A']/tmp;     
  for(i=0;i<maxn;i++)
    {
	pnvectret[i] = pck(pnvectret[i]*tmp,881);
    }
}
return(retval);
}

int srenorm(length1, ndata1,ndata,calc)
double length1;
Renormdata *ndata1, *ndata;
char calc;
{
int maxn1, maxn;
double p1, omp1, sum1,gamma, tmp;
char char1;
int i, retval=0;
double *pcharvect1, *pcharvectret, *pnvect1, *pnvectret;
double *virtcts1, *virtctsret;
double *virtctsp1, *virtctspret;
double negvirtct1, negvirtctp1;

pcharvect1 = (*ndata1).pcharvect;
pcharvectret = (*ndata).pcharvect;

virtcts1 = (*ndata1).virtcts;
virtctsret = (*ndata).virtcts;

virtctsp1 = (*ndata1).virtctsp;
virtctspret = (*ndata).virtctsp;

pnvectret = (*ndata).data + (*ndata).maxn['Z'-'A'];
(*ndata).pnvect[Char-'A']=pnvectret;

maxn1 = (*ndata1).maxn[Char-'A'];
maxn = maxn1;
(*ndata).maxn['Z'-'A'] += maxn;
if((*ndata).maxn['Z'-'A'] > 2*NoSeq) 
{
  printf("srenorm: index out of bounds\n");
  early_exit(1544);
}
(*ndata).maxn[Char-'A'] = maxn;

if(maxn1>0 && !(((*ndata1).pnvect[Char-'A'] < (*ndata1).data+2*NoSeq) && 
     ((*ndata1).pnvect[Char-'A'] >= (*ndata1).data)))
  {
    fprintf(stderr,"srenorm error:\n");
    fprintf(stderr,"problem passing pointer, or pnvect out of bounds!\n");
    printf("Char = %c\n",Char);
    early_exit(1554);
  }

if(maxn1>0)
     pnvect1 = (*ndata1).pnvect[Char-'A'];

gamma = rateequil[Pos];          /* rate of equilibration= # of codons*/

p1 = pck(exp(-gamma*length1),-1);
omp1 = -expm1(-gamma*length1); /* 1-p1 */

 if(omp1 > 0.0) /* previously (gamma*length1 > 0) */
       {
	 negvirtct1 = -gamma*length1*p1/omp1;
	 negvirtctp1 = (gamma*length1*p1*(omp1-gamma*length1)/omp1)/omp1;
       }
    else
       {
	 negvirtct1 = - 1.0;
	 negvirtctp1 = - 0.0;
       }

sum1 = 1.0;

if(calc == 'n') /*not yet calculated pcharvect*/
{
  for(char1='A';char1<='Z';char1++)
  {
    pcharvectret[char1-'A'] = p1*pcharvect1[char1-'A'];
    virtctsret[char1-'A'] = 
                          (virtcts1[char1-'A']+gamma*length1);
    virtctspret[char1-'A'] = 
			   (virtctsp1[char1-'A']-gamma*length1);/*squares cancel*/
  }
  pcharvectret[1] += (omp1)*sum1;  /* this is where we keep uninherited ch*/
  virtctsret[1] *= p1*pcharvect1[char1-'A'];
  virtctsret[1] += omp1*sum1*(virtcts1[25]+negvirtct1);
  if(pcharvect1[1] > 0.0)
    virtctsret[1] /= pcharvect1[char1-'A'];
  virtctspret[1] *= p1*pcharvect1[char1-'A'];
  virtctspret[1] += omp1*sum1*
              (square(virtcts1[25]+negvirtct1) +
	       (virtctsp1[25]+negvirtctp1));
  if(pcharvect1[1] > 0.0)
    virtctspret[1] /= pcharvectret[1];
  virtctspret[1] -= square(virtctsret[1]);
  pcharvectret[25] = pcharvect1[25];
  virtctsret[25] = virtcts1[25];
  virtctspret[25] = virtctsp1[25];
}


for(i=0;i<maxn1;i++) 
    pnvectret[i] = p1*pnvect1[i];

/* just make sure it worked; remove the following later:*/
tmp = 0.0;
for(i=0;i<maxn;i++)
     tmp += pnvectret[i];
if(fabs(tmp-pcharvectret[Char-'A'])> .00001) 
{
  printf("error in srenorm\n");
  printf("tmp= %g pcharvret = %g\n",tmp, pcharvectret[Char-'A']);
  printf("Pos= %d, Char= %c, maxn= %d\n",Pos,Char,maxn);
  early_exit(1618);
}

return(retval);

}


void updateleaf(nodeptr)
Treeptr nodeptr;
{
int i, seqno, *maxn;
double *pcharvectret, sum=0.0;
Renormdata *ndata;

ndata = &((*nodeptr).mrdata);
seqno=(*nodeptr).seqno;
pcharvectret = (*ndata).pcharvect;
maxn = (*ndata).maxn;

if((*nodeptr).mcalc['Z'-'A']!='y')
{
  (*nodeptr).mcalc['Z'-'A'] = 'y';
  for(i=0;i<25;i++)
    {
      pcharvectret[i] = 0.0;
      (*ndata).virtcts[i] = 0.0;
      (*ndata).virtctsp[i] = 0.0;
    }
  (*ndata).virtcts[25] =0.0;
  (*ndata).virtctsp[25] =0.0;

  if((*nodeptr).inc == 'y')
   {

    if(seq[seqno][Pos] != 'B' && seq[seqno][Pos] != 'Z')
      {
	if(oldprob[Pos][seq[seqno][Pos]-'A']==0.0) 
	  printf("warning: zero prob (char = %c)\n",seq[seqno][Pos]);
	pcharvectret[25] = 
	  log(logck(pck(oldprob[Pos][seq[seqno][Pos]-'A'],25),25)); 
	pcharvectret[seq[seqno][Pos]-'A'] = 1.0;
      }
    if(seq[seqno][Pos] == 'B')
      {
	sum = pck(oldprob[Pos]['N'-'A'] + oldprob[Pos]['D'-'A'],26); 
	pcharvectret['N'-'A'] = 
	  oldprob[Pos]['N'-'A']/sum; 
	pcharvectret['D'-'A'] = 
	  oldprob[Pos]['D'-'A']/sum; 
	pcharvectret[25] = log(logck(sum,265));
      }
    if(seq[seqno][Pos] == 'Z')
      {
	sum = pck(oldprob[Pos]['Q'-'A'] + oldprob[Pos]['E'-'A'],27); 
	pcharvectret['Q'-'A'] = 
	  oldprob[Pos]['Q'-'A']/sum; 
	pcharvectret['E'-'A'] = 
	  oldprob[Pos]['E'-'A']/sum; 
	pcharvectret[25] = log(logck(sum,275));
      }
   }
  else  /* 'x' character or some such uninformative case */
   {
     pcharvectret[1] = 1.0;
     pcharvectret[25] = 0.0;
   }
}

if((*nodeptr).mcalc[Char-'A']!='y'&& (*nodeptr).inc=='y')
  {
    (*nodeptr).mcalc[Char-'A']='y';
    if(pcharvectret[Char-'A'] > 0.0)
      {
	(*ndata).pnvect[Char-'A'] = (*ndata).data + maxn['Z'-'A'];
	maxn['Z'-'A'] += 1;
	maxn[Char-'A'] = 1;
	(*ndata).pnvect[Char-'A'][0]=pcharvectret[Char-'A']; 
      }
    else maxn[Char-'A'] = 0;
    
  }
if((*nodeptr).inc == 'y' && pcharvectret[Char-'A'] > 0.0 &&
(!(((*ndata).pnvect[Char-'A'] < (*ndata).data+2*NoSeq) && ((*ndata).pnvect[Char-'A'] >= (*ndata).data)))
 )
  {

    fprintf(stderr,
	    "updateleaf: problem passing pointer, or pnvect out of bounds!\n");
    /*    fprintf(stderr,
	    "Pos= %d, Char = %c, maxn= %d, pnvect[char]= %g data = %p\n",
	    Pos, Char, maxn['Z'-'A'], *(ndata->pnvect+(Char-'A')),
	    ndata->data); */
    fprintf(stderr,
	    "maxn[char]= %d, pcharvectret[char] = %g, seqno=%d\n",
	   maxn[Char-'A'],pcharvectret[Char-'A'],seqno);
    early_exit(1714);
  }

}

Renormdata *finddata(nodeptr,dir)
Treeptr nodeptr;
char dir;
{
switch (dir) 
  {
  case 'm':
    return(&((*nodeptr).mrdata));
    break;
  case 'l':
    return(&((*nodeptr).lrdata));
    break;
  case 'r':
    return(&((*nodeptr).rrdata));
    break;
  }
 printf("unknown direction!\n");
 exit(1);
}

int nodeupdate(nodeptr,dir)
Treeptr nodeptr;
char dir;
{
double  length1, length2;
 int  seqno, retval;
Treeptr nodeptr1, nodeptr2;
char dir1, dir2, calc;
Renormdata *ndata, *ndata1, *ndata2;
double p1,p2,gamma;

retval = 0;

if(nodeptr == NULL)
  {
    printf("at empty node!");
    return(1);
  }

/*printf("at node: %d, dir= %c inc= %c \n",(*nodeptr).seqno,dir,(*nodeptr).inc, *maxnptr); */

seqno = (*nodeptr).seqno;

if(seqno >= 0)  
  {                          /* at a leaf */
    if(dir!='m') 
      {
	printf("error:leaves have no children!");
	return(-1);
      }
    updateleaf(nodeptr);

  }
else
  {
    switch (dir) 
      {
      case 'm':
	if((*nodeptr).mcalc[Char-'A'] == 'y')
	  {
	    return(retval);
	  }
	else
	  {
	    (*nodeptr).mcalc[Char-'A'] = 'y';
	    calc = (*nodeptr).mcalc[25];
	    (*nodeptr).mcalc[25] = 'y';
	    nodeptr1 = (*nodeptr).left;
	    dir1 = 'm';
	    nodeptr2 = (*nodeptr).right;
	    dir2 = 'm';
	  }
	break;
      case 'l':
	if((*nodeptr).lcalc[Char-'A'] == 'y')
	  {
	    return(retval);
	  }
	else
	  {
	    (*nodeptr).lcalc[Char-'A'] = 'y';
	    calc = (*nodeptr).lcalc[25];
	    (*nodeptr).lcalc[25] = 'y';
	    nodeptr1 = (*nodeptr).mother;
	    dir1 = (*nodeptr).side;
	    nodeptr2 = (*nodeptr).right;
	    dir2 = 'm';
	  }
	break;
      case 'r':
	if((*nodeptr).rcalc[Char-'A'] == 'y')
	  {
	    return(retval);
	  }
	else
	  {
	    (*nodeptr).rcalc[Char-'A'] = 'y';
	    calc = (*nodeptr).rcalc[25];
	    (*nodeptr).rcalc[25] = 'y';
	    nodeptr1 = (*nodeptr).mother;
	    dir1 = (*nodeptr).side;
	    nodeptr2 = (*nodeptr).left;
	    dir2 = 'm';
	  }
	break;
      default:
	printf("nrenorm: unknown side-type");
	break;
      }

    if((ndata = finddata(nodeptr,dir))==NULL)
      printf("nodeupdate: bad dir");
    if((ndata1 = finddata(nodeptr1,dir1))==NULL)
      {
	printf("nodeupdate: bad dir1");
	early_exit(1833);
      }
    if((ndata2 = finddata(nodeptr2,dir2))==NULL)
      printf("nodeupdate: bad dir2");

    if(nodeupdate(nodeptr1,dir1)<0)
      {
	printf("bad nodeupdate1");
	retval = -1;
      }

    if(nodeupdate(nodeptr2,dir2)<0)
      {
	printf("bad nodeupdate2");
	retval = -1;
      }

    if(nodeptr1 == (*nodeptr).mother)
      {
	length1 = (*nodeptr).length;
      }
    else
      {
	length1 = (*nodeptr1).length;
      }

    if(nodeptr2 == (*nodeptr).mother)
      {
	length2 = (*nodeptr).length;
      }
    else
      {
	length2 = (*nodeptr2).length;
      }

	if(vrenorm(nonnegck(length1,1911), nonnegck(length2,1911),ndata1,ndata2,ndata,calc) < 0)
	  {
	    gamma = rateequil[Pos]; 

	    p1 = pck(exp(-gamma*length1),-1);
	    p2 = pck(exp(-gamma*length2),-2);

	    printf("impossible tree! nodes %d, %d, %d\n",
		   (*nodeptr).seqno,(*nodeptr1).seqno,(*nodeptr2).seqno);
	    printf("length1=%g, length2=%g, Pos= %d\n",length1,length2,Pos);
	    printf("p1= %g, p2= %g gamma= %g\n",p1,p2,gamma);
	    early_exit(1879);
	  }
	(*nodeptr).inc = 'y';

  }

if(retval==0 && (*nodeptr).inc != 'y') retval=1;

return(retval);
}
  

Double_with_error *porig(nodeptr,dir)
Treeptr nodeptr;
char dir;
{
Treeptr nodeptr1;
double p1, omp1, *pcharvect1, *pnvect1, length1, sum1, tmp, errsum, anssum,gamma;
double *virtcts1;
double *virtctsp1;
double negvirtct1, negvirtctp1;
char  dir1;
int maxn1, i;
Renormdata *ndata1;

switch (dir) 
  {
  case 'm':
    nodeptr1 = (*nodeptr).mother;
    dir1 = (*nodeptr).side;
    break;
/*  case 'l':
    nodeptr1 = (*nodeptr).left;
    dir1 = 'm';
    break;
  case 'r':
    nodeptr1 = (*nodeptr).right;
    dir1 = 'm';
    break;      these don't apply*/
  default:
    printf("porig: unknown side-type");
    break;
  }


if(nodeptr1 == (*nodeptr).mother)
  length1 = (*nodeptr).length;
else
  length1 = (*nodeptr1).length;

/*printf("start at %d %c\n",(*nodeptr1).seqno,dir1);*/

gamma = rateequil[Pos];          /* rate of equilibration*/

if(nodeupdate(nodeptr1,dir1)==0)
  {
    ndata1 = finddata(nodeptr1,dir1);
    pcharvect1 = (*ndata1).pcharvect;
    maxn1 = (*ndata1).maxn[Char-'A'];
    if(maxn1>0)
      pnvect1 = (*ndata1).pnvect[Char-'A'];
    virtcts1 = (*ndata1).virtcts;
    virtctsp1 = (*ndata1).virtctsp;

    p1 = pck(exp(-gamma*length1),77);
    omp1 = -expm1(-gamma*length1); /* 1-p1 */
    if( omp1 > 0.0)  /*previously gamma*length1 > 0 */
      {
	negvirtct1 = -gamma*length1*p1/omp1;
	negvirtctp1 = (gamma*length1*p1*(omp1-gamma*length1)/omp1)/omp1;
      }
    else
      {
	negvirtct1 = - 1.0;
	negvirtctp1 = - 0.0;
      }
                   
/* sum1 is not a good name, but: */

    sum1 = (omp1)*oldprob[Pos][Char-'A'] /* pcharvect1[] sum to 1.0 */
      + pcharvect1[1]*p1*oldprob[Pos][Char-'A'];
    anssum = pck(sum1,77);              /* character unshared by inheritance*/

    Virtcts = 
      omp1*
	(virtcts1[25]+negvirtct1) +
      p1*pcharvect1[1]*
	(virtcts1[1]+gamma*length1)          +
      (p1*pcharvect1[Char-'A']/oldprob[Pos][Char-'A'])*
	(virtcts1[Char-'A']+gamma*length1);
    Virtcts /= (omp1+p1*pcharvect1[1]+p1*pcharvect1[Char-'A']/
		oldprob[Pos][Char-'A']);
    Virtctsp = 
      omp1*
	(square(virtcts1[25]+negvirtct1) +
	 virtctsp1[25]+negvirtctp1)
	  +
      p1*pcharvect1[1]*
	  (square(virtcts1[1]+gamma*length1)+
	   (virtctsp1[1]-gamma*length1))
	   +
      (p1*pcharvect1[Char-'A']/oldprob[Pos][Char-'A'])*
	(square(virtcts1[Char-'A']+gamma*length1)+
	 (virtctsp1[Char-'A']-gamma*length1));
    Virtctsp /= (omp1+p1*pcharvect1[1]+p1*pcharvect1[Char-'A']/
		oldprob[Pos][Char-'A']);
    Virtctsp -= square(Virtcts);
    sum1 += pck(pcharvect1[Char-'A']*p1,44);
    Posloglike = log(logck(pck(sum1,66),66)) + pcharvect1[25];
    if(Posloglike > 0.0000001)
      {
	printf("error: likelihood greater than 1! = %g\n",Posloglike);
      }
if(!(Virtctsp <= 0.0 || Virtctsp > 0.0))
  {
    printf("ill-defined Virtctsp = %g\n",Virtctsp);
    printf("Virtcts = %g, p1= %g omp1= %g length1=%g\n",
	   Virtcts, p1, omp1, length1);
    printf("virtctsp1[25]=%g,pcharvect1[1]=%g",
	   virtctsp1[25], pcharvect1[1]);
  }
if(sum1 == 0.0)
  {
    printf("zero likelihood tree, sum1 = 0.0\n");
    printf("Pos = %d, Char = %c, p1 = %f, oldprob = %f, pchar[25] = %f\n",
	   Pos, Char, p1, oldprob[Pos][Char-'A'],pcharvect1[25]);
  }
if(p1==1.0 && sum1 == 0.0)
  {
    printf("at node %d:",(*nodeptr).seqno);
    printf(" zero likelihood at pos=%d\n",Pos);

  }
    anssum = pck(anssum/sum1,55);                 /* normalize */
    
    errsum = anssum*(1.0-anssum);

    for(i=0;i<maxn1;i++)
      { 
/*debug:	if(p1*pnvect1[i]/(sum1) < 1.001)
	  {
	    p1 = p1;
	  }
	else
	  {
	    printf("i= %d maxn1 = %d p1= %g pnvect1[i]= %g sum1=%g\n",
		   i,maxn1,p1,pnvect1[i],sum1);
	  }
*/
	tmp = pck(p1*pnvect1[i]/(sum1),88);
                                                /*shared with i+1 others*/
	anssum += tmp/((double) (i+2));
	errsum += (tmp*(1.0-tmp))/((double) (i+2));
      }
    answer.value = pck(anssum, 99);
    if(!(errsum>=0.0))
      printf("negative? %g ans.val = %g sum1= %g errsum = %g maxn1=%d\n",
	      errsum, answer.value, sum1,
	     errsum,maxn1);
    answer.error = nonnegck(errsum,2077);
    return(&answer);
  }
else
  {
/*    printf("no active sequences\n");*/
    answer.value = 1.0;
    answer.error = 0.0;
    return(&answer);
  }

}

void cleartreedata()
{
int i,j,maxi;
if(NoSeq >=3)
  maxi = 2*NoSeq-2;
else
  maxi = 4;
for(i=0;i<maxi;i++)
  {
    for(j=0;j<26;j++)
      {
	(Top[Treeit][i]).mcalc[j]='n';
	(Top[Treeit][i]).rcalc[j]='n';
	(Top[Treeit][i]).lcalc[j]='n';

	((Top[Treeit][i]).mrdata).maxn[j]=0;
	((Top[Treeit][i]).lrdata).maxn[j]=0;
	((Top[Treeit][i]).rrdata).maxn[j]=0;
	(Top[Treeit][i]).inc = '?';
      }
  }
}

void cleartreechanges()
{
int i,maxi;
if(NoSeq >=3)
  maxi = 2*NoSeq-2;
else
  maxi = 4;
for(i=0;i<maxi;i++)
  {
    (Top[Treeit][i]).changes= 0.0;
    changenorms[i] = 0.0;
  }
}

void poschangesupdate()
{
int i,j,maxi;
double gamma, length, pchange, pnochange, pdontcare;
double *pcharvect1, *pcharvect2;
Renormdata *ndata1, *ndata2;
Treeptr node1, node2;
double tmp;
char calc;

gamma = nonnegck(rateequil[Pos],650);

if(NoSeq >=3)
  maxi = 2*NoSeq-2;
else
  maxi = 4;

for(i=0;i<maxi;i++)
  {
    node1 = &(Top[Treeit][i]);
    if((*node1).mcalc['Z'-'A'] == 'y')
      {
	length = nonnegck((*node1).length,665);
	ndata1 = &((*node1).mrdata);
	pcharvect1 = (*ndata1).pcharvect;
	node2 = (*node1).mother;
	if(node2==NULL)
	  {
	    printf("could not find mother of node!");
	    early_exit(2118);
	  }
	if((*node1).side == 'l')
	  {
	    ndata2 = &((*node2).lrdata);
	    calc = (*node2).lcalc['Z'-'A'];
	  }
	else if((*node1).side == 'r')
	  {
	    ndata2 = &((*node2).rrdata);
	    calc = (*node2).rcalc['Z'-'A'];
	  }
	else if((*node1).side == 'm')
	  {
	    ndata2 = &((*node2).mrdata);
    	    calc = (*node2).mcalc['Z'-'A'];
	  }
	else
	  {
	    printf("unknown side type ! ( %c )\n",(*node1).side);
	    early_exit(2138);
	  }
	  
       if(calc == 'y')
       {
	pcharvect2 = (*ndata2).pcharvect;
	pnochange = 0.0;

	for(j=0;j<25;j++)
	  {
	    if(j != 1 && pcharvect1[j] > 0.0 )
	      {
		if(oldprob[Pos][j]==0.0)
		  {
		    printf("zero prob char at 760: pos %d char %c pchr %g\n",
			   Pos,'A'+j,pcharvect1[j]);
		    early_exit(2154);
		  }
		pnochange += 
		  nonnegck(pcharvect1[j]*pcharvect2[j]/oldprob[Pos][j],760);
	      }
	  }
	pchange = pck((1.0 - pcharvect1[1])*(1.0-pcharvect2[1]),761);
	pdontcare = pck(pcharvect1[1]+pcharvect2[1] -
			pcharvect1[1]*pcharvect2[1],762);
	tmp = pchange*(-expm1(-gamma*length))+ pnochange*exp(-gamma*length);
	if(tmp + pdontcare <=0.0)
	  {
	    fprintf(stderr,"impossible tree!? p= %g\n",
		    pchange*(-expm1(-gamma*length))+ pnochange*exp(-gamma*length));
	    printf("pchange=%g nochange=%g gamma*length= %g\n",
		   pchange, pnochange, gamma*length);
	    printf("-expm1(-gamm*l)= %g ,exp(-gam*l)= %g\n",
		   -expm1(-gamma*length), exp(-gamma*length));
	    printf("p dont care = %g",pdontcare);
	    fflush(stdout);
		
		early_exit(2175);
	  }

/* following is true EM: */
/*	(*node1).changes += nonnegck((pchange+pdontcare)*gamma*length/
	  (pchange*(-expm1(-gamma*length))+ pnochange*exp(-gamma*length)+
	   pdontcare),770);
	changenorms[i] += gamma;
*/
/* following hopefully converges faster */
/*	if(pchange>pnochange)
	  rdontcaremin = pdontcare/(pchange+pdontcare);
	else
	  rdontcaremin = pdontcare/(pnochange+pdontcare);

	(*node1).changes += (pchange+pdontcare)*gamma*length/
	  (pchange*(-expm1(-gamma*length))+ pnochange*exp(-gamma*length)+
	   pdontcare)-rdontcaremin*gamma*length ;
	changenorms[i] += gamma*(1.0-rdontcaremin) ;
*/
/* even faster? */
	(*node1).changes += nonnegck(pchange*gamma*length/
	  (pchange*(-expm1(-gamma*length))+ pnochange*exp(-gamma*length)+
	   pdontcare),770);
	changenorms[i] += nonnegck(gamma*
	  (pchange*(-expm1(-gamma*length))+ pnochange*exp(-gamma*length))/
	  (pchange*(-expm1(-gamma*length))+ pnochange*exp(-gamma*length)+
	   pdontcare),770);

       }				     
      }
  }

}

void poscountsupdate()
{
int pos1, i,k;
char char1;
double norm, p,dp;
double p1,dp1,p2,dp2,poslike1,poslike2;
double virtcts1, virtcts2;
double virtctsp1, virtctsp2;
Double_with_error *tporig, newcounts[26];
double gamma;

pos1 = Pos;
gamma = rateequil[pos1];
cleartreedata();


norm = 20.0*Pscount;    /*pseudocount*/
for(char1='A';char1<='Z';char1++)
{
  newcounts[char1-'A'].value = Pscount;
  newcounts[char1-'A'].error = 0.0;
}

for(i=0;i<NoSeq;i++)
{
  if(treeloc[Treeit][i]==NULL) 
    {
      printf("lost some sequences!!!\n");
      printf("seq i= %d has null tree location\n",i);
      printf("NoSeq= %d\n",NoSeq);
      printf("Treeit= %d\n",Treeit);
      printf("Iter= %d\n",Iter);
      printf("NoTrees= %d\n",NoTrees);
      early_exit(2243);
    }
  else 
    {
      if(seq[i][pos1] != 'X')
	(*treeloc[Treeit][i]).inc = 'y';
      else
	(*treeloc[Treeit][i]).inc = 'n';
    }
}

Virtcts = 0.0; /* default if no active sequences */
Virtctsp = -1.0;
Posloglike = 0.0;

for(k=0;k<NoSeq;k++)
{
  p = 0.0;
  Char = seq[k][Pos];
  if((*treeloc[Treeit][k]).inc == 'y' && Char != 'B' && Char != 'Z') 
    {
      tporig = porig(treeloc[Treeit][k],'m');
      p = (*tporig).value;
      dp = (*tporig).error;
      newcounts[seq[k][pos1]-'A'].value += p;
      newcounts[seq[k][pos1]-'A'].error += dp;
      norm += p;
    }
  if(Char == 'B')
    {
      Char = 'N';
      tporig = porig(treeloc[Treeit][k],'m');
      p1 = (*tporig).value;
      dp1 = (*tporig).error;
      poslike1 = Posloglike;
      virtcts1 = Virtcts;
      virtctsp1 = Virtctsp;
      Char = 'D';
      tporig = porig(treeloc[Treeit][k],'m');
      p2 = (*tporig).value;
      dp2 = (*tporig).error;
      poslike2 = Posloglike;
      virtcts2 = Virtcts;
      virtctsp2 = Virtctsp;
      p = (p1*poslike1+p2*poslike2)/(poslike1+poslike2);
      dp = (dp1*poslike1+dp2*poslike2)/(poslike1+poslike2);
      newcounts['N'-'A'].value += p1*poslike1/(poslike1+poslike2);
      newcounts['N'-'A'].error +=dp1*poslike1/(poslike1+poslike2);
      newcounts['D'-'A'].value += p2*poslike2/(poslike1+poslike2);
      newcounts['D'-'A'].error +=dp2*poslike2/(poslike1+poslike2);
      norm += p;
      if(poslike1 > poslike2)
	{
	  Posloglike = poslike1 + log(logck(1.0+exp(poslike2-poslike1),700));
	}
      else
	{
	  Posloglike = poslike2 + log(logck(1.0+exp(poslike1-poslike2),701));
	}
      Virtcts = (exp(poslike1)*virtcts1 + exp(poslike2)*virtcts2)/
	(exp(poslike1) + exp(poslike2));
      Virtctsp = (exp(poslike1)*virtctsp1 + exp(poslike2)*virtctsp2)/
	(exp(poslike1) + exp(poslike2));
    }
  if(Char == 'Z')
    {
      Char = 'Q';
      tporig = porig(treeloc[Treeit][k],'m');
      p1 = (*tporig).value;
      dp1 = (*tporig).error;
      poslike1 = Posloglike;
      virtcts1 = Virtcts;
      Char = 'E';
      tporig = porig(treeloc[Treeit][k],'m');
      p2 = (*tporig).value;
      dp2 = (*tporig).error;
      poslike2 = Posloglike;
      virtcts2 = Virtcts;
      p = (p1*poslike1+p2*poslike2)/(poslike1+poslike2);
      dp = (dp1*poslike1+dp2*poslike2)/(poslike1+poslike2);
      newcounts['Q'-'A'].value += p1*poslike1/(poslike1+poslike2);
      newcounts['Q'-'A'].error +=dp1*poslike1/(poslike1+poslike2);
      newcounts['E'-'A'].value += p2*poslike2/(poslike1+poslike2);
      newcounts['E'-'A'].error +=dp2*poslike2/(poslike1+poslike2);
      norm += p;
      if(poslike1 > poslike2)
	{
	  Posloglike = poslike1 + log(logck(1.0+exp(poslike2-poslike1),710));
	}
      else
	{
	  Posloglike = poslike2 + log(logck(1.0+exp(poslike1-poslike2),711));
	}
      Virtcts = (exp(poslike1)*virtcts1 + exp(poslike2)*virtcts2)/
	(exp(poslike1) + exp(poslike2));
      Virtctsp = (exp(poslike1)*virtctsp1 + exp(poslike2)*virtctsp2)/
	(exp(poslike1) + exp(poslike2));
      
    }
}

Virtctsp *= gamma; /* change of variables */

for(char1='A';char1<='Z';char1++) 
{
  treecounts[Treeit][pos1][char1-'A'].value = 
    newcounts[char1-'A'].value;
  treecounts[Treeit][pos1][char1-'A'].error = 
    newcounts[char1-'A'].error;
}    



virttreecounts[Treeit][pos1] = Virtcts;
virttreecountsp[Treeit][pos1] = Virtctsp;

}


void treechangesnormalize()
{
int  i, maxi;

if(NoSeq >=3)
  maxi = 2*NoSeq-2;
else
  maxi = 4;

for(i=0;i<maxi;i++)
  {
    if(changenorms[i] > 0.0)
      (Top[Treeit][i]).changes /= changenorms[i];
  }

}

void updatelengths(treeit)
int treeit;
{
int i, maxi;

Treeit = treeit;
if(NoSeq >=3)
  maxi = 2*NoSeq-2;
else
  maxi = 4;

    for(i=0;i<maxi;i++)
      {
/*	if(((Top[Treeit][i]).length < 0.1/((double) NoSeq*LengthSeq) &&
	   (Top[Treeit][i]).changes < (Top[Treeit][i]).length) ||
	   Top[Treeit][i].length == 0.0)
	  {
	    (Top[Treeit][i]).length = 0.0;
	  }
	else
*/
	  (Top[Treeit][i]).length = nonnegck((Top[Treeit][i]).changes,771);

/*printf("%d:%g",i,(Top[Treeit][i]).length);*/
      }
}

void treeout(p,donetop,raw) /* originally from Phylip */
Treeptr p;
char donetop;
int raw;
{
  /* write out file with representation of final tree */
  long  w;
  double x;
  int attop=0,j,k;
  char donetopl;

donetopl = donetop;

if(p == (*(*p).mother).mother && donetop == 0)
  {
    attop=1;
    donetopl = 1;
  }


  if ((*p).seqno>=0) 
    {
      k=0;
      for(j=0;j<NAMELENGTH;j++)
	if(names[(*p).seqno][j] != ' ') k=j;
      for(j=0;j<=k;j++)
	if(names[(*p).seqno][j] == ' ')
	  fprintf(treefile,"_");
	else
	  fprintf(treefile,"%c",names[(*p).seqno][j]);
      Col += NAMELENGTH;
    }
  else 
    {
      putc('(', treefile);
      Col++;
      treeout((*p).left,donetopl,raw);
    putc(',', treefile);
    Col++;
    if (Col > 55) {
      putc('\n', treefile);
      Col = 0;
    }
    treeout((*p).right,donetopl,raw);
    if (attop==1)
      {
	putc(',', treefile);
	treeout((*p).mother,donetopl,raw);
    }
    putc(')', treefile);
    Col++;
  }
  x = ((*p).length)*Rateequil*(1.0-Selfchange);
  if(raw==1) x = (*p).length;
  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (attop == 1)
    fprintf(treefile, ";\n");
  else {
    fprintf(treefile, ":%*.5f", (int)(w + 7), x); 
    Col += w + 8;

  }
}  /* treeout */


void updatelocaltopology(treeit)
int treeit;
{
int i, maxi;
Treeptr nodeo, nodex, tmpptr;

Treeit = treeit;
if(NoSeq >=3)
  maxi = 2*NoSeq-2;
else
  maxi = 4;

    for(i=0;i<maxi;i++)
      {
	if((Top[Treeit][i]).length < .02/LengthSeq 
	   && (Top[Treeit][i]).seqno < 0
	   && ((Top[Treeit][i]).mother)->seqno < 0 
	   && (   ((Top[Treeit][i]).left)->length > 1.0/LengthSeq
	       || ((Top[Treeit][i]).right)->length > 1.0/LengthSeq)
	   && (Top[Treeit][i]).length <= 
	   (1.0-exp(- (((Top[Treeit][i]).left)->length)))*
	   (1.0-exp(- (((Top[Treeit][i]).right)->length))))
	  {
	    nodeo = &(Top[Treeit][i]);
	    nodex = (Top[Treeit][i]).mother;
	    printf("changing topology tree %d nodes %d and %d",Treeit,nodeo->seqno,nodex->seqno);

	    treeout(Top[Treeit],0,0);
	    
	    if(nodeo->side != 'm')
	      {
		printf("(a)\n");
		(nodeo->right)->mother = nodex;
		if(nodeo->side == 'l')
		  nodex->left = nodeo->right;
		else if(nodeo->side == 'r')
		  nodex->right = nodeo->right;
		else
		  printf("unknown side %c!\n",nodeo->side);
		(nodeo->right)->side = nodeo->side;

		nodeo->mother = nodex->mother;
		if(nodex->side == 'l')
		  (nodex->mother)->left = nodeo;
		else if(nodex->side == 'r')
		  (nodex->mother)->right = nodeo;
		else if(nodex->side == 'm')
		  (nodex->mother)->mother = nodeo;
		else
		  printf("unknown side %c!\n",nodex->side);
		nodeo->side = nodex->side;

		nodex->mother = nodeo;
		nodeo->right = nodex;
		nodex->side = 'r';
		
		nodeo->length = nodex->length;
		if(nodeo->length >= 1.0/LengthSeq)
		  nodeo->length-=.01/LengthSeq;
		nodex->length = .03/LengthSeq;
		if((nodeo->left)->length >= 1.0/LengthSeq)
		  (nodeo->left)->length -= .01/LengthSeq;

		if(&(Top[Treeit][0]) == nodex) /* have to move 'root' */
		  {
		    nodeo->right = nodeo->mother;
		    nodeo->mother = nodex;
		    (nodeo->right)->side = 'r';
		    nodex->side = 'm';
		  }


	      }
	    else
	      {
		printf("(b)\n");
		(nodeo->right)->mother = nodex;
		tmpptr = nodex->right;
		nodex->right = nodeo->right;
		
		(nodex->left)->mother = nodeo;
		nodeo->right = nodex->left;
		(nodeo->right)->side = 'r';

		tmpptr->mother = nodex;
		nodex->left = tmpptr;
		(nodex->left)->side = 'l';

		nodeo->length = .03/LengthSeq;
		nodex->length = .03/LengthSeq;
	      }
	  }
      }
}

double calccountschange()
{
int pos1,char1;
double maxchange;

maxchange = 0.0;
for(pos1=0;pos1 < LengthSeq;pos1++)
      {
	for(char1=0;char1<26;char1++) 
	  {
	    if(fabs(counts[pos1][char1].value - 
		    oldcounts[pos1][char1])  >
	       fabs(maxchange))
	      {
	        maxchange = 
		  counts[pos1][char1].value - 
		    oldcounts[pos1][char1];
	      }
	    oldcounts[pos1][char1]=counts[pos1][char1].value;
	  }
      }
return(maxchange);
}

void average()
{  /* average stored counts over all trees */
int i,pos1,char1;
double tmp, maxlike,sumlike,tweight;

maxlike = loglikevect[0];
for(i=0; i < NoTrees; i++)
  if(loglikevect[i]> maxlike) maxlike=loglikevect[i];

sumlike = 0.0;
for(i=0; i < NoTrees; i++)
  if(loglikevect[i]>maxlike+log(EPS))
     sumlike += exp(loglikevect[i]-maxlike);

for(pos1=0;pos1<LengthSeq;pos1++)
  {
    for(char1=0;char1<26;char1++)
      {
	counts[pos1][char1].value = 0.0;
	counts[pos1][char1].error = 0.0;
      }
    virtcounts[pos1] = 0.0;
    virtcountsp[pos1] = 0.0;
  }

for(i=0; i < NoTrees; i++)
  {
    tweight = 0.0;
    if(loglikevect[i]-maxlike> log(EPS))
      tweight = exp(loglikevect[i]-maxlike)/sumlike;
    
    if(tweight>0.0)
      {
	for(pos1=0;pos1<LengthSeq;pos1++)
	  {
	    for(char1=0;char1<26;char1++)
	      {
		tmp = treecounts[i][pos1][char1].value;
		counts[pos1][char1].value += tweight*tmp;
		counts[pos1][char1].error += tweight*tmp*tmp;
		counts[pos1][char1].error += 
		  tweight*tweight*treecounts[i][pos1][char1].error;
		if(treecounts[i][pos1][char1].error < 0) printf("AHA!!\n"); 
	      }
	    virtcounts[pos1] += tweight*virttreecounts[i][pos1];
	    virtcountsp[pos1] += tweight*virttreecountsp[i][pos1];
	  }

      }
  }
if(Convgd==1) printf("counts=%f, err=%f\n",counts[0][22].value,counts[0][22].error);


}

void applyvirtcts()
{
int pos1;
int ntied;
char char1;
double maxcts, smaxcts, ctstoadd, ctssofar;
double lambda, gamma, effectivevirtcts,norm;
Double_with_error *newcounts;

for(pos1=0;pos1<LengthSeq;pos1++)
  {
    newcounts = counts[pos1];
    maxcts = 0.0;
    norm = 0.0;
    for(char1='A';char1<='Z';char1++) 
      {
	if(char1 != 'B')
	  norm += counts[pos1][char1-'A'].value;
      }
    for(char1='A';char1<='Z';char1++) 
      {
	if(nonnegck((newcounts[char1-'A']).value,111) > maxcts) 
	  maxcts = (newcounts[char1-'A']).value;
      }

 if(maxcts > 0.0){
    ntied = 0;
    for(char1='A';char1<='Z';char1++) 
      {
	if((newcounts[char1-'A']).value == maxcts) ntied ++;
      }    
    if(ntied == 0)
      {
	printf("ntied == 0!\n");
	early_exit(3083);
      }
    Virtcts = virtcounts[pos1];
    Virtctsp = virtcountsp[pos1];

    gamma = rateequil[pos1];

    lambda = norm + Virtcts + Virtctsp*(1.0/((double) ntied)-1.0/gamma);

while(square(lambda)-
   4.0*(Virtctsp/((double) ntied))*
   (norm-maxcts*((double) ntied)) < 0.0)
  {
    Virtctsp /= 2.0;  /* total hack!*/
    lambda = norm + Virtcts + Virtctsp*(1.0/((double) ntied)-1.0/gamma);
  }
     lambda = 
     (lambda + sqrt(nonnegck(square(lambda)-
			     4.0*(Virtctsp/((double) ntied))*
			     (norm-maxcts*((double) ntied))     ,911)))/2.0;
/*
else
{
  printf("\n bad input to sqrt, %g, ntied= %d, gamma = %g, Virtctsp=%g",
	 square(lambda)-
	 4.0*(Virtctsp/((double) ntied))*
	 (norm-maxcts*((double) ntied)), ntied, gamma, Virtctsp);
  printf("Pos = %d, Char= %c", Pos, Char);
  printf("norm = %g, maxcts= %g, Virtcts= %g\n", norm, maxcts, Virtcts);
}
*/

effectivevirtcts = 
((lambda/(lambda-Virtctsp/((double) ntied)))*
 (maxcts+(Virtcts-Virtctsp/gamma)/((double) ntied)) - maxcts)*
((double) ntied);

ctssofar = 0.0;

if(effectivevirtcts >= 0.0)
{
  
  for(char1='A';char1<='Z';char1++) 
    {
      if(newcounts[char1-'A'].value == maxcts) 
	{
	  newcounts[char1-'A'].error += 
	    (2.0*newcounts[char1-'A'].value+effectivevirtcts/(double)ntied)*
	      effectivevirtcts/(double) ntied;
	  newcounts[char1-'A'].value += 
	    effectivevirtcts/(double) ntied;
	}
    }
}
else
{
  /*find second largest counts:*/
  smaxcts = 0.0;
  for(char1='A';char1<='Z';char1++) 
    {
      if(newcounts[char1-'A'].value > smaxcts &&
	 newcounts[char1-'A'].value < maxcts )
	smaxcts = newcounts[char1-'A'].value;
    }

  ctstoadd = Virtcts;
  while(effectivevirtcts/((double) ntied) + maxcts < smaxcts && smaxcts > 0.0)
    {
      ntied = 0;
      for(char1='A';char1<='Z';char1++) 
	{
	  if(newcounts[char1-'A'].value >= smaxcts)  
	    {
	      ctstoadd -= (smaxcts-newcounts[char1-'A'].value);
	      norm += (smaxcts-newcounts[char1-'A'].value);
	      newcounts[char1-'A'].error += 
		(2.0*newcounts[char1-'A'].value+(smaxcts-newcounts[char1-'A'].value))*
		  (smaxcts-newcounts[char1-'A'].value);
	      newcounts[char1-'A'].value = smaxcts;
	      ntied ++;
	    }
	}
      if(ntied ==0)
	{
	  printf("ntied == 0!\n");
	  early_exit(3117);
	}
      ctssofar = Virtcts-ctstoadd;
      maxcts = smaxcts;

      lambda = norm + ctstoadd + Virtctsp*(1.0/((double) ntied)-1.0/gamma);
      lambda = 
	(lambda + sqrt(nonnegck(square(lambda)-
				4.0*(Virtctsp/((double) ntied))*
				(norm-maxcts*((double) ntied))    ,913)))/2.0;

      effectivevirtcts = 
	((lambda/(lambda-Virtctsp/((double) ntied)))*
	 (maxcts+(ctstoadd-Virtctsp/gamma)/((double) ntied)) - maxcts)*
	   ((double) ntied);

      smaxcts = 0.0;
      for(char1='A';char1<='Z';char1++) 
	{
	  if(newcounts[char1-'A'].value > smaxcts &&
	     newcounts[char1-'A'].value < maxcts )
	    smaxcts = newcounts[char1-'A'].value;
	}

    } /* end while*/
  if(effectivevirtcts <= 0.0)
    {
      for(char1='A';char1<='Z';char1++) 
	{
	  if(newcounts[char1-'A'].value == maxcts)  
	    {
	      if(smaxcts > 0.0) 
		{
		  newcounts[char1-'A'].error += 
		    (2.0*newcounts[char1-'A'].value+effectivevirtcts/(double)ntied)*
		      effectivevirtcts/(double) ntied;
		  newcounts[char1-'A'].value += 
		    effectivevirtcts/(double) ntied;
		}
	      else
		{
		  virtcounts[pos1] -= effectivevirtcts/(double)ntied;
 		}
	    }
	}
    }
  else
    {
      printf("\n bistable virtcts: Virtcts= %g, Virtctsp = %g, eff=%g\n",
	     Virtcts,Virtctsp,effectivevirtcts);
      printf("lambda = %g, ntied = %d \n",lambda, ntied);
      early_exit(3168);
    }
} /*end else*/

} /*end if */

} /*end for*/

}

void savecounts()
{
int pos1, char1;
double tmp;

for(pos1=0;pos1<LengthSeq;pos1++)
  {
    for(char1=0;char1<26;char1++)
      {
	tmp = counts[pos1][char1].value;
	prevcounts[pos1][char1] =
	  tmp;
	counts[pos1][char1].error -= tmp*tmp;
	if(counts[pos1][char1].error < 0.0)
	  {
	    if(counts[pos1][char1].error < -0.000001)
	      {
		printf("negative variance!\n");
		printf("char1=%d, counts=%f, esq =%f\n",char1,tmp,
		       counts[pos1][char1].error+tmp*tmp);
		early_exit(3194);
	      }
	    else
	      counts[pos1][char1].error = 0.0;
	  }
      }
  }
}


void calcprobs()
{
int pos1, char1;
double p,pct,norm;

for(pos1=0;pos1<LengthSeq;pos1++)
  {
    norm = 0.0;
    for(char1=0;char1<25;char1++)
      {
	norm += counts[pos1][char1].value;
      }
    norm += 20.0; /*pseudocounts */
    for(char1=0;char1<25;char1++)
      {
	pct = 1.0; /* pseudo count */
	if(char1 == 'J' -'A') pct = 0.0;
	p = (counts[pos1][char1].value + pct)/(norm+EPS);

	probs[pos1][char1].value = p;
	probs[pos1][char1].error = 
	  (counts[pos1][char1].error + 
	   (counts[pos1][char1].value+pct)*
	   (norm -(counts[pos1][char1].value+pct))/ norm)/
	     (norm*norm);
      }
    probs[pos1]['B'-'A'].value = 0.0; 
    probs[pos1]['O'-'A'].value = 0.0;
    probs[pos1]['U'-'A'].value = 0.0;
    probs[pos1]['X'-'A'].value = 0.0; 
    probs[pos1]['Z'-'A'].error = 0.0;
    probs[pos1]['B'-'A'].error = 0.0; 
    probs[pos1]['O'-'A'].error = 0.0;
    probs[pos1]['U'-'A'].error = 0.0;
    probs[pos1]['X'-'A'].error = 0.0; 
    probs[pos1]['Z'-'A'].error = 0.0;

  }

}

void calcentropy()
{
int pos1, char1;
double p,pct,entro,norm;

for(pos1=0;pos1<LengthSeq;pos1++)
  {
    entro = 0.0;
    norm = 0.0;
    for(char1=0;char1<25;char1++)
      norm += counts[pos1][char1].value;
    for(char1=0;char1<26;char1++)
      {
	pct = 0.0; /* no pseudo count for entropy*/
	if(char1 == 'J' -'A') pct = 0.0;
	p = (counts[pos1][char1].value + pct)/(norm+EPS);

	if(p > 0.0 && p < 1.0) entro -= p*log(logck(p,222))/log(2.0);
      }
    entropy[pos1] = entro;
  }

}

void printprobs()
{
FILE *fp;
int pos1, i;

fp = Fopen("results","w");

for(i=0;i<26;i++)
  fprintf(fp,"%-24c ",'A'+i);
fprintf(fp,"\n");
for(pos1=0;pos1<LengthSeq;pos1++)
  {
    fprintf(fp,"position %d:\n",pos1);
    for(i=0;i<25;i++)
      fprintf(fp,"%-10.6f +- %-10.6f ",probs[pos1][i].value,sqrt(nonnegck(probs[pos1][i].error,813)));
      fprintf(fp,"%-10.6f +- %-10.6f ",0.0,0.0); /* for Z */
    fprintf(fp,"\n");
  }

fclose(fp);

}	    

void printmodel()
{
FILE *fp;
int pos1, i;

fp = Fopen("model","w");

for(i=0;i<26;i++)
  {
    if(i!=1 && i!=9 && i!= 14 && i!= 20 && i!= 23)
      fprintf(fp,"%-10c ",'A'+i);
  }
fprintf(fp,"\n");
for(pos1=0;pos1<LengthSeq;pos1++)
  {
    fprintf(fp,"position %d : gamma= %g gapprob= %g oldprobs:\n",pos1,rateequil[pos1],oldprob[pos1][9]);
    for(i=0;i<25;i++)
      {
	if(i!=1 && i!=9 && i!= 14 && i!= 20 && i!= 23)
	  fprintf(fp,"%-10.6f ",oldprob[pos1][i]);
      }
    fprintf(fp,"\n");
  }

fclose(fp);

}	    

void printent()
{
FILE *fp;
 int pos1;

fp = Fopen("entropy","w");

for(pos1=0;pos1<LengthSeq;pos1++)
  {
    fprintf(fp,"%3d %-10.6f ",pos1,entropy[pos1]); /* for Z */
    fprintf(fp,"\n");
  }

fclose(fp);

}	    

void printvirtcts()
{
FILE *fp;
 int pos1;

fp = Fopen("virtcts","w");


for(pos1=0;pos1<LengthSeq;pos1++)
  {
    fprintf(fp,"%3d     ",pos1); 
    fprintf(fp,"%10g ",virtcounts[pos1]);
    fprintf(fp,"\n");
  }

fclose(fp);

}

void printcounts()
{
FILE *fp;
int pos1, i;

fp = Fopen("counts","w");
/*fprintf(fp,"counts: (no virtual or pseudo counts included)\n");*/
for(i=0;i<25;i++)
  {
    if(i!=1 && i!=9 && i!= 14 && i!= 20 && i!= 23)
       fprintf(fp,"%-15c ",'A'+i);
  }
fprintf(fp,"\n");
for(pos1=0;pos1<LengthSeq;pos1++)
  {
    fprintf(fp,"position %d :\n",pos1);
    for(i=0;i<25;i++)
      {
	if(i!=1 && i!=9 && i!= 14 && i!= 20 && i!= 23)
	  fprintf(fp,"%-6.3f +- %-5.3f ",counts[pos1][i].value,sqrt(nonnegck(counts[pos1][i].error,811)));
      }

    fprintf(fp,"\n");
  }
fclose(fp);
}

void printgaps()
{
FILE *fp;
int pos1;

fp = Fopen("gaps","w");

for(pos1=0;pos1<LengthSeq;pos1++)
  fprintf(fp,"%-6.3f +- %-5.3f\n",counts[pos1][9].value,sqrt(nonnegck(counts[pos1][9].error,811)));

fclose(fp);
}


void modelupdate(ltol,improve) /*  updates oldprob and changes */
double ltol, improve;
{
int pos1, i;  
char char1;
double norm, likechange;
double loglikelihood;
int numit;
int lasttime=0;
double oldloglike, maxloglike, sumloglike,countschange, loglikeinit;
/*double oldposloglike[MaxLengthSeq];*/

printf("modelupdate\n");

numit = 0;
likechange = ltol + 1.0;
countschange = COUNTTOL + 1.0;
oldcountsinit();
sumloglike = loglikeinit = 0.0;

while(((likechange > ltol && 
       fabs(countschange) > COUNTTOL && 
       sumloglike - loglikeinit < improve &&
       numit < MAXUPIT) ||
      numit < 2) || lasttime ==0)
  {
    if(!((likechange > ltol && 
       fabs(countschange) > COUNTTOL && 
       sumloglike - loglikeinit < improve &&
       numit < MAXUPIT) ||
      numit < 2)) lasttime = 1;
    numit ++;
    for(Treeit = 0; Treeit < NoTrees; Treeit++)
      {
	loglikelihood = 0.0;
	cleartreechanges();

	for(pos1=0;pos1<LengthSeq;pos1++) 
	  {   
	    Pos = pos1;
	    poscountsupdate();
	    loglikelihood += Posloglike;
/*	    if(Treeit ==1 && numit > 1)
	      {
		if(Posloglike+.000001 < oldposloglike[pos1])
		  {
		    printf("likelihood decrease:Pos= %d, old= %g, new= %g\n",
			   Pos,oldposloglike[pos1],Posloglike);
		  }
	      }
	    if(Treeit == 1)
	      oldposloglike[pos1] = Posloglike;
*/	    
	    poschangesupdate();
	  }   

	loglikevect[Treeit] = loglikelihood;
	treechangesnormalize();

      }  /*done all trees */


    for(Treeit = 0; Treeit < NoTrees; Treeit++)
      {
	if(UserTree < 3)
	  updatelengths(Treeit);
	if(UserTree < 2)
	  updatelocaltopology(Treeit);
      }
    average(); /* average olprob over trees */
    if(lasttime ==1 && Convgd ==1)
      {
	savecounts(); /* adjust errors */
	printcounts();
	printgaps();
      }
    applyvirtcts();  /* take virt cts into account*/
    

    for(pos1=0;pos1<LengthSeq;pos1++)
      {
	norm = 20.0*Pscount;
	for(char1='A';char1<='Z';char1++) 
	  {
	    norm += counts[pos1][char1-'A'].value;
	  }
	norm = nonnegck(norm,6805);
	if(norm==0.0)
	  {
	    if(numit==1 && Iter==1)
	      {
		printf("warning, no counts at pos %d\n",pos1);
		printf("virtcounts: %g, virtctsp = %g\n",virtcounts[pos1],
		       virtcountsp[pos1]);
	      }
	    norm=1.0;
	  }
	for(char1='A';char1<='Z';char1++) 
	  {
	    oldprob[pos1][char1-'A'] = (counts[pos1][char1-'A'].value+Pscount)/norm;
	  }
	oldprob[pos1][1] = 1.0; /*uninheritance has no gain */
      }

    
    rateequilinit();
    
    
    maxloglike = loglikevect[0];
    for(Treeit = 0; Treeit < NoTrees; Treeit++)
      if(loglikevect[Treeit] > maxloglike) maxloglike = loglikevect[Treeit];
    sumloglike = 0.0;
    for(Treeit = 0; Treeit < NoTrees; Treeit++)
      if(loglikevect[Treeit] - maxloglike > log(EPS))
	sumloglike += exp((loglikevect[Treeit]-maxloglike));
    sumloglike = log(sumloglike) + maxloglike;
    countschange = calccountschange();
    
	    printf("sumloglike= %g dcounts = %g\n loglikes=",
		   sumloglike,countschange);
	    for(i=0;i<NoTrees &&i<TREEMAX;i++)
	      printf(" %g, ",loglikevect[i % TREEMAX]);
	    printf("\n ");

    
    if(numit >= 2)
      {
	if( sumloglike + LTOL < oldloglike && Pscount == 0.0)
	  {
	    printf("danger: likelihood decrease in EM!\n");
	    printf("sumloglike %g, oldloglike %g\n",sumloglike,oldloglike);
	    printf("numit = %d\n", numit);
	    printf("\n loglikes=");
	    for(i=0;i<NoTrees &&i<TREEMAX;i++)
	      printf(" %g, ",loglikevect[i % TREEMAX]);
	    printf("\n ");
	  }
	likechange = fabs(sumloglike - oldloglike);
      }
    if(numit == 1) loglikeinit = sumloglike;
    oldloglike = sumloglike;
    
  } /* end while */

printf(" total loglike = %g (%d its)\n",sumloglike, numit);


}

void printrenormdata(ndata)
Renormdata *ndata;
{
int char1, i;  

printf("pcharvect:\n");
 for(char1=0;char1<25;char1++) 
  if((*ndata).pcharvect[char1-'A']>0.0)
    {
      printf("%c: %f \n",char1,(*ndata).pcharvect[char1-'A']);
      for(i=0;i< (*ndata).maxn[char1-'A'];i++)
	printf("%d: %g ",i+1,(*ndata).pnvect[char1-'A'][i]);
      printf("\n");
    }
printf("\n");
}


void printlike()
{
FILE *fp;
int i;
double maxlike;

maxlike = loglikevect[0];
for(i=0; i < NoTrees; i++)
  if(loglikevect[i]> maxlike) maxlike=loglikevect[i];

fp = Fopen("likelihood","w");

fprintf(fp,"loglikelihood= %g\n",maxlike);

fclose(fp);
}

int comparetrees(tree1, tree2)  
Treeptr tree1, tree2; /*will find identical trees if numbered the same */
{
int j;
int numnodes;

numnodes = 2*NoSeq - 2;

j=0;
while(j < numnodes &&
      tree1[j].seqno == tree2[j].seqno && 
      (*tree1[j].mother).seqno == (*tree2[j].mother).seqno) 
  j++;

if(j== numnodes) return(1); /* they match */
return(0);

}

int checktreenovelty()
{
int ttree;

for(ttree=(Treeit+1)%NoTrees;ttree != Treeit; ttree=(ttree+1)%NoTrees)
  {
    if(comparetrees(Top[ttree],Top[Treeit]) == 1)
      {
	printf("\n exact match:");
	printf(" (same tree as %d)\n",ttree);
	return(ttree);
      }
  }
return(-1);
}

int checktreeconvergence()
{
int pos1, char1,  poschange, charchange;
double maxchange, minmaxchange;
int ttree, tmp;

/* see if tree we just made matches any previous tree */
maxchange = 1.0 + FINALTOL;
minmaxchange = 1000.0;

tmp = checktreenovelty();
if(tmp >= 0)
  return(tmp);

for(ttree=(Treeit+1)%NoTrees;ttree != Treeit; ttree=(ttree+1)%NoTrees)
  {
    maxchange = 0.0;
    pos1=0;
    while(pos1 < LengthSeq && fabs(maxchange) < FINALTOL)
      {
	for(char1=0;char1<26;char1++) 
	  {
	    if(fabs(treecounts[Treeit][pos1][char1].value - 
		    treecounts[ttree][pos1][char1].value)  >
	       fabs(maxchange))
	      {
		maxchange = 
		  treecounts[Treeit][pos1][char1].value - 
		  treecounts[ttree][pos1][char1].value;
		poschange = pos1;
		charchange = char1;
	      }
	  }
	if(fabs(virttreecounts[Treeit][pos1] - 
		virttreecounts[ttree][pos1])  >
	   fabs(maxchange))
	  {
	    maxchange = 
	      virttreecounts[Treeit][pos1] - 
		virttreecounts[ttree][pos1];
	    poschange = pos1;
	    charchange = 1;
	  }
	
	pos1 ++;
      }
    if(fabs(maxchange) < minmaxchange)
      minmaxchange = fabs(maxchange);
  }

if(fabs(maxchange) < FINALTOL) /*converged!*/
  {
    printf("\n close to tree %d \n",ttree);
    return(ttree);
  }
else
  {
    printf("minmaxchange=%g \n",minmaxchange);
    return(-1);
  }
}

int checkconvergence()
{
int pos1, char1;
double maxchange, minmaxchange;

/* see if counts matches previous iteration */

maxchange = 1.0 + FINALTOL;
minmaxchange = 0.0;

pos1=0;
maxchange = 0.0;
while(pos1 < LengthSeq && fabs(maxchange) < FINALTOL)
      {
	for(char1=0;char1<26;char1++) 
	  {
	    if(fabs(counts[pos1][char1].value - 
		    prevcounts[pos1][char1])  >
	       fabs(maxchange))
	      {
	        maxchange = 
		  counts[pos1][char1].value - 
		  prevcounts[pos1][char1];
	      }
	  }
	pos1 ++;
      }


if(fabs(maxchange) < FINALTOL) /*converged!*/
  {
    printf("\n converged! \n");
    return(1);
  }
else
  {
    printf("\n maxchange = %g",maxchange);
  return(-1);
  }
}

void improvenewtree(ltol)
double ltol;
{
int pos1, i;  
double  likechange;
double loglikelihood, bestloglike;
int numit;
double oldloglike, extraploglike;
/*double oldposloglike[MaxLengthSeq];*/

numit = 0;
likechange = ltol + 1.0;

bestloglike = loglikevect[0];
for(i=0;i<NoTrees;i++)
  if(i != Treeit && loglikevect[i] > bestloglike)
    bestloglike = loglikevect[i];
loglikelihood = 2.0*bestloglike - 1.0;
extraploglike = bestloglike;

while(likechange > ltol && numit < MAXUPIT && loglikelihood < bestloglike &&
      (extraploglike > bestloglike - LIKECUTOFF || numit < MINTREEIMPIT))
  {
    numit ++;

    loglikelihood = 0.0;
    cleartreechanges();
    
    for(pos1=0;pos1<LengthSeq;pos1++) 
      {   
	Pos = pos1;
	poscountsupdate();
	loglikelihood += Posloglike;
	poschangesupdate();
      }   

    loglikevect[Treeit] = loglikelihood;
    treechangesnormalize();

    if(UserTree < 3)
      updatelengths(Treeit);
    if(UserTree < 2)
      updatelocaltopology(Treeit);
    printf("numit=%d, loglike=%16.10g\n",numit,loglikelihood);
    if(numit >= 2)
      {
	if( loglikelihood + LTOL < oldloglike && Pscount == 0.0)
	  {
	    printf("warning: likelihood decrease in EM!\n");
/*	    printf("loglike %g, oldloglike %g\n",loglikelihood,oldloglike);
	    printf("numit = %d\n", numit);
	    printf("\n loglikes=");
	    for(i=0;i<NoTrees;i++)
	      printf(" %g, ",loglikevect[i]);
	    printf("\n ");
	    early_exit(2774);
*/
	    oldloglike = loglikelihood;
	  }
	likechange = loglikelihood - oldloglike;
	extraploglike = loglikelihood + likechange*(double)(MAXUPIT - numit);
	if(Pscount != 0.0)
	  likechange = fabs(likechange);
      }
    oldloglike = loglikelihood;
    
  } /* end while */

if(likechange <= ltol || loglikelihood >= bestloglike)
  printf(" new tree %d optimized in %d its\n", Treeit, numit);
else
  printf(" gave up improving new tree %d after %d its\n", Treeit, numit);

}

void newtree(tol)
double tol;
{
int i, worsti;
double minlike;
Treeptr tmpptr;
char command[70], numberstring[22];

tmpptr=allocatetree();

Treeit = NoTrees;

if(tmpptr==NULL)
  {
    if(Iter==0)
      {
	printf("not enough memory for tree!\n");
	printf("exiting...\n");
	early_exit(2808);
      }
    else if(Iter<2)
      {
	printf("not enough memory for two trees!\n");
	printf("try reducing MaxNoSeq or MaxLengthSeq\n");
	printf("exiting...\n");
	early_exit(2809);
      }
    else
      {
	printf("out of memory, replacing worst tree\n");
      }
  }
if(tmpptr==NULL||Treeit==TREEMAX)
  {
    minlike = loglikevect[0];
    worsti = 0;
    
    for(i=0; i < NoTrees; i++)
      {
	if(loglikevect[i]< minlike) 
	  {
	    minlike=loglikevect[i];
	    worsti = i;
	  }
      }
    Treeit = worsti;
  }
else
  {
    Top[NoTrees] = tmpptr;
    NoTrees ++;
  }

if(UserTree==0)
  writedistances();

if(NoTrees==1 && UserTree==0)
  {
    if(TreeScript == 0)
      {
	printf("\nWill use weighbor to build trees\n");
	system("which weighbor");
      }
    else if(NoSeq > MaxSlowSeqs )
      {
	printf("\nMore than MaxSlowSeqs seqs; will use treescript=:\n");
	system("cat treescript");
      }
    else if(zerodist == 1 )
      {
	printf("\nBecause of identical seqs, will use treescript=:\n");
	system("cat treescript");
      }
    else if( NoSeq < 4 )
      {
	printf("\nBecause less than 4 seqs, will use treescript=:\n");
	system("cat treescript");
      }
    else
      {
	printf("\nFewer than MaxSlowSeqs seqs; will use slowtreescript=:\n");
	system("cat slowtreescript");
      }
  }

if(TreeScript == 0 && UserTree==0)
  {
    strncpy(command,"weighbor -b ",70);
    sprintf(numberstring,"%20.14lg ",1/Selfchange);
    strncat(command,numberstring,21);
    strncat(command,"-L ",15);
    sprintf(numberstring,"%20.14lg ",Eff_Length);
    strncat(command,numberstring,21);
    strcat(command,"< infile > treefile ");
    
    
    printf("command = %s\n", command);
    if(system(command) != 0)
      {
	printf("could not execute weighbor\n");
	early_exit(3148);
      }
  }
else if(UserTree==0 && (NoSeq > 16 || NoSeq < 4 || zerodist == 1))
  {
    if(system("treescript")!= 0)
      {
	if(NoTrees==1)
	  printf("could not execute treescript; trying neighbor\n");
	if(system("echo y | nnneighbor > treescript.out")!=0)
	  {
	    if(NoTrees==1)
	      printf("nnneighbor didn't work; trying neighbor\n");
	    if(system("echo y | neighbor > treescript.out")!=0)
	      {
		if(NoTrees==1)
		  printf("Hopefully you already have a valid tree installed in treefile.\n");
	      }
	  }
      }
  }
else if(UserTree==0)
  {
    if(system("slowtreescript")!= 0)
      {
	if(NoTrees==1)
	  printf("could not execute slowtreescript; trying treescript\n");
	if(system("treescript")!= 0)
	  {
	    if(NoTrees==1)
	      printf("could not execute slowtreescript; trying fitch\n");
	    if(system("echo y | fitch > slowtreescript.out")!=0)
	      {
		if(NoTrees==1)
		  printf("trying nnneighbor\n");
		if(system("echo y | nnneighbor > slowtree")!=0)
		  {
		    if(NoTrees==1)
		      printf("Hopefully you already have a valid tree installed in treefile.\n");
		  }
	      }
	  }
      }
  }

if((treefile=Fopen("treefile","r"))==NULL)
  {
    printf("can't open treefile for reading!\n");
    early_exit(3013);
  }

for(i=0;i<NoSeq;i++)
      treeloc[Treeit][i] = NULL;


if(ReadTree(Top[Treeit], treefile)!=0) 
  {
    printf("couldn't parse tree!\n");
    early_exit(3023);
  }
printf("\n");

fclose(treefile);

if(UserTree == 0 && checktreenovelty() < 0 && Iter > 1)
  {
    improvenewtree(tol);
  }

}

void printbesttree()
{
int i, besti;
double maxlike;

treefile = Fopen("besttree","w");
Col = 0;
maxlike = loglikevect[0];
besti = 0;

for(i=0; i < NoTrees; i++)
  {
    if(loglikevect[i]> maxlike) 
      {
	maxlike=loglikevect[i];
	besti = i;
      }
  }

treeout(Top[besti],0,0);

fclose(treefile);

treefile = Fopen("besttree.raw","w");
treeout(Top[besti],0,1);
fclose(treefile);

}

void printbesttrees()
{
int i, besti, sbesti;
double maxlike;

treefile = Fopen("besttree","w");
Col = 0;

maxlike = loglikevect[0];
besti = 0;

for(i=0; i < NoTrees; i++)
  {
    if(loglikevect[i]> maxlike) 
      {
	maxlike=loglikevect[i];
	besti = i;
      }
  }

treeout(Top[besti],0,0);

fclose(treefile);

treefile = Fopen("besttree.raw","w");
treeout(Top[besti],0,1);
fclose(treefile);

/* now 2nd best: */

treefile = Fopen("2ndbesttree","w");
Col = 0;

if(besti == 0)
  {
    maxlike = loglikevect[1];
    sbesti = 1;
  }
else
  {
    maxlike = loglikevect[0];
    sbesti = 0;
  }    

for(i=0; i < NoTrees; i++)
  {
    if(loglikevect[i]> maxlike && i != besti) 
      {
	maxlike=loglikevect[i];
	sbesti = i;
      }
  }

treeout(Top[sbesti],0,0);

fclose(treefile);
}

void cpnode(fromnode,tonode)
Treeptr fromnode, tonode;
{
(*tonode).seqno = (*fromnode).seqno;
(*tonode).length = (*fromnode).length ;
if((*tonode).seqno <0)
{
  (*tonode).left = tonode + ((*fromnode).left-fromnode);
  (*tonode).right = tonode + ((*fromnode).right-fromnode);
}
else
{
  (*tonode).left = NULL;
  (*tonode).right = NULL;
}
(*tonode).mother = tonode + ((*fromnode).mother-fromnode);
(*tonode).side = (*fromnode).side ;
}


void cptree(fromtree,totree)
Treeptr fromtree, totree;
{
int j;
int numnodes;

numnodes = 2*NoSeq - 2;

for(j=0;j<numnodes;j++)
  {
    cpnode(&(fromtree[j]),&(totree[j]));
  }
for(j=0;j<numnodes;j++)
  {
    if((*fromtree[j].mother).seqno != (*totree[j].mother).seqno)
      {
	printf("cptree failed\n");
	early_exit(3539);
      }
  }

}


int iterate()
{
double  bestlike, worstlike;
int  i, convgd, worsti;

/* should enter with Iter = 2 */

convgd = -1;

printf("iterating: (max %d iterations)",MAINITMAX);
fflush(stdout);

while(convgd < 0 &&  Iter  < MAINITMAX)
{
printf("\n%d: ", Iter+1);

newtree(exp(-.1*(double)(Iter*Iter+1)));

convgd = checktreeconvergence();

if(convgd <0)  /*if new tree same as an old one will forget new tree */
  {
    Iter ++;

    modelupdate(log(1.0+1.0/(double)Iter), REALIMPROVEMENT);

    printbesttree();

    convgd = checkconvergence();

    bestlike = loglikevect[0];
    worstlike = bestlike;
    worsti = 0;
    for(i=1; i < NoTrees; i++)
      {
	if(loglikevect[i] > bestlike) 
	  bestlike = loglikevect[i];
	if(loglikevect[i] < worstlike)
	  {
	    worstlike = loglikevect[i];
	    worsti = i;
	  }
      }
    if(NoTrees > TREEMIN && worstlike - bestlike < -LIKECUTOFF)
      {            
	if(worsti != NoTrees - 1) /* discard  worst tree */
	  {
	    loglikevect[worsti] = loglikevect[NoTrees-1];
	    cptree(Top[NoTrees - 1], Top[worsti]);

	    for(i=0;i<NoSeq; i++)
	      treeloc[worsti][i] = Top[worsti]+
		(treeloc[NoTrees-1][i]-Top[NoTrees-1]);
	  }
	if(worsti == NoTrees -1) /* new tree has no effect, so stop */
	  convgd = NoTrees;

	NoTrees --;  
	printf("removed tree %d, leaving %d trees\n",worsti,NoTrees);
      }
  }
else /*if new tree same as an old one will remove new tree */
  {
    NoTrees --;
    printf("removed last tree, leaving %d trees\n",NoTrees);
  }

fflush(stdout);

}   /* end of while */


return(convgd);
}



void oldprobinit()
{

oldprobinit1();

Iter = 0;
NoTrees = 0;

if(UserTree==0)
  {
    printf("\nmaking first tree ");

    fflush(stdout);

    newtree(log(2.0)*(double) LengthSeq);

    fflush(stdout);
  }
else
  {
    printf("reading user tree from treefile\n");
    newtree(log(2.0)*(double) LengthSeq);
  }

Iter = 1;

if(UserTree==0)
  {

    oldprobinit2();

    printf("\nmaking second tree ");

    fflush(stdout);

    newtree(log(2.0)*(double) LengthSeq);

    printf("\n");
    fflush(stdout);

    /* now have 2 trees ready */
    Iter = 2;

    oldprobiniteq();
    modelupdate(log(2.0), REALIMPROVEMENT);

    Treeit = 0;
    improvenewtree(log(2.0));
    Treeit = 1;
    improvenewtree(log(2.0));

    fflush(stdout);

    printbesttree();
  }
}





void usage(char *argv)
{
printf("unknown option %s\n",argv);
printf("allowed options: \n");
printf("-s     execute \"treescript\" to build tree (instead of weighbor)\n");
printf("-v     read data verbosely\n");
printf("-u     read tree from treefile; do not calculate distance matrix\n");
printf("-U     same as -u but also disallow local rearrangements of tree\n");
printf("-T     same as -U but also disallows branch length changes\n");
printf("-P cts apply cts pseudocounts to each residue\n");
exit(1);
}


int main(int argc, char *argv[])
{
extern char *optarg;
extern int optind;
extern int optopt;
extern int opterr;
extern int optreset;

char *fn = "data";
int convgd,i,ch;
printf("rind version 6.9.5 2009\n");

  /*
  command-line options:
   */
Verbose = 0;
UserTree = 0;
TreeScript = 0;
while ((ch = getopt(argc, argv, "svuUTP:")) != -1) {
  switch(ch) {
  case 's':
    TreeScript = 1;
    break;
  case 'v':
    Verbose = 1;
    break;
  case 'u': /* read tree from treefile; do local rearrangements */
    UserTree = 1; 
    break;
  case 'U': /* read tree from treefile; no rearrangments */
    UserTree = 2; 
    break;
  case 'T': /* read tree from treefile; no rearrangments, fixed lengths */
    UserTree = 3; 
    break;
  case 'P':
    Pscount = atof(optarg);
    if (isnan(Pscount) || !isfinite(Pscount)) {
      fprintf(stderr, "Improper number of pseudocounts: %s\n", optarg);
      fprintf(stderr, "Option -%c requires a numerical argument.\n", ch);
      exit(1);
    }
    break;
  default:
    usage((char *)&ch);
  }
}
for (i = optind; i < argc; i++) {
  fn = argv[i];
}

Convgd=0;
readdata(fn);


printf("initializing...\n");
oldprobinit();
prevcountsinit();

if(UserTree==0)
  {
    convgd = iterate();  /* all the action is here */
  }

Convgd = 1;
modelupdate(LTOL,BIGNO);

if(UserTree==0)
  {
    newtree(LTOL);

    convgd = checktreeconvergence();

    if(convgd < 0)
      {
	Iter ++;
	modelupdate(LTOL,BIGNO);
      }
    else
      {
	NoTrees --;
      }
  }
printf("convergence complete\n");

if(UserTree==0)
  printbesttrees();
else
  printbesttree();

printmodel();

calcentropy();

printent();

printvirtcts();

printlike();

for(i=0;i< NoTrees;i++)
free(Top[i]);

return 0;

}
