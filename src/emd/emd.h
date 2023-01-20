#ifndef _EMD_H
#define _EMD_H
/*
    emd.h

    Last update: 3/24/98

    An implementation of the Earth Movers Distance.
    Based of the solution for the Transportation problem as described in
    "Introduction to Mathematical Programming" by F. S. Hillier and
    G. J. Lieberman, McGraw-Hill, 1990.

    Copyright (C) 1998 Yossi Rubner
    Computer Science Department, Stanford University
    E-Mail: rubner@cs.stanford.edu   URL: http://vision.stanford.edu/~rubner

    Retrieved from http://robotics.stanford.edu/~rubner/emd/default.htm
    Dec 2022.
*/


/* DEFINITIONS */
#define MAX_SIG_SIZE   5000
#define MAX_ITERATIONS 2000
#define INFINITY_      1e20
#define EPSILON        1e-6

/*****************************************************************************/
/* feature_t SHOULD BE MODIFIED BY THE USER TO REFLECT THE FEATURE TYPE      */
typedef int feature_t;
/*****************************************************************************/


typedef struct
{
  int n;                /* Number of features in the signature */
  feature_t *Features;  /* Pointer to the features vector */
  double *Weights;       /* Pointer to the weights of the features */
} signature_t;


typedef struct
{
  int from;             /* Feature number in signature 1 */
  int to;               /* Feature number in signature 2 */
  double amount;         /* Amount of flow from "from" to "to" */
} flow_t;


double emdwrap(int n, int mp, double *p, double *wp,
	       int mq, double *q, double *wq, int df);

double emd(signature_t *Signature1, signature_t *Signature2,
	  double (*func)(feature_t *, feature_t *),
	  flow_t *Flow, int *FlowSize);

#endif
