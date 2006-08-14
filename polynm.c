#include <stdio.h>
#include <stdlib.h>

#include "floatmath.h"
#include "util.h"

extern void dsolve (double a[50][50], double b[50], int m);

int polynm (float *xbuf, float *ybuf, float *wtbuf, long npt,
	    float *coeff, long ncoeff, float *chisq_r, char *errstr) {
  double a[50][50], b[50];
  long pt, k, j;

  float tmp, wt, ymod, chisq;

  /* Sanity check sizes of arrays, etc. */
  if(npt < ncoeff) {
    report_syserr(errstr, "too few data points (%ld) for degree %ld", npt, ncoeff);
    goto error;
  }

  if(ncoeff > 50) {
    report_syserr(errstr, "degree of polynomial (%ld) too large, max 50", ncoeff);
    goto error;
  }

  /* Clear arrays */
  for(k = 0; k < ncoeff; k++) {
    for(j = 0; j < ncoeff; j++)
      a[k][j] = 0.0;

    b[k] = 0.0;
  }

  /* Sum */
  for(pt = 0; pt < npt; pt++) {
    wt = (wtbuf ? wtbuf[pt] : 1.0);

    for(k = 0; k < ncoeff; k++) {
      for(j = 0; j <= k; j++)
	a[k][j] += powf(xbuf[pt], j+k) * wt;

      b[k] += ybuf[pt] * powf(xbuf[pt], k) * wt;
    }
  }

  /* Fill in duplicates */
  for(k = 0; k < ncoeff-1; k++)
    for(j = k+1; j < ncoeff; j++)
      a[k][j] = a[j][k];

  /* Solve */
  dsolve(a, b, ncoeff);

  /* Copy out answer */
  if(coeff)
    for(k = 0; k < ncoeff; k++)
      coeff[k] = b[k];

  if(chisq_r) {
    /* Calculate chisq */
    chisq = 0.0;
    
    for(pt = 0; pt < npt; pt++) {
      wt = (wtbuf ? wtbuf[pt] : 1.0);
      
      ymod = 0.0;
      for(k = 0; k < ncoeff; k++)
	ymod += b[k] * powf(xbuf[pt], k);
      
      tmp = ymod - ybuf[pt];
      chisq += wt * tmp*tmp;
    }

    *chisq_r = chisq;
  }

  return(0);

 error:
  return(1);
}
