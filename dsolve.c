#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"

/* gauss elimination to solve ax=b */

void dsolve (double *a, double *b, int m) {
  double temp, big, pivot, rmax;
  int i, iu, j, k, l = 0, jl, ib, ir;

  iu = m-1;
  for(i = 0; i < iu; i++) {
    big = 0.0;

    /* find largest remaining term in ith column for pivot */
    for(k = i; k < m; k++) {
      rmax = fabsf(a[i*m+k]);
      if(rmax > big) {
	big = rmax;
	l = k;
      }
    }

    /* check for non-zero term */
    if(big == 0.0) {
      for(ib = 0; ib < m; ib++) b[ib] = 0.0;
/*        fprintf(stderr, "solve: Zero determinant\n"); */
      return;
    }

    if(i != l) {
      /* switch rows */
      for(j = 0; j < m; j++) {
	temp    = a[j*m+i];
	a[j*m+i] = a[j*m+l];
	a[j*m+l] = temp;
      }
      temp = b[i];
      b[i] = b[l];
      b[l] = temp;
    }

    /* pivotal reduction */
    pivot = a[i*m+i];
    jl = i+1;

    for(j = jl; j < m; j++) {
      temp = a[i*m+j]/pivot;
      b[j] -= temp*b[i];
      for(k = i; k < m; k++) a[k*m+j] -= temp*a[k*m+i];
    }
  }

  /* back substitution for solution */
  for(i = 0; i < m; i++) {
    ir = m-1-i;
    if(a[ir*m+ir] != 0.0) {
      temp = b[ir];
      if(ir != m-1) {
	for(j = 1; j <= i; j++) {
	  k = m-j;
	  temp -= a[k*m+ir]*b[k];
	}
      }
      b[ir] = temp/a[ir*m+ir];
    }
    else
      b[ir] = 0.0;
  }
}
