#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util.h"

void dmatinv (double *a, int m) {
  double b[m*m];

  double temp, big, pivot, rmax;
  int i, iu, j, k, l = 0, jl, ib, ir, kk;
  
  /* Initialise b as a unit matrix */
  for(j = 0; j < m; j++) {
    for(i = 0; i < m; i++)
      b[j*m+i] = 0.0;

    b[j*m+j] = 1.0;
  }

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
      for(ib = 0; ib < m; ib++) b[ib*m+ib] = 0.0;
/*        fprintf(stderr, "matinv: Zero determinant\n"); */
      return;
    }

    if(i != l) {
      /* switch rows */
      for(j = 0; j < m; j++) {
	temp    = a[j*m+i];
	a[j*m+i] = a[j*m+l];
	a[j*m+l] = temp;
	temp    = b[j*m+i];
	b[j*m+i] = b[j*m+l];
	b[j*m+l] = temp;
      }
    }

    /* pivotal reduction */
    pivot = a[i*m+i];
    jl = i+1;

    for(j = jl; j < m; j++) {
      temp = a[i*m+j]/pivot;
      for(k = i; k < m; k++)
	a[k*m+j] -= temp*a[k*m+i];

      for(k = 0; k < m; k++)
	b[k*m+j] -= temp*b[k*m+i];
    }
  }

  /* back substitution for solution */
  for(j = 0; j < m; j++)
    for(i = 0; i < m; i++) {
      ir = m-1-i;
      temp = b[j*m+ir];
      if(ir != m-1) {
	for(k = 1; k <= i; k++) {
	  kk = m - k;
	  temp -= a[kk*m+ir] * b[j*m+kk];
	}
      }

      b[j*m+ir] = temp / a[ir*m+ir];
    }

  /* copy it back into a */
  memcpy(a, b, sizeof(b));
}

