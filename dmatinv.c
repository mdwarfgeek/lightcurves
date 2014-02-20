#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util.h"

void dmatinv (double a[50][50], int m) {
  double b[50][50];

  double temp, big, pivot, rmax;
  int i, iu, j, k, l = 0, jl, ib, ir, kk;
  
  /* Initialise b as a unit matrix */
  for(j = 0; j < m; j++) {
    for(i = 0; i < m; i++)
      b[j][i] = 0.0;

    b[j][j] = 1.0;
  }

  iu = m-1;
  for(i = 0; i < iu; i++) {
    big = 0.0;

    /* find largest remaining term in ith column for pivot */
    for(k = i; k < m; k++) {
      rmax = fabsf(a[i][k]);
      if(rmax > big) {
	big = rmax;
	l = k;
      }
    }

    /* check for non-zero term */
    if(big == 0.0) {
      for(ib = 0; ib < m; ib++) b[ib][ib] = 0.0;
/*        fprintf(stderr, "matinv: Zero determinant\n"); */
      return;
    }

    if(i != l) {
      /* switch rows */
      for(j = 0; j < m; j++) {
	temp    = a[j][i];
	a[j][i] = a[j][l];
	a[j][l] = temp;
	temp    = b[j][i];
	b[j][i] = b[j][l];
	b[j][l] = temp;
      }
    }

    /* pivotal reduction */
    pivot = a[i][i];
    jl = i+1;

    for(j = jl; j < m; j++) {
      temp = a[i][j]/pivot;
      for(k = i; k < m; k++)
	a[k][j] -= temp*a[k][i];

      for(k = 0; k < m; k++)
	b[k][j] -= temp*b[k][i];
    }
  }

  /* back substitution for solution */
  for(j = 0; j < m; j++)
    for(i = 0; i < m; i++) {
      ir = m-1-i;
      temp = b[j][ir];
      if(ir != m-1) {
	for(k = 1; k <= i; k++) {
	  kk = m - k;
	  temp -= a[kk][ir] * b[j][kk];
	}
      }

      b[j][ir] = temp / a[ir][ir];
    }

  /* copy it back into a */
  memcpy(a, b, sizeof(b));
}

