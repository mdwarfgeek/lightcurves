#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lightcurves.h"

#include "util.h"

void medsig (float *a, long n, float *median_r, float *sigma_r) {
  long m;
  float median, mad;

  if(n == 0) {
    /* Special case - can't do anything */
    *median_r = 0.0;
    if(sigma_r)
      *sigma_r = 0.0;
    return;
  }

  fquicksort(a, n);
  median = (n % 2 == 0 ? 0.5 * (a[n/2-1] + a[n/2]) : a[n/2]);

  *median_r = median;

  if(sigma_r) {
    for(m = 0; m < n; m++)
      a[m] = fabsf(a[m] - median);
    
    fquicksort(a, n);
    mad = (n % 2 == 0 ? 0.5 * (a[n/2-1] + a[n/2]) : a[n/2]);

    *sigma_r = 1.48 * mad;
  }
}
