#include <stdlib.h>

#include "util.h"

/* performs linear filtering on array xbuf */

int hanning (float xbuf[], int npt, char *errstr) {
  float *ybuf = (float *) NULL;
  float sum = 0.0, xmns, xmnf;
  int nfilt = 7, i, il, ilow, nelem;

  if(npt <= nfilt)
    return(0);

  /* set first and last edges equal */
  il   = nfilt/2;
  ilow = MAX(3,nfilt/4);
  ilow = (ilow/2)*2 + 1;

  for(i = 0; i < ilow; i++)
    sum += xbuf[i];

  xmns = sum/((float) ilow);

  sum=0.0;
  for(i = 0; i < ilow; i++)
    sum += xbuf[npt-1-i];

  xmnf = sum/((float) ilow);

  /* allocate ybuf array */
  nelem = npt + nfilt;  /* Max. number of elements req'd */

  ybuf = (float *) malloc(nelem * sizeof(float));
  if(!ybuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* reflect edges before filtering */
  for(i = 0; i < il; i++) {
    ybuf[i] = 2.0 * xmns - xbuf[il+ilow-1-i];
    ybuf[npt+i+il] = 2.0 * xmnf - xbuf[npt-i-ilow-1];
  }

  for(i = 0; i < npt; i++)
    ybuf[i+il] = xbuf[i];

  /* do linear filtering on rest */
  for(i = 0; i < npt; i++)
    xbuf[i] = (ybuf[i] + 6.0 * ybuf[i+1] + 15.0 * ybuf[i+2] + 20.0 * ybuf[i+3] +
	       15.0 * ybuf[i+4] + 6.0 * ybuf[i+5] + ybuf[i+6]) / 64.0;

  //xbuf[i] = (ybuf[i] + 4.0 * ybuf[i+1] + 6.0 * ybuf[i+2] +
  //	       4.0 * ybuf[i+3] + ybuf[i+4]) / 16.0;

  free((void *) ybuf);
  ybuf = (float *) NULL;

  return(0);

 error:
  if(ybuf)
    free((void *) ybuf);

  return(1);
}
