#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>

#include "lightcurves.h"

#include "cvtunit.h"
#include "util.h"

int main (int argc, char *argv[]) {
  struct buffer_info buf;
  struct lc_point data[3], newdata[3];
  char errstr[ERRSTR_LEN];
  int i, j, k;

  /* Create disk buffer */
  if(buffer_init(&buf, errstr))
    fatal(1, "buffer_init: %s", errstr);

  /* Get disk buffer */
  if(buffer_alloc(&buf, 6, 3, errstr))
    fatal(1, "buffer_alloc: %s", errstr);

  for(i = 0; i < 3; i++) {
    data[i].aper[0].flux = i;
    data[i].aper[0].fluxerr = i*2;
    data[i].satur = i % 2;
  }

  /* Write out */
  if(buffer_put_object(&buf, data, 0, 3, 2, errstr))
    fatal(1, "buffer_put_object: %s", errstr);

  /* Read back */
  for(k = 0; k < 6; k++) {
    if(buffer_fetch_object(&buf, newdata, 0, 3, k, errstr))
      fatal(1, "buffer_fetch_object: %s", errstr);

    for(j = 0; j < NFLUX; j++) {
      for(i = 0; i < 3; i++)
	printf("%d %d %f %f %d\n", k, j,
	       newdata[i].aper[j].flux, newdata[i].aper[j].fluxerr, newdata[i].satur);
    }
  }

  /* Release disk buffer */
  buffer_close(&buf);

  return(0);
}

