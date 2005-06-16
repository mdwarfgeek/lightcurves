#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>

#include "lightcurves.h"

#include "cvtunit.h"
#include "floatmath.h"
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
  if(buffer_alloc(&buf, 6, 3, 2, errstr))
    fatal(1, "buffer_alloc: %s", errstr);

  for(i = 0; i < 3; i++) {
    data[i].flux = i;
    data[i].fluxerr = i*2;
    data[i].satur = i % 2;
  }

  /* Write out */
  if(buffer_put_object(&buf, data, 0, 3, 2, 1, errstr))
    fatal(1, "buffer_put_object: %s", errstr);

  /* Read back */
  for(k = 0; k < 6; k++)
    for(j = 0; j < 2; j++) {
      if(buffer_fetch_object(&buf, newdata, 0, 3, k, j, errstr))
	fatal(1, "buffer_fetch_object: %s", errstr);

      for(i = 0; i < 3; i++)
	printf("%d %d %f %f %d\n", k, j,
	       newdata[i].flux, newdata[i].fluxerr, newdata[i].satur);
    }

  /* Release disk buffer */
  buffer_close(&buf);

  return(0);
}

