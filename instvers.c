#include <sys/types.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>

#include "lightcurves.h"

#include "util.h"

int read_instvers (char *filename, struct instvers **instverslist_r, int *ninstvers_r,
		   char *errstr) {
  FILE *fp;
  char line[16384], *p;
  int rv;

  /* These make the syntax a bit more sane */
  struct instvers *instverslist = *instverslist_r;
  int ninstvers = *ninstvers_r;

  /* Open file */
  fp = fopen(filename, "r");
  if(!fp) {
    report_syserr(errstr, "open: %s", filename);
    goto error;
  }

  /* Read lines */
  while(fgets(line, sizeof(line), fp)) {
    p = sstrip(line);

    if(*p == '\0')
      continue;

    ninstvers++;
    instverslist = (struct instvers *) realloc(instverslist,
					       ninstvers * sizeof(struct instvers));
    if(!instverslist) {
      report_syserr(errstr, "realloc");
      goto error;
    }

    rv = sscanf(p, "%d %ld",
		&(instverslist[ninstvers-1].iver),
		&(instverslist[ninstvers-1].date));
    if(rv != 2) {
      report_err(errstr, "could not understand: %s", p);
      goto error;
    }
  }

  if(ferror(fp)) {
    report_syserr(errstr, "read");
    goto error;
  }

  fclose(fp);

  *instverslist_r = instverslist;
  *ninstvers_r = ninstvers;

  return(0);

 error:
  return(1);  /* XXX - memory leak */
}
