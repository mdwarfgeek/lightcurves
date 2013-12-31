#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "lightcurves.h"
#include "util.h"

char **read_file_list (int argc, char **argv, long *nf_r, char *errstr) {
  char **fnlist = (char **) NULL, **fnp;

  FILE *fp;
  char line[16384], *p;

  long op, nf, a, f;

  op = 0;
  nf = 0;
  for(a = 1; a < argc; a++) {
    if(*(argv[a]) == '@') {
      /* @list form */
      fp = fopen(argv[a] + 1, "r");
      if(!fp) {
	report_syserr(errstr, "open: %s", argv[a] + 1);
	goto error;
      }
      
      /* Count number of lines */
      while(fgets(line, sizeof(line), fp)) {
	p = sstrip(line);
	
	if(*p != '\0')
	  nf++;
      }
      
      if(ferror(fp))
	error(1, "%s: read", argv[a] + 1);
      
      rewind(fp);

      /* Allocate buffer space */
      fnlist = (char **) realloc(fnlist, nf * sizeof(char *));
      if(!fnlist)
	error(1, "realloc");

      fnp = fnlist + op;

      f = 0;
      while(fgets(line, sizeof(line), fp)) {
	p = sstrip(line);
	
	if(*p != '\0') {
	  *fnp = strdup(p);
	  if(!*fnp)
	    error(1, "strdup");
	  
	  fnp++;
	  f++;
	}
      }

      op += f;

      if(ferror(fp))
	error(1, "%s: read", argv[a] + 1);

      if(op != nf)
	fatal(1, "unexpected number of lines in %s: expected %d, got %d", argv[a]+1, nf, f);

      fclose(fp);
    }
    else {
      /* Single filename */
      nf++;

      fnlist = (char **) realloc(fnlist, nf * sizeof(char *));
      if(!fnlist)
	error(1, "realloc");

      fnp = fnlist + op;

      *fnp = strdup(argv[a]);
      if(!*fnp)
	error(1, "strdup");

      op++;
    }
  }

  *nf_r = nf;

  return(fnlist);

 error:
  /* we bomb anyway do no need to worry about leaks */

  return((char **) NULL);
}
