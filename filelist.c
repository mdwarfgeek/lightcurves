#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifndef _WIN32
/* No glob on Win32, but we can get it done for us by linking with
   setargv.obj.  On my MinGW install this (or some equivalent) seems
   to be the default anyway. */
#include <glob.h>
#endif

#include "lightcurves.h"
#include "util.h"

char **read_file_list (int argc, char **argv, int *nf_r, char *errstr) {
  char **fnlist = (char **) NULL;

  FILE *fp = (FILE *) NULL;
  char line[16384], *av, *p, *s;

  int a, f, nf = 0;

#ifndef _WIN32
  glob_t gbuf;
  int rv, op;

  /* Init */
  gbuf.gl_pathv = (char **) NULL;
#endif

  /* Loop over arguments */
  for(a = 0; a < argc; a++) {
    av = argv[a];

    if(*av == '@') {  /* @list form */
      av++;  /* skip past the @ */

      fp = fopen(av, "r");
      if(!fp) {
	report_syserr(errstr, "open: %s", av);
	goto error;
      }
      
      while(fgets(line, sizeof(line), fp)) {
	p = sstrip(line);

        if(*p) {
          fnlist = (char **) realloc(fnlist, (nf+1) * sizeof(char *));
          if(!fnlist) {
            report_syserr(errstr, "realloc");
            goto error;
          }
          
          s = strdup(p);
          if(!s) {
            report_syserr(errstr, "strdup");
            goto error;
          }
          
          fnlist[nf] = s;
          nf++;
        }
      }

      if(ferror(fp)) {
        report_syserr(errstr, "read: %s", av);
        goto error;
      }

      fclose(fp);
      fp = (FILE *) NULL;
    }
    else {  /* glob or single filename */
#ifndef _WIN32
      rv = glob(av, 0, NULL, &gbuf);
      if(rv == 0) {  /* succeeded */
        /* Allocate block */
        fnlist = (char **) realloc(fnlist,
                                   (nf+gbuf.gl_pathc) * sizeof(char *));
        if(!fnlist) {
          report_syserr(errstr, "realloc");
          goto error;
        }

        /* Zero out the new pointers */
        memset(fnlist+nf, 0, gbuf.gl_pathc * sizeof(char *));

        /* Record so we know to free them */
        op = nf;

        nf += gbuf.gl_pathc;
        
        for(f = 0; f < gbuf.gl_pathc; f++) {
          s = strdup(gbuf.gl_pathv[f]);
          if(!s) {
            report_syserr(errstr, "strdup");
            goto error;
          }

          fnlist[op] = s;
          op++;
        }

        globfree(&gbuf);
        gbuf.gl_pathv = (char **) NULL;
      }
      else if(rv == GLOB_NOMATCH) {  /* no match */
#endif  /* _WIN32 */

        /* Assume it's a single entry. */
        fnlist = (char **) realloc(fnlist, (nf+1) * sizeof(char *));
        if(!fnlist) {
          report_syserr(errstr, "realloc");
          goto error;
        }
        
        s = strdup(av);
        if(!s) {
          report_syserr(errstr, "strdup");
          goto error;
        }
        
        fnlist[nf] = s;
        nf++;
#ifndef _WIN32
      }
      else {
        report_err(errstr, "glob error: %s", av);
        goto error;
      }
#endif
    }
  }

  *nf_r = nf;

  return(fnlist);

 error:
  if(fp)
    fclose(fp);

#ifndef _WIN32
  if(gbuf.gl_pathv)
    globfree(&gbuf);
#endif

  if(fnlist) {
    for(f = 0; f < nf; f++)
      if(fnlist[f])
        free((void *) fnlist[f]);

    free((void *) fnlist);
  }

  return((char **) NULL);
}
