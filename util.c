#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <ctype.h>
#include <errno.h>

#include <fitsio.h>

#include "util.h"

char progname[PATH_MAX+1] = { '\0' };

void setprogname (const char *name) {
  (void) strncpy(progname, name, PATH_MAX);
  progname[PATH_MAX] = '\0';
}

void fatal (int code, const char *fmt, ...) {
  va_list ap;

  if(progname[0])
    (void) fprintf(stderr, "%s: ", progname);
  va_start(ap, fmt);
  if(fmt)
    (void) vfprintf(stderr, fmt, ap);
  va_end(ap);
  (void) fprintf(stderr, "\n");

  exit(code);
}

void warning (const char *fmt, ...) {
  va_list ap;

  if(progname[0])
    (void) fprintf(stderr, "%s: ", progname);
  va_start(ap, fmt);
  if(fmt)
    (void) vfprintf(stderr, fmt, ap);
  va_end(ap);
  (void) fprintf(stderr, "\n");
}

void error (int code, const char *fmt, ...) {
  va_list ap;

  if(progname[0])
    (void) fprintf(stderr, "%s: ", progname);
  va_start(ap, fmt);
  if(fmt) {
    (void) vfprintf(stderr, fmt, ap);
    (void) fprintf(stderr, ": ");
  }
  va_end(ap);
  (void) fprintf(stderr, "%s\n", strerror(errno));

  exit(code);
}

void fitsio_err (char *errstr, int status, const char *fmt, ...) {
  char errmsg[FLEN_STATUS];
  va_list ap;
  int rv;

  ffgerr(status, errmsg);
  
  va_start(ap, fmt);
  rv = vsnprintf(errstr, ERRSTR_LEN, fmt, ap);
  va_end(ap);
  
  if(rv != -1 && rv < (ERRSTR_LEN - 1))
    (void) snprintf(errstr + rv, ERRSTR_LEN - rv, ": %s", errmsg);
}          

void report_err (char *errstr, const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  (void) vsnprintf(errstr, ERRSTR_LEN, fmt, ap);
  va_end(ap);
}

void report_syserr (char *errstr, const char *fmt, ...) {
  va_list ap;
  int rv;

  va_start(ap, fmt);
  rv = vsnprintf(errstr, ERRSTR_LEN, fmt, ap);
  va_end(ap);
  
  if(rv != -1 && rv < ERRSTR_LEN)
    (void) snprintf(errstr + rv, ERRSTR_LEN - rv, ": %s", strerror(errno));
}

char *sstrip (char *str) {
  char *p;

  /* Stop at the first comment character */
  p = strchr(str, '#');
  if(p)
    *p = '\0';

  /* First remove whitespace from start of string */
  while(*str != '\0' && isspace((unsigned char) *str))
    str++;

  if(*str == '\0')
    return(str);

  /* Remove whitespace from end of string */
  p = str + strlen(str) - 1;

  while(p > str && isspace((unsigned char) *p))
    p--;

  if(p == str && isspace((unsigned char) *p))
    *p = '\0';
  else
    *(p+1) = '\0';

  return(str);
}

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
