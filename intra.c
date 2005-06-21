#include <sys/types.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>

#include "lightcurves.h"

#include "floatmath.h"
#include "util.h"

static char *sstrip (char *str);

int read_intra (char *filename, struct intra *intralist, int nmefs,
		char *errstr) {
  FILE *fp;
  char line[16384], *p, *ep;

  int reading_header = 1;
  int cur_row = 0, bin;

  int mef = 0;
  long nbin = 0;
  float *map = (float *) NULL;

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

    if(reading_header) {
      mef = (int) strtol(p, &ep, 10);
      if(!isspace(*ep) || mef < 1) {
	report_err(errstr, "could not understand header: %s", p);
	goto error;
      }

      p = ep;

      nbin = strtol(p, &ep, 10);
      if(*ep != '\0' || nbin < 1) {
	report_err(errstr, "could not understand header: %s", p);
	goto error;
      }

      /* Check MEF */
      if(mef > nmefs) {
	report_err(errstr, "inconsistent number of MEFs: %d %d", mef, nmefs);
	goto error;
      }

      /* Allocate array */
      map = (float *) malloc(nbin*nbin * sizeof(float));
      if(!map) {
	report_syserr(errstr, "malloc");
	goto error;
      }

      reading_header = 0;
      cur_row = 0;
    }
    else {
      /* Read in this line of values */
      for(bin = 0; bin < nbin; bin++) {
	/* Read one */
	map[cur_row*nbin+bin] = strtod(p, &ep);
	if((bin == nbin-1 && *ep != '\0') ||
	   (bin < nbin-1 && !isspace(*ep))) {
	  report_err(errstr, "could not understand line: %s");
	  goto error;
	}

	p = ep;
      }

      cur_row++;
      if(cur_row >= nbin) {
	/* Finished */
	reading_header = 1;

	intralist[mef-1].map = map;
	intralist[mef-1].nbin = nbin;
	intralist[mef-1].binsize = 1.0 / nbin;
      }
    }
  }

  if(ferror(fp)) {
    report_syserr(errstr, "read");
    goto error;
  }

  fclose(fp);

  return(0);

 error:
  return(1);  /* XXX - memory leak */
}

float calc_intra (float x, float y, struct intra *corr) {
  float fracx, fracy, delx, dely, val1, val2, val;
  int xbin, ybin, arg1, arg2, arg3, arg4;

  /* Calculate fractional parts of pixel values */
  fracx = x - NINT(x);
  fracy = y - NINT(y);

  /* Find bins */
  xbin = (int) ((fracx + 0.5) / corr->binsize);
  ybin = (int) ((fracy + 0.5) / corr->binsize);

  assert(xbin >= 0 && xbin < corr->nbin);
  assert(ybin >= 0 && ybin < corr->nbin);

  /* Bilinear interpolation, use the centre of the bin */
  arg1 = ybin * corr->nbin + xbin;
  arg2 = ybin * corr->nbin + xbin+1;
  arg3 = (ybin+1) * corr->nbin + xbin;
  arg4 = (ybin+1) * corr->nbin + xbin+1;

  delx = fracx - ((xbin+0.5) * corr->binsize - 0.5);
  dely = fracy - ((ybin+0.5) * corr->binsize - 0.5);

  val1 = (1.0-delx) * corr->map[arg1] + delx * corr->map[arg2];
  val2 = (1.0-delx) * corr->map[arg3] + delx * corr->map[arg4];
  val = (1.0-dely) * val1 + dely * val2;

  return(val);
}

static char *sstrip (char *str) {
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
