#include <sys/types.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <unistd.h>

#include <fitsio.h>
#include <cpgplot.h>

#include "lightcurves.h"
#include "hjd.h"

#include "sla.h"
#include "cvtunit.h"
#include "floatmath.h"
#include "util.h"

struct table_section {
  long instart;
  long outstart;
  long ncopy;
};

static int update_lc (fitsfile *reff, fitsfile *fits,
		      struct buffer_info *buf, struct lc_mef *mefinfo,
		      int outcls, int wantoutcls,
		      char *errstr);

/* Getopt stuff */
extern char *optarg;
extern int optind;

/* Verbose flag */
int verbose = 1;

static void usage (char *av) {
  fprintf(stderr, "Usage:\t%s [options] reffile file [...]\n\n", av);
  fprintf(stderr,
	  "Lightcurve processing:\n"
	  "         -i file   Apply intrapixel correction from 'file'.\n"
	  "         -s level  Override saturation level to 'level'.\n"
	  "Output:\n"
	  "         -c cls    Write out only class==cls (e.g. to select just targets).\n"
	  "         -o file   Writes updated lightcurves to 'file'.\n"
	  "         -u        Updates lightcurves in-place.\n"
	  "         -p        Disables plots.\n"
	  "         -q        Decreases the verbosity level of the program.\n"
  	  "         -v        Increases the verbosity level of the program.\n");
  exit(1);
}  

int main (int argc, char *argv[]) {
  char *pn = (char *) NULL, *avzero, *ep;
  int c;

  char errstr[ERRSTR_LEN];

  char *refname, **fnlist = (char **) NULL;
  struct lc_mef *meflist = (struct lc_mef *) NULL;
  struct intra *intralist = (struct intra *) NULL;
  struct buffer_info buf;

  long f, nf = 0;

  fitsfile *inf, *outf, *tmplf;
  int status = 0, ext, mef, nmefs;

  char outfile[FLEN_FILENAME-1], tmpfile[FLEN_FILENAME-1], fnbuf[FLEN_FILENAME];
  int dooutput = 0;
  int doreplace = 0;
  int outcls = 0;
  int wantoutcls = 0;
  int fd = -1, rv;

  char intrafile[FLEN_FILENAME];
  int dointra = 0;

  int len, maxflen, fspc;
  float *medbuf1 = (float *) NULL, *medbuf2, medsat, medlim;
  long nmedsat, nmedlim, nstartot;

  int diffmode = 0;
  int noplots = 0;

  float satlev = -1.0;
  float sysulim = -1.0, sysllim = -1.0;

  /* Set the program name for error reporting */
  if(argv[0])
    pn = basename(argv[0]);
  
  if(!pn)
    pn = "update";

  setprogname(pn);

  avzero = argv[0];

  /* Extract command-line arguments */
  while((c = getopt(argc, argv, "c:i:o:s:upqv")) != -1)
    switch(c) {
    case 'c':
      outcls = (int) strtol(optarg, &ep, 0);
      if(*ep != '\0')
	fatal(1, "invalid class flag: %s", optarg);
      wantoutcls = 1;
      break;
    case 'i':
      strncpy(intrafile, optarg, sizeof(intrafile)-1);
      intrafile[sizeof(intrafile)-1] = '\0';
      dointra = 1;
      break;
    case 'o':
      strncpy(outfile, optarg, sizeof(outfile)-1);
      outfile[sizeof(outfile)-1] = '\0';
      dooutput = 1;
      break;
    case 's':
      satlev = (float) strtod(optarg, &ep);
      if(*ep != '\0' || satlev < 0)
	fatal(1, "invalid satlev value: %s", optarg);
      break;
    case 'u':
      doreplace = 1;
      break;
    case 'p':
      noplots++;
      break;
    case 'q':
      verbose--;
      break;
    case 'v':
      verbose++;
      break;
    case '?':
    default:
      usage(avzero);
    }

  argc -= optind;
  argv += optind;

  if(verbose < 0)
    verbose = 0;

  /* Get arguments */
  if(argc < 2)
    usage(avzero);

  refname = argv[0];

  fnlist = read_file_list(argc, argv, &nf, errstr);
  if(!fnlist)
    fatal(1, "%s", errstr);

  /* Make stdout unbuffered */
  setvbuf(stdout, (char *) NULL, _IONBF, 0);

  /* Find longest filename */
  maxflen = 0;

  for(f = 0; f < nf; f++) {
    len = strlen(fnlist[f]);
    if(len > maxflen)
      maxflen = len;
  }

  /* How much space do we need to display file numbers? */
  fspc = ceil(log10f(nf));

  /* Open reference file */
  ffopen(&inf, refname, READONLY, &status);
  if(status) {
    fitsio_err(errstr, status, "ffopen: %s", refname);
    fatal(1, "%s", errstr);
  }

  /* Check if we were told which extension to use */
  ffextn(refname, &ext, &status);
  if(status) {
    fitsio_err(errstr, status, "ffextn: %s", refname);
    fatal(1, "%s", errstr);
  }

  if(ext == -99) {
    /* Nope, get number of HDUs and process all of them */
    ffthdu(inf, &nmefs, &status);
    if(status) {
      fitsio_err(errstr, status, "ffthdu: %s", refname);
      fatal(1, "%s", errstr);
    }

    nmefs--;  /* no PHDU */
  }
  else
    nmefs = 1;

  /* Special case when updates of separate output file are requested */
  tmplf = inf;

  if(doreplace && dooutput) {
    /* First, check if it exists, if not it's just the normal output case */
    rv = access(outfile, F_OK);
    if(rv != 0)
      doreplace = 0;
    else {
      /* It's there, better open it then */
      ffopen(&tmplf, outfile, READONLY, &status);
      if(status) {
	fitsio_err(errstr, status, "ffopen: %s", refname);
	fatal(1, "%s", errstr);
      }
    }
  }

  if(doreplace) {
    /* Create temporary file for in-place edit */
    snprintf(tmpfile, sizeof(tmpfile), "%s_XXXXXX", progname);

    fd = mkstemp(tmpfile);
    if(fd == -1)
      error(1, "mkstemp: %s", tmpfile);

    /* Form output name */
    fnbuf[0] = '!';
    strncpy(&(fnbuf[1]), tmpfile, sizeof(fnbuf)-1);
    fnbuf[sizeof(fnbuf)-1] = '\0';
  }
  else {
    /* Form normal output name */
    fnbuf[0] = '!';
    strncpy(&(fnbuf[1]), outfile, sizeof(fnbuf)-1);
    fnbuf[sizeof(fnbuf)-1] = '\0';
  }

  /* Create output file */
  if(dooutput || doreplace) {
    /* Create file */
    ffinit(&outf, fnbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffinit: %s", &(fnbuf[1]));
      fatal(1, "%s", errstr);
    }

    /* Close file descriptor if it was temporary - we have it now */
    if(doreplace)
      close(fd);

    /* Copy in PHDU */
    ffcopy(inf, outf, 0, &status);
    if(status) {
      fitsio_err(errstr, status, "ffcopy");
      fatal(1, "%s", errstr);
    }
  }

  /* Allocate arrays */
  meflist = (struct lc_mef *) malloc(nmefs * sizeof(struct lc_mef));
  intralist = (struct intra *) malloc(nmefs * sizeof(struct intra));
  medbuf1 = (float *) malloc(2 * nmefs * sizeof(float));
  if(!meflist || !intralist || !medbuf1)
    error(1, "malloc");

  medbuf2 = medbuf1 + nmefs;

  /* Initialise intrapixel list */
  for(mef = 0; mef < nmefs; mef++)
    intralist[mef].map = (float *) NULL;

  /* Load intrapixel list if requested */
  if(dointra) {
    if(read_intra(intrafile, intralist, nmefs, errstr))
      fatal(1, "read_intra: %s", errstr);
  }

  /* Create disk buffer */
  if(buffer_init(&buf, errstr))
    fatal(1, "buffer_init: %s", errstr);

  nmedsat = 0;
  nmedlim = 0;
  nstartot = 0;

  /* Initialise PGPLOT for diagnostics */
#ifdef DEBUG
  cpgopen("?");
  cpgsubp(3, 1);
  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);
  cpgsch(1.4);
#endif

  for(mef = 0; mef < nmefs; mef++) {
    if(verbose)
      printf("Processing MEF %d of %d\n", mef+1, nmefs);

    /* Read in reference info */
    meflist[mef].stars = (struct lc_star *) NULL;
    meflist[mef].nstars = 0;

    meflist[mef].avsigma = 0.0;
    meflist[mef].avapcor = 0.0;
    meflist[mef].avscint = 0.0;

    meflist[mef].frames = (struct lc_frame *) NULL;
    meflist[mef].nf = nf;

    if(ext == -99) {
      /* Move there */
      ffmahd(inf, mef+2, (int *) NULL, &status);
      if(tmplf != inf)
	ffmahd(tmplf, mef+2, (int *) NULL, &status);
      if(status) {
	fitsio_err(errstr, status, "ffmahd: %s: HDU %d", refname, mef+2);
	fatal(1, "%s", errstr);
      }
    }

    /* Read it in */
    if(read_lc(inf, &(meflist[mef]), errstr))
      fatal(1, "read_lc: HDU %d: %s", mef+2, errstr);

    /* Get disk buffer */
    if(buffer_alloc(&buf, meflist[mef].nstars, nf, NFLUX, errstr))
      fatal(1, "buffer_alloc: %s", errstr);

    /* Allocate buffer for frame info */
    meflist[mef].frames = (struct lc_frame *) malloc(nf * sizeof(struct lc_frame));
    if(!meflist[mef].frames)
      error(1, "malloc");

    /* Read frames */
    for(f = 0; f < nf; f++) {
      if(verbose && isatty(1))
	printf("\r Reading %*s (%*ld of %*ld)", maxflen, fnlist[f], fspc, f+1, fspc, nf);

      if(read_cat(fnlist[f], f, mef, &(meflist[mef]), &buf,
		  dointra, &(intralist[mef]), diffmode, satlev, errstr))
	fatal(1, "read_cat: %s: %s", fnlist[f], errstr);
    }

    if(verbose && isatty(1))
      printf("\n");

    /* Sort out averages */
    meflist[mef].avsigma = sqrtf(meflist[mef].avsigma / nf);
    meflist[mef].avapcor /= nf;
    meflist[mef].avscint /= nf;

    /* Change MJD to be relative to the first frame */
    for(f = 0; f < nf; f++)
      meflist[mef].frames[f].mjd -= meflist[mef].mjdref;

    if(verbose) {
      printf("  Number of objects:    %ld\n", meflist[mef].nstars);

      if(meflist[mef].satmag != -999.0)
	printf("  Saturation level:     %.1f\n", meflist[mef].satmag);
      else
	printf("  Saturation level:     undetermined\n");

      printf("  5-sigma limit:        %.1f\n",
	     meflist[mef].refflim);
    }

    if(meflist[mef].satmag != -999.0) {
      medbuf1[nmedsat] = meflist[mef].satmag;
      nmedsat++;
    }

    medbuf2[nmedlim] = meflist[mef].refflim;
    nmedlim++;

    nstartot += meflist[mef].nstars;

    /* Call into the main part of the program */
    if(lightcurves_append(&buf, &(meflist[mef]), errstr))
      fatal(1, "%s", errstr);

    /* Calculate average extinction */
    if(meflist[mef].degree >= 0) {
      meflist[mef].avextinc = 0.0;
      
      for(f = 0; f < nf; f++)
	meflist[mef].avextinc += powf(10.0, -0.4 * meflist[mef].frames[f].extinc);
      
      meflist[mef].avextinc /= nf;
    }
    else
      meflist[mef].avextinc = 1.0;

    /* Write out lightcurves for this MEF if requested */
    if(dooutput || doreplace) {
      if(verbose)
	printf(" Writing %s\n", doreplace ? tmpfile : outfile);

      /* Append a new HDU */
      ffcrhd(outf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffcrhd");
	fatal(1, "%s", errstr);
      }

      if(update_lc(tmplf, outf, &buf, &(meflist[mef]), outcls, wantoutcls, errstr))
	fatal(1, "write_lc: %s", errstr);

      /* Flush */
      ffflus(outf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffflus");
	fatal(1, "%s", errstr);
      }
    }

    sysulim = meflist[mef].sysulim;  /* kludge */
    sysllim = meflist[mef].sysllim;  /* kludge */
  }

#ifdef DEBUG
  cpgclos();
#endif

  /* Release disk buffer */
  buffer_close(&buf);

  /* Close reference */
  if(tmplf != inf)
    ffclos(tmplf, &status);
  ffclos(inf, &status);
  if(status) {
    fitsio_err(errstr, status, "ffclos");
    fatal(1, "%s", errstr);
  }

  /* Close output file */
  if(dooutput || doreplace) {
    ffclos(outf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffclos");
      fatal(1, "%s", errstr);
    }

    if(doreplace) {
      if(dooutput) {
	/* Overwrite original output with temporary file */
	rv = unlink(outfile);
	if(rv == -1)
	  error(1, "unlink: %s", outfile);
	
	rv = rename(tmpfile, outfile);
	if(rv == -1)
	  error(1, "rename: %s to %s", tmpfile, outfile);
      }
      else {
	/* Overwrite original input with temporary file */
	rv = unlink(refname);
	if(rv == -1)
	  error(1, "unlink: %s", refname);
	
	rv = rename(tmpfile, refname);
	if(rv == -1)
	  error(1, "rename: %s to %s", tmpfile, refname);
      }
    }
  }

  medsig(medbuf1, nmedsat, &medsat, (float *) NULL);
  medsig(medbuf2, nmedlim, &medlim, (float *) NULL);

  free((void *) medbuf1);
  medbuf1 = (float *) NULL;

  if(verbose)
    printf("\n"
	   "Total objects:           %ld\n"
	   "Median saturation level: %.1f\n"
	   "Median 5-sigma limit:    %.1f\n",
	   nstartot, medsat, medlim);

  /* Do diagnostic plots */
  if(!noplots) {
    if(do_plots(meflist, nmefs, medsat, medlim,
		sysulim < 0.0 ? medsat : sysulim,
		sysllim < 0.0 ? (sysulim < 0.0 ? medsat : sysulim)+USEMAG : sysllim,
		errstr))
      fatal(1, "do_plots: %s");
  }

  /* Free MEF info */
  for(mef = 0; mef < nmefs; mef++) {
    if(meflist[mef].stars)
      free((void *) meflist[mef].stars);
    if(meflist[mef].frames)
      free((void *) meflist[mef].frames);
    if(intralist[mef].map)
      free((void *) intralist[mef].map);
  }

  free((void *) meflist);
  meflist = (struct lc_mef *) NULL;
  free((void *) intralist);
  intralist = (struct intra *) NULL;

  /* Free file list */
  free((void *) fnlist);
  fnlist = (char **) NULL;
  
  return(0);
}

static int update_lc (fitsfile *reff, fitsfile *fits,
		      struct buffer_info *buf, struct lc_mef *mefinfo,
		      int outcls, int wantoutcls, 
		      char *errstr) {
  int status = 0, col, ncoluse, anynull;

  char *colnames[14] = { "medflux", "rms", "chisq", "nchisq", "pointer",
			 "hjd", "flux", "fluxerr", "xlc", "ylc", "airmass", "ha",
                         "weight", "flags" };
  int gcols[14];

  char kbuf[FLEN_KEYWORD];
  char cbuf[FLEN_COMMENT];

  long pt, star;

  long nmeasexist, nmeasout;

  struct lc_point *lcbuf = (struct lc_point *) NULL;

  double *epos = (double *) NULL;

  float *fluxbuf = (float *) NULL, *fluxerrbuf, *xlcbuf, *ylcbuf, *airbuf, *habuf, *wtbuf;
  double *hjdbuf = (double *) NULL;
  unsigned char *flagbuf = (unsigned char *) NULL;

  int ikey, nkeys, keylen;
  char card[FLEN_CARD], key[FLEN_KEYWORD];

  long satflag;
  unsigned char flags;

  int ncols, tcol;
  long previnstart, prevoutstart, previnpos, prevoutpos;
  int intype, outtype;
  long inrpt, inwidth, inbytes, outrpt, outwidth, outbytes;
  int found;

  struct table_section *copysect = (struct table_section *) NULL;
  long sect, ncopysect = 0;

  unsigned char *rawbuf = (unsigned char *) NULL;
  long rowsize;

  float *medlist, tmp;
  float medflux, sigflux, chisq;
  long nchisq, nmed;

  long starin, starout, nstarin, nstarout = 0, pointer;

  /* Get existing number of measurements */
  ffgkyj(reff, "NMEAS", &nmeasexist, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkyj: NMEAS");
    goto error;
  }

  /* Get number of rows and columns */
  ffgnrw(reff, &nstarin, &status);
  ffgncl(reff, &ncols, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgncl");
    goto error;
  }

  /* Get column numbers */
  ncoluse = sizeof(colnames) / sizeof(colnames[0]);
  
  for(col = 0; col < ncoluse; col++) {
    ffgcno(reff, CASEINSEN, colnames[col], &(gcols[col]), &status);
    if(status == COL_NOT_UNIQUE)
      status = 0;  /* ignore */
    else if(status) {
      fitsio_err(errstr, status, "ffgcno: %s", colnames[col]);
      goto error;
    }
  }

  /* Copy in header keywords */
  ffghsp(reff, &nkeys, (int *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffghsp");
    goto error;
  }

  for(ikey = 0; ikey < nkeys; ikey++) {
    ffgrec(reff, ikey+1, card, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgrec: card %d", ikey+1);
      goto error;
    }

    ffgknm(card, key, &keylen, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgknm: %s", card);
      goto error;
    }

    /* Change NAXIS2 to zero so we can manipulate the table empty */
    if(!strcmp(key, "NAXIS2")) {
      ffpkyj(fits, "NAXIS2", 0, "number of rows in table", &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: NAXIS2");
	goto error;
      }
    }
    else {
      /* Everything else gets copied */
      ffprec(fits, card, &status);
      if(status) {
	fitsio_err(errstr, status, "ffprec: card %d", ikey+1);
	goto error;
      }
    }
  }

  /* Allocate buffer for earth positions */
  epos = (double *) malloc(3 * mefinfo->nf * sizeof(double));
  if(!epos) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Update number of measurements and number of updates */
  nmeasout = nmeasexist+mefinfo->nf;

  ffukyj(fits, "NMEAS", nmeasout,
	 "Number of points in each lightcurve", &status);
  ffukyj(fits, "NUPDATE", mefinfo->nupdate+1,
	 "Number of times file has been appended to", &status);
  if(status) {
    fitsio_err(errstr, status, "ffkpy: frame info");
    goto error;
  }

  /* Add new keywords */
  for(pt = 0; pt < mefinfo->nf; pt++) {
    snprintf(kbuf, sizeof(kbuf), "TV%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Time value for datapoint %ld", nmeasexist+pt+1);
    ffpkyg(fits, kbuf, mefinfo->frames[pt].mjd, 7, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyg: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "OFF%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Frame offset for datapoint %ld", nmeasexist+pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].offset, 4, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "RMS%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Frame RMS for datapoint %ld", nmeasexist+pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].rms, 4, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "EXTC%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Extinction correction for datapoint %ld", nmeasexist+pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].extinc, 4, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "SEE%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Seeing for datapoint %ld", nmeasexist+pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].seeing, 3, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "FANG%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Field angle for datapoint %ld", nmeasexist+pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].fang, 6, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "IANG%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Field angle modulo pi for datapoint %ld", nmeasexist+pt+1);
    ffpkyj(fits, kbuf, mefinfo->frames[pt].iang, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    if(mefinfo->frames[pt].tamb != -999) {
      snprintf(kbuf, sizeof(kbuf), "TAMB%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "[degC] Ambient temp for datapoint %ld", nmeasexist+pt+1);
      ffpkyf(fits, kbuf, mefinfo->frames[pt].tamb, 1, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].humid != -999) {
      snprintf(kbuf, sizeof(kbuf), "HUM%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "[%%] Humidity for datapoint %ld", nmeasexist+pt+1);
      ffpkyf(fits, kbuf, mefinfo->frames[pt].humid, 1, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].press != -999) {
      snprintf(kbuf, sizeof(kbuf), "PRES%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "[hPa] Pressure for datapoint %ld", nmeasexist+pt+1);
      ffpkyf(fits, kbuf, mefinfo->frames[pt].press, 1, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].skytemp != -999) {
      snprintf(kbuf, sizeof(kbuf), "TSKY%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "[degC] Sky temp for datapoint %ld", nmeasexist+pt+1);
      ffpkyf(fits, kbuf, mefinfo->frames[pt].skytemp, 1, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].split_nexp >= 0) {
      snprintf(kbuf, sizeof(kbuf), "IEXP%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Exposure index for datapoint %ld", nmeasexist+pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].split_iexp, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }

      snprintf(kbuf, sizeof(kbuf), "NEXP%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Exposure count for datapoint %ld", nmeasexist+pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].split_nexp, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }
    }

    snprintf(kbuf, sizeof(kbuf), "IUPD%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Update number when datapoint %ld was added", nmeasexist+pt+1);
    ffpkyj(fits, kbuf, mefinfo->nupdate+1, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
      goto error;
    }

    /* Calculate Earth's heliocentric position at this MJD */
    getearth(mefinfo->mjdref + mefinfo->frames[pt].mjd, epos + 3*pt);
  }

  /* Expand vectors */
  for(col = 5; col < ncoluse; col++) {
    ffmvec(fits, gcols[col], nmeasout, &status);
    if(status) {
      fitsio_err(errstr, status, "ffmvec");
      goto error;
    }
  }

  /* Decide which bytes need copying */
  previnstart = 1;
  prevoutstart = 1;
  previnpos = 1;
  prevoutpos = 1;

  for(col = 1; col <= ncols; col++) {
    ffgtcl(reff, col, &intype, &inrpt, &inwidth, &status);
    ffgtcl(fits, col, &outtype, &outrpt, &outwidth, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgtcl: %d", col);
    }
    
    if(intype == TSTRING)
      inwidth = 1;
    if(outtype == TSTRING)
      outwidth = 1;
    
    inbytes = inrpt * inwidth;
    outbytes = outrpt * outwidth;
    
    found = 0;
    for(tcol = 0; tcol < ncoluse; tcol++)
      if(col == gcols[tcol])
	found = 1;

    if(found) {
      /* Skip this one */
      copysect = (struct table_section *) realloc(copysect,
						  (ncopysect+1) * sizeof(struct table_section));
      if(!copysect) {
	report_syserr(errstr, "realloc");
	goto error;
      }

      copysect[ncopysect].instart = previnstart;
      copysect[ncopysect].outstart = prevoutstart;
      copysect[ncopysect].ncopy = previnpos-previnstart;
      ncopysect++;

      previnstart = previnpos+inbytes;
      prevoutstart = prevoutpos+outbytes;
    }
    
    previnpos += inbytes;
    prevoutpos += outbytes;
  }

  copysect = (struct table_section *) realloc(copysect,
					      (ncopysect+1) * sizeof(struct table_section));
  if(!copysect) {
    report_syserr(errstr, "realloc");
    goto error;
  }
  
  copysect[ncopysect].instart = previnstart;
  copysect[ncopysect].outstart = prevoutstart;
  copysect[ncopysect].ncopy = previnpos-previnstart;
  ncopysect++;

  /* If we're filtering, figure out the correct size */
  if(wantoutcls) {
    nstarout = 0;

    for(starin = 0; starin < nstarin; starin++) {
      /* Read pointer to figure out original star number */
      ffgcvj(reff, gcols[4], starin + 1, 1, 1, 0, &star, &anynull, &status);
      if(status) {
	fitsio_err(errstr, status, "ffgcv");
	goto error;
      }
      
      star--;  /* 1-base to zero-base */

      if(mefinfo->stars[star].cls == outcls)
	nstarout++;
    }
  }
  else
    nstarout = nstarin;

  /* Expand the table to the correct size */
  ffukyj(fits, "NAXIS2", nstarout, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffpkyj: NAXIS2");
    goto error;
  }

  /* Get row size */
  ffgkyj(fits, "NAXIS1", &rowsize, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkyj: NAXIS1");
    goto error;
  }

  /* Allocate input buffer */
  lcbuf = (struct lc_point *) malloc(mefinfo->nf * sizeof(struct lc_point));
  if(!lcbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Allocate output buffers */
  fluxbuf = (float *) malloc(8 * nmeasout * sizeof(float));
  hjdbuf = (double *) malloc(nmeasout * sizeof(double));
  flagbuf = (unsigned char *) malloc(nmeasout * sizeof(unsigned char));
  rawbuf = (unsigned char *) malloc(rowsize * sizeof(unsigned char));
  if(!fluxbuf || !hjdbuf || !flagbuf || !rawbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  fluxerrbuf = fluxbuf + nmeasout;
  xlcbuf = fluxbuf + 2 * nmeasout;
  ylcbuf = fluxbuf + 3 * nmeasout;
  airbuf = fluxbuf + 4 * nmeasout;
  habuf = fluxbuf + 5 * nmeasout;
  wtbuf = fluxbuf + 6 * nmeasout;

  medlist = fluxbuf + 7 * nmeasout;

  /* Loop through all stars */
  starout = 0;

  for(starin = 0; starin < nstarin; starin++) {
    /* Read pointer to figure out original star number */
    ffgcvj(reff, gcols[4], starin + 1, 1, 1, 0, &pointer, &anynull, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgcv");
      goto error;
    }

    star = pointer-1;

    /* Skip the ones without the correct class if we're doing that */
    if(wantoutcls && mefinfo->stars[star].cls != outcls)
      continue;  /* I'm too lazy to reindent the rest of the loop */

    /* Read existing lightcurve info */
    ffgcvd(reff, gcols[5], starin + 1, 1, nmeasexist, -999.0, hjdbuf, &anynull, &status);
    ffgcve(reff, gcols[6], starin + 1, 1, nmeasexist, -999.0, fluxbuf, &anynull, &status);
    ffgcve(reff, gcols[7], starin + 1, 1, nmeasexist, -999.0, fluxerrbuf, &anynull, &status);
    ffgcve(reff, gcols[8], starin + 1, 1, nmeasexist, -999.0, xlcbuf, &anynull, &status);
    ffgcve(reff, gcols[9], starin + 1, 1, nmeasexist, -999.0, ylcbuf, &anynull, &status);
    ffgcve(reff, gcols[10], starin + 1, 1, nmeasexist, -999.0, airbuf, &anynull, &status);
    ffgcve(reff, gcols[11], starin + 1, 1, nmeasexist, -999.0, habuf, &anynull, &status);
    ffgcve(reff, gcols[12], starin + 1, 1, nmeasexist, -999.0, wtbuf, &anynull, &status);
    ffgcvb(reff, gcols[13], starin + 1, 1, nmeasexist, 0, flagbuf, &anynull, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgcv");
      goto error;
    }

    /* Get lightcurve */
    if(buffer_fetch_object(buf, lcbuf, 0, mefinfo->nf, star, 0, errstr))
      goto error;

    /* Fill in buffer */
    satflag = 0;
    for(pt = 0; pt < mefinfo->nf; pt++) {
      flags = 0;

      if(lcbuf[pt].flux != 0.0) {
	fluxbuf[nmeasexist+pt] = mefinfo->zp - lcbuf[pt].flux;

	if(abs(lcbuf[pt].flux > 20)) {
	  printf("Warning: daft-looking flux for star %ld point %ld: %.2g\n",
		 star+1, pt+1, lcbuf[pt].flux);
	}

	/* Unset the all saturated flag if not saturated */
	if(lcbuf[pt].conf)
	  flags |= FLAG_CONF;
	if(lcbuf[pt].satur) {
	  flags |= FLAG_SATUR;
	  satflag++;
	}

	if(lcbuf[pt].fluxerrcom > 0.0)
	  fluxerrbuf[nmeasexist+pt] = lcbuf[pt].fluxerrcom;
	else
	  fluxerrbuf[nmeasexist+pt] = -999.0;
      }
      else {
	fluxbuf[nmeasexist+pt] = -999.0;
	fluxerrbuf[nmeasexist+pt] = -999.0;
	flags |= FLAG_NODP;
      }

      xlcbuf[nmeasexist+pt] = lcbuf[pt].x;
      ylcbuf[nmeasexist+pt] = lcbuf[pt].y;
      airbuf[nmeasexist+pt] = lcbuf[pt].airmass;
      habuf[nmeasexist+pt] = lcbuf[pt].ha;
      wtbuf[nmeasexist+pt] = lcbuf[pt].wt;

      flagbuf[nmeasexist+pt] = flags;

      /* Calculate HJD (as UTC) */
      hjdbuf[nmeasexist+pt] = mefinfo->mjdref + mefinfo->frames[pt].mjd +
	                      hjdcorr(epos + 3*pt,
				      mefinfo->stars[star].ra,
				      mefinfo->stars[star].dec);
    }

    /* Raw-copy over the bits that didn't change */
    for(sect = 0; sect < ncopysect; sect++) {
      if(copysect[sect].ncopy <= 0)
	continue;

      ffgtbb(reff, starin+1, copysect[sect].instart, copysect[sect].ncopy, rawbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffgtbb");
	goto error;
      }

      ffptbb(fits, starout+1, copysect[sect].outstart, copysect[sect].ncopy, rawbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffptbb");
	goto error;
      }
    }

    /* Recalculate median flux */
    nmed = 0;
    
    for(pt = 0; pt < nmeasout; pt++) {
      if(fluxbuf[pt] != -999.0 && fluxerrbuf[pt] != -999.0) {
	medlist[nmed] = fluxbuf[pt];
	nmed++;
      }
    }
    
    if(nmed > 0) {
      medsig(medlist, nmed, &medflux, &sigflux);
    
      if(nmed == 1)
	sigflux = -999.0;

      /* Calculate chi-squared */
      chisq = 0.0;
      nchisq = 0;
      
      for(pt = 0; pt < nmeasout; pt++) {
	if(fluxbuf[pt] != -999.0 && fluxerrbuf[pt] != -999.0) {
	  tmp = fluxbuf[pt] - medflux;
	  
	  chisq += tmp*tmp / (fluxerrbuf[pt] * fluxerrbuf[pt]);
	  nchisq++;
	}
      }
    }
    else {
      medflux = -999.0;
      sigflux = -999.0;
      chisq = -999.0;
      nchisq = 0;
    }

    /* Write in modified data */
    ffpcne(fits, gcols[0], starout+1, 1, 1, &medflux, -999.0, &status);
    ffpcne(fits, gcols[1], starout+1, 1, 1, &sigflux, -999.0, &status);
    ffpcne(fits, gcols[2], starout+1, 1, 1, &chisq, -999.0, &status);
    ffpcnj(fits, gcols[3], starout+1, 1, 1, &nchisq, -999, &status);
    ffpclj(fits, gcols[4], starout+1, 1, 1, &pointer, &status);
    ffpcnd(fits, gcols[5], starout+1, 1, nmeasout, hjdbuf, -999.0, &status);
    ffpcne(fits, gcols[6], starout+1, 1, nmeasout, fluxbuf, -999.0, &status);
    ffpcne(fits, gcols[7], starout+1, 1, nmeasout, fluxerrbuf, -999.0, &status);
    ffpcne(fits, gcols[8], starout+1, 1, nmeasout, xlcbuf, -999.0, &status);
    ffpcne(fits, gcols[9], starout+1, 1, nmeasout, ylcbuf, -999.0, &status);
    ffpcne(fits, gcols[10], starout+1, 1, nmeasout, airbuf, -999.0, &status);
    ffpcne(fits, gcols[11], starout+1, 1, nmeasout, habuf, -999.0, &status);
    ffpcne(fits, gcols[12], starout+1, 1, nmeasout, wtbuf, -999.0, &status);
    ffpclb(fits, gcols[13], starout+1, 1, nmeasout, flagbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpcl");
      goto error;
    }

    starout++;
  }

  free((void *) lcbuf);
  lcbuf = (struct lc_point *) NULL;
  free((void *) epos);
  epos = (double *) NULL;
  free((void *) copysect);
  copysect = (struct table_section *) NULL;
  free((void *) fluxbuf);
  fluxbuf = (float *) NULL;
  free((void *) hjdbuf);
  hjdbuf = (double *) NULL;
  free((void *) flagbuf);
  flagbuf = (unsigned char *) NULL;
  free((void *) rawbuf);
  rawbuf = (unsigned char *) NULL;

  return(0);

 error:
  if(lcbuf)
    free((void *) lcbuf);
  if(epos)
    free((void *) epos);
  if(copysect)
    free((void *) copysect);
  if(fluxbuf)
    free((void *) fluxbuf);
  if(hjdbuf)
    free((void *) hjdbuf);
  if(flagbuf)
    free((void *) flagbuf);
  if(rawbuf)
    free((void *) rawbuf);

  return(1);
}

