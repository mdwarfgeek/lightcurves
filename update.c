#include <sys/types.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#ifdef DEBUG
#include <cpgplot.h>
#endif

#include "lightcurves.h"

#ifdef HJD
#include "hjd.h"
#endif

#include "cvtunit.h"
#include "fitsutil.h"
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
	  "         -V file   Use 'instrument version' table from 'file' (MEarth only).\n\n"
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

  int f, nf = 0;

  fitsfile *inf, *outf, *tmplf;
  int status = 0, ext, mef, nmefs;

  char outfile[FLEN_FILENAME-1], fnbuf[FLEN_FILENAME];
  char tmpbasebuf[FLEN_FILENAME-1], *tmpbase;

#ifdef _WIN32
  char tmpfile[MAX_PATH];
#else
  char tmpfile[FLEN_FILENAME-1];
#endif

  int dooutput = 0;
  int doreplace = 0;
  int outcls = 0;
  int wantoutcls = 0;
  int fd = -1, rv;

  char intrafile[FLEN_FILENAME];
  int dointra = 0;

  char instversfile[FLEN_FILENAME];
  int doinstvers = 0;
  struct instvers *instverslist = (struct instvers *) NULL;
  int ninstvers = 0;

  struct dtai_table dtab, *dtptr = NULL;
  struct iers_table itab, *itptr = NULL;
  struct jpleph_table jtab, ttab, *jtptr = NULL, *ttptr = NULL;

  int len, maxflen, fspc;
  float *medbuf1 = (float *) NULL, *medbuf2, medsat, medlim;
  long nmedsat, nmedlim, nstartot;

  int diffmode = 0;
  int noplots = 0;

  float satlev = -1.0;
#ifdef PLOTS
  float sysulim = -1.0, sysllim = -1.0;
#endif

  /* Set the program name for error reporting */
  if(argv[0])
    pn = basename(argv[0]);
  
  if(!pn)
    pn = "update";

  setprogname(pn);

  avzero = argv[0];

  /* Extract command-line arguments */
  while((c = getopt(argc, argv, "c:i:o:s:upqvV:")) != -1)
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
    case 'V':
      strncpy(instversfile, optarg, sizeof(instversfile)-1);
      instversfile[sizeof(instversfile)-1] = '\0';
      doinstvers = 1;
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

  fnlist = read_file_list(argc-1, argv+1, &nf, errstr);
  if(!fnlist)
    fatal(1, "%s", errstr);

  /* Setup Earth orientation data and JPL ephemerides */
  rv = dtai_read(&dtab, (char *) NULL);
  if(rv == -2) {
    printf("Could not find leap second file, continuing without\n");
    dtptr = NULL;
  }
  else if(rv)
    fatal(1, "dtai_open: error %d", rv);
  else
    dtptr = &dtab;

  rv = iers_open(&itab, &dtab, (char *) NULL);
  if(rv == -2) {
    printf("Could not find UT1-UTC file, continuing without\n");
    itptr = NULL;
  }
  else if(rv)
    fatal(1, "iers_open: %d", rv);
  else
    itptr = &itab;

  rv = jpleph_open(&jtab, 0, (char *) NULL);
  if(rv == -2) {
    printf("Could not find JPL ephemeris, continuing without\n");
    jtptr = NULL;
    ttptr = NULL;
  }
  else if(rv)
    fatal(1, "jpleph_open: %d", rv);
  else {
    jtptr = &jtab;

    if(!jtab.has_time) {
      rv = jpleph_open(&ttab, 1, (char *) NULL);
      if(rv == -2) {
        printf("Time ephemeris problem, continuing without\n");
        ttptr = NULL;
      }
      else if(rv)
        fatal(1, "jpleph_open: %d", rv);
      else
        ttptr = &ttab;
    }
  }

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
	fitsio_err(errstr, status, "ffopen: %s", outfile);
	fatal(1, "%s", errstr);
      }
    }
  }

  if(doreplace) {
    /* First figure out where to create it - must be same device as output
     * for rename().  Copy first because dirname may modify its argument.
     */
    strncpy(tmpbasebuf, outfile, sizeof(tmpbasebuf)-1);
    tmpbasebuf[sizeof(tmpbasebuf)-1] = '\0';

    tmpbase = dirname(tmpbasebuf);

    /* Create temporary file for in-place edit */
#ifdef _WIN32
    /* Obtain temporary file name */
    if(GetTempFileName(tmpbase, 
                       progname,
                       0,
                       tmpfile) == 0)
      fatal(1, "GetTempFileName failed");
#else
    snprintf(tmpfile, sizeof(tmpfile), "%s/%s_XXXXXX", tmpbase, progname);

    fd = mkstemp(tmpfile);
    if(fd == -1)
      error(1, "mkstemp: %s", tmpfile);
#endif

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

  /* Load "instrument version" table if requested */
  if(doinstvers) {
    if(read_instvers(instversfile, &instverslist, &ninstvers, errstr))
      fatal(1, "read_instvers: %s", errstr);
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
    meflist[mef].nrows = 0;

    meflist[mef].warned = 0;

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
    if(buffer_alloc(&buf, meflist[mef].nstars, nf, errstr))
      fatal(1, "buffer_alloc: %s", errstr);

    /* Allocate buffer for frame info */
    meflist[mef].frames = (struct lc_frame *) malloc(nf * sizeof(struct lc_frame));
    if(!meflist[mef].frames)
      error(1, "malloc");

    /* Read frames */
    for(f = 0; f < nf; f++) {
      if(verbose && isatty(1))
	printf("\r Reading %*s (%*d of %*d)",
               maxflen, fnlist[f], fspc, f+1, fspc, nf);

      if(read_cat(fnlist[f], f, mef, &(meflist[mef]), &buf,
		  dtptr, itptr, jtptr, ttptr,
		  dointra, &(intralist[mef]),
		  doinstvers, instverslist, ninstvers,
		  diffmode, satlev, errstr))
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
    if(lightcurves_append(&buf, &(meflist[mef]),
			  meflist[mef].nstars != meflist[mef].nrows, errstr))
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

#ifdef PLOTS
    sysulim = meflist[mef].sysulim;  /* kludge */
    sysllim = meflist[mef].sysllim;  /* kludge */
#endif
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
	rv = rename(tmpfile, outfile);
	if(rv == -1)
	  error(1, "rename: %s to %s", tmpfile, outfile);
      }
      else {
	/* Overwrite original input with temporary file */
	rv = rename(tmpfile, refname);
	if(rv == -1)
	  error(1, "rename: %s to %s", tmpfile, refname);
      }
    }
  }

  fmedsig(medbuf1, nmedsat, &medsat, (float *) NULL);
  fmedsig(medbuf2, nmedlim, &medlim, (float *) NULL);

  free((void *) medbuf1);
  medbuf1 = (float *) NULL;

  if(verbose)
    printf("\n"
	   "Total objects:           %ld\n"
	   "Median saturation level: %.1f\n"
	   "Median 5-sigma limit:    %.1f\n",
	   nstartot, medsat, medlim);

#ifdef PLOTS
  /* Do diagnostic plots */
  if(!noplots) {
    if(do_plots(meflist, nmefs, medsat, medlim,
		sysulim < 0.0 ? medsat : sysulim,
		sysllim < 0.0 ? (sysulim < 0.0 ? medsat : sysulim)+USEMAG : sysllim,
		errstr))
      fatal(1, "do_plots: %s");
  }
#endif

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
  int status = 0, icol, ntmpl, ncoluse, anynull;

  char *colnames[] = { "pointer",
		       "bjd",
#ifdef HJD
                       "hjd",
#endif
                       "flux", "fluxerr", "xlc", "ylc", "airmass", "ha",
		       "weight", "sky", "peak", "flags",
		       "medflux", "rms", "chisq", "nchisq" };
  unsigned char tpap[] = { 0,
			   0,
#ifdef HJD
			   0,
#endif
                           1, 1, 0, 0, 0, 0,
			   1, 0, 0, 0,
			   1, 1, 0, 0 };
  unsigned char tvec[] = { 0,
			   1,
#ifdef HJD
			   1,
#endif
                           1, 1, 1, 1, 1, 1,
			   1, 1, 1, 1,
			   0, 0, 0, 0 };

  char cn[FLEN_VALUE+1];
  int *gcols = (int *) NULL;

  char kbuf[FLEN_KEYWORD];
  char cbuf[FLEN_COMMENT];

  long pt, star, soff, istar;

  int ap, ap1, ap2, napcol;

  long nmeasexist, nupdate, nmeasout;

  struct lc_point *lcbuf = (struct lc_point *) NULL;

#ifdef HJD
  double *epos = (double *) NULL;
#endif

  float *medfluxbuf = (float *) NULL, *rmsbuf;
  float *fluxbuf = (float *) NULL, *fluxerrbuf, *wtbuf;
  float *airbuf = (float *) NULL, *habuf, *locskybuf, *peakbuf;
  double *bjdbuf = (double *) NULL, *xlcbuf, *ylcbuf;
#ifdef HJD
  double *hjdbuf;
#endif
  unsigned char *flagbuf = (unsigned char *) NULL;

  int ikey, nkeys, keylen;
  char card[FLEN_CARD], key[FLEN_KEYWORD];

  long satflag;
  unsigned char flags;

  int ncolsfile, tcol;
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
  int allast;

  /* Get existing number of measurements and number of updates */
  ffgkyj(reff, "NMEAS", &nmeasexist, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkyj: NMEAS");
    goto error;
  }

  ffgkyj(reff, "NUPDATE", &nupdate, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    nupdate = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkyj: NUPDATE");
    goto error;
  }

  /* Get number of rows and columns */
  ffgnrw(reff, &nstarin, &status);
  ffgncl(reff, &ncolsfile, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgncl");
    goto error;
  }

  /* Figure out apertures */
  if(mefinfo->aperture) {
    ap1 = mefinfo->aperture-1;
    ap2 = ap1;
    napcol = 1;
  }
  else {
    if(mefinfo->apselmode & APSEL_ALL) {
      ap1 = 0;
      ap2 = NFLUX;
      napcol = ap2-ap1 + 1;
    }
    else {
      ap1 = 0;
      ap2 = 0;
      napcol = 1;
    }
  }

  /* Decide how many columns */
  ntmpl = sizeof(colnames) / sizeof(colnames[0]);

  ncoluse = 0;
  for(tcol = 0; tcol < ntmpl; tcol++) {
    ncoluse++;

    if(tpap[tcol])
      ncoluse += ap2-ap1;
  }

  /* Allocate arrays */
  gcols = (int *) malloc(ncoluse * sizeof(int));
  if(!gcols) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Get column numbers */
  icol = 0;
  for(tcol = 0; tcol < ntmpl; tcol++) {
    ffgcno(reff, CASEINSEN, colnames[tcol], gcols+icol, &status);
    if(status == COL_NOT_UNIQUE)
      status = 0;  /* ignore */
    else if(status) {
      fitsio_err(errstr, status, "ffgcno: %s", colnames[tcol]);
      goto error;
    }

    icol++;

    if(tpap[tcol])
      for(ap = ap1; ap < ap2; ap++) {
	snprintf(cn, sizeof(cn), "%s%d", colnames[tcol], ap+1);
	ffgcno(reff, CASEINSEN, cn, gcols+icol, &status);
	if(status == COL_NOT_UNIQUE)
	  status = 0;  /* ignore */
	else if(status) {
	  fitsio_err(errstr, status, "ffgcno: %s", cn);
	  goto error;
	}

	icol++;
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

#ifdef HJD
  /* Allocate buffer for earth positions */
  epos = (double *) malloc(3 * mefinfo->nf * sizeof(double));
  if(!epos) {
    report_syserr(errstr, "malloc");
    goto error;
  }
#endif

  /* Update number of measurements and number of updates */
  nmeasout = nmeasexist+mefinfo->nf;

  ffukyj(fits, "NMEAS", nmeasout,
	 "Number of points in each lightcurve", &status);
  ffukyj(fits, "NUPDATE", nupdate+1,
	 "Number of times file has been appended to", &status);
  if(status) {
    fitsio_err(errstr, status, "ffkpy: frame info");
    goto error;
  }

  /* Read astrom flag */
  ffgkyl(fits, "ASTONLY", &allast, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    allast = -1;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkyl: ASTONLY");
    goto error;
  }

  /* Add new keywords */
  for(pt = 0; pt < mefinfo->nf; pt++) {
    snprintf(kbuf, sizeof(kbuf), "TV%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Time value for datapoint %ld", nmeasexist+pt+1);
    ffpkyg(fits, kbuf, mefinfo->frames[pt].mjd, 14, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyg: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "TEXP%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Exposure time for datapoint %ld", nmeasexist+pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].exptime, 3, cbuf, &status);
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

    snprintf(kbuf, sizeof(kbuf), "ELL%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Ellipticity for datapoint %ld", nmeasexist+pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].ellipt, 3, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "SKY%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Sky level for datapoint %ld", nmeasexist+pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].skylev, 2, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "NOIS%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Sky noise for datapoint %ld", nmeasexist+pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].skynoise, 2, cbuf, &status);
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

    if(mefinfo->frames[pt].rtstat >= 0) {
      snprintf(kbuf, sizeof(kbuf), "RTST%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Trigger status for datapoint %ld", nmeasexist+pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].rtstat, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].schpri >= 0) {
      snprintf(kbuf, sizeof(kbuf), "SPRI%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Scheduling priority for datapoint %ld", nmeasexist+pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].schpri, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].schcad >= 0) {
      snprintf(kbuf, sizeof(kbuf), "SCAD%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Scheduling cadence for datapoint %ld", nmeasexist+pt+1);
      ffpkye(fits, kbuf, mefinfo->frames[pt].schcad, 6, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkye: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].schtype[0]) {
      snprintf(kbuf, sizeof(kbuf), "STYP%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Observation type code for datapoint %ld", nmeasexist+pt+1);
      ffpkys(fits, kbuf, mefinfo->frames[pt].schtype, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkys: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].cadencenum != -999) {
      snprintf(kbuf, sizeof(kbuf), "CADN%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Cadence number for datapoint %ld", nmeasexist+pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].cadencenum, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }
    }

    snprintf(kbuf, sizeof(kbuf), "ISEG%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Segment number for datapoint %ld", nmeasexist+pt+1);
    ffpkyj(fits, kbuf, mefinfo->frames[pt].iseg+1, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
      goto error;
    }

    if(mefinfo->frames[pt].instvers) {
      snprintf(kbuf, sizeof(kbuf), "IVER%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Instrument version for datapoint %ld", nmeasexist+pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].instvers->iver, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }

      snprintf(kbuf, sizeof(kbuf), "IDAT%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Inst last change date for datapoint %ld", nmeasexist+pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].instvers->date, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }
    }

    if(allast >= 0 && strcmp(mefinfo->frames[pt].schtype, "a"))
      allast = 0;

    snprintf(kbuf, sizeof(kbuf), "IUPD%ld", nmeasexist+pt+1);
    snprintf(cbuf, sizeof(cbuf), "Update number when datapoint %ld was added", nmeasexist+pt+1);
    ffpkyj(fits, kbuf, nupdate+1, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
      goto error;
    }

    /* Suppress output of frame transformation if we only had the
       photometric comparison stars.  It really needs everything. */
    if(mefinfo->nstars == mefinfo->nrows) {
      snprintf(kbuf, sizeof(kbuf), "LXX%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Linear transformation to ref. for frame %ld", nmeasexist+pt+1);
      ffpkyd(fits, kbuf, mefinfo->frames[pt].tr[0], 12, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyd: %s", kbuf);
	goto error;
      }
      
      snprintf(kbuf, sizeof(kbuf), "LXY%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Linear transformation to ref. for frame %ld", nmeasexist+pt+1);
      ffpkyd(fits, kbuf, mefinfo->frames[pt].tr[1], 12, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyd: %s", kbuf);
	goto error;
      }
      
      snprintf(kbuf, sizeof(kbuf), "LXD%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Linear transformation to ref. for frame %ld", nmeasexist+pt+1);
      ffpkyd(fits, kbuf, mefinfo->frames[pt].tr[2], 12, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyd: %s", kbuf);
	goto error;
      }
      
      snprintf(kbuf, sizeof(kbuf), "LYY%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Linear transformation to ref. for frame %ld", nmeasexist+pt+1);
      ffpkyd(fits, kbuf, mefinfo->frames[pt].tr[3], 12, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyd: %s", kbuf);
	goto error;
      }
      
      snprintf(kbuf, sizeof(kbuf), "LYX%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Linear transformation to ref. for frame %ld", nmeasexist+pt+1);
      ffpkyd(fits, kbuf, mefinfo->frames[pt].tr[4], 12, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyd: %s", kbuf);
	goto error;
      }
      
      snprintf(kbuf, sizeof(kbuf), "LYD%ld", nmeasexist+pt+1);
      snprintf(cbuf, sizeof(cbuf), "Linear transformation to ref. for frame %ld", nmeasexist+pt+1);
      ffpkyd(fits, kbuf, mefinfo->frames[pt].tr[5], 12, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }
    }

#ifdef HJD
    /* Calculate Earth's heliocentric position at this MJD */
    getearth(mefinfo->mjdref + mefinfo->frames[pt].mjd, epos + 3*pt);
#endif
  }

  if(allast >= 0) {
    ffukyl(fits, "ASTONLY", allast, "Are all observations for astrometry?", &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyl: ASTONLY");
      goto error;
    }
  }

  /* Expand vectors */
  icol = 0;
  for(tcol = 0; tcol < ntmpl; tcol++) {
    if(tvec[tcol]) {
      ffmvec(fits, gcols[icol], nmeasout, &status);
      icol++;

      if(tpap[tcol])
	for(ap = ap1; ap < ap2; ap++) {	
	  ffmvec(fits, gcols[icol], nmeasout, &status);
	  icol++;
	}
    }
    else {
      icol++;
      
      if(tpap[tcol])
	icol += ap2-ap1;
    }
  }

  if(status) {
    fitsio_err(errstr, status, "ffmvec");
    goto error;
  }

  /* Decide which bytes need copying */
  previnstart = 1;
  prevoutstart = 1;
  previnpos = 1;
  prevoutpos = 1;

  for(icol = 1; icol <= ncolsfile; icol++) {
    ffgtcl(reff, icol, &intype, &inrpt, &inwidth, &status);
    ffgtcl(fits, icol, &outtype, &outrpt, &outwidth, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgtcl: %d", icol);
    }
    
    if(intype == TSTRING)
      inwidth = 1;
    if(outtype == TSTRING)
      outwidth = 1;
    
    inbytes = inrpt * inwidth;
    outbytes = outrpt * outwidth;
    
    found = 0;
    for(tcol = 0; tcol < ncoluse; tcol++)
      if(icol == gcols[tcol])
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
      ffgcvj(reff, gcols[0], starin + 1, 1, 1, 0, &pointer, &anynull, &status);
      if(status) {
	fitsio_err(errstr, status, "ffgcv");
	goto error;
      }
      
      star = -1;
      for(istar = 0; istar < mefinfo->nstars; istar++)
	if(mefinfo->stars[istar].ptr == pointer) {
	  star = istar;
	  break;
	}

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
  medfluxbuf = (float *) malloc(2 * napcol * sizeof(float));
  fluxbuf = (float *) malloc(3 * nmeasout * napcol * sizeof(float));
  airbuf = (float *) malloc(5 * nmeasout * sizeof(float));
#ifdef HJD
  bjdbuf = (double *) malloc(4 * nmeasout * sizeof(double));
#else
  bjdbuf = (double *) malloc(3 * nmeasout * sizeof(double));
#endif
  flagbuf = (unsigned char *) malloc(nmeasout * sizeof(unsigned char));
  rawbuf = (unsigned char *) malloc(rowsize * sizeof(unsigned char));
  if(!medfluxbuf || !fluxbuf || !airbuf || !bjdbuf || !flagbuf || !rawbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  rmsbuf = medfluxbuf + napcol;

  fluxerrbuf = fluxbuf + napcol * nmeasout;
  wtbuf = fluxbuf + 2 * napcol * nmeasout;

  habuf = airbuf + nmeasout;
  locskybuf = airbuf + 2 * nmeasout;
  peakbuf = airbuf + 3 * nmeasout;
  medlist = airbuf + 4 * nmeasout;

  xlcbuf = bjdbuf + nmeasout;
  ylcbuf = bjdbuf + 2 * nmeasout;
#ifdef HJD
  hjdbuf = bjdbuf + 3 * nmeasout;
#endif

  /* Loop through all stars */
  starout = 0;

  for(starin = 0; starin < nstarin; starin++) {
    /* Read pointer to figure out original star number */
    ffgcvj(reff, gcols[0], starin + 1, 1, 1, 0, &pointer, &anynull, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgcv");
      goto error;
    }

    star = -1;
    for(istar = 0; istar < mefinfo->nstars; istar++)
      if(mefinfo->stars[istar].ptr == pointer) {
	star = istar;
	break;
      }

    /* Skip the ones without the correct class if we're doing that */
    if(wantoutcls && mefinfo->stars[star].cls != outcls)
      continue;  /* I'm too lazy to reindent the rest of the loop */

    /* Read existing lightcurve info */
    icol = 1;

    ffgcvd(reff, gcols[icol++], starin + 1, 1, nmeasexist, -999.0, bjdbuf, &anynull, &status);
#ifdef HJD
    ffgcvd(reff, gcols[icol++], starin + 1, 1, nmeasexist, -999.0, hjdbuf, &anynull, &status);
#endif

    for(ap = 0; ap < napcol; ap++)
      ffgcve(reff, gcols[icol++], starin + 1, 1, nmeasexist, -999.0, fluxbuf+ap*nmeasout, &anynull, &status);

    for(ap = 0; ap < napcol; ap++)
      ffgcve(reff, gcols[icol++], starin + 1, 1, nmeasexist, -999.0, fluxerrbuf+ap*nmeasout, &anynull, &status);

    ffgcvd(reff, gcols[icol++], starin + 1, 1, nmeasexist, -999.0, xlcbuf, &anynull, &status);
    ffgcvd(reff, gcols[icol++], starin + 1, 1, nmeasexist, -999.0, ylcbuf, &anynull, &status);
    ffgcve(reff, gcols[icol++], starin + 1, 1, nmeasexist, -999.0, airbuf, &anynull, &status);
    ffgcve(reff, gcols[icol++], starin + 1, 1, nmeasexist, -999.0, habuf, &anynull, &status);

    for(ap = 0; ap < napcol; ap++)
      ffgcve(reff, gcols[icol++], starin + 1, 1, nmeasexist, -999.0, wtbuf+ap*nmeasout, &anynull, &status);

    ffgcve(reff, gcols[icol++], starin + 1, 1, nmeasexist, -999.0, locskybuf, &anynull, &status);
    ffgcve(reff, gcols[icol++], starin + 1, 1, nmeasexist, -999.0, peakbuf, &anynull, &status);
    ffgcvb(reff, gcols[icol++], starin + 1, 1, nmeasexist, 0, flagbuf, &anynull, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgcv");
      goto error;
    }

    /* Get lightcurve */
    if(buffer_fetch_object(buf, lcbuf, 0, mefinfo->nf, star, errstr))
      goto error;

    /* Fill in buffer */
    satflag = 0;
    for(pt = 0; pt < mefinfo->nf; pt++) {
      flags = 0;

      if(lcbuf[pt].aper[mefinfo->stars[star].iap].flux != 0.0) {
	fluxbuf[nmeasexist+pt] = mefinfo->zp - lcbuf[pt].aper[mefinfo->stars[star].iap].flux;

	if(abs(lcbuf[pt].aper[mefinfo->stars[star].iap].flux > 20)) {
	  printf("Warning: daft-looking flux for star %ld point %ld: %.2g\n",
		 star+1, pt+1, lcbuf[pt].aper[mefinfo->stars[star].iap].flux);
	}

	/* Unset the all saturated flag if not saturated */
	if(lcbuf[pt].conf)
	  flags |= FLAG_CONF;
	if(lcbuf[pt].satur) {
	  flags |= FLAG_SATUR;
	  satflag++;
	}

	if(lcbuf[pt].aper[mefinfo->stars[star].iap].fluxvarcom > 0.0)
	  fluxerrbuf[nmeasexist+pt] = sqrtf(lcbuf[pt].aper[mefinfo->stars[star].iap].fluxvarcom);
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
      wtbuf[nmeasexist+pt] = lcbuf[pt].aper[mefinfo->stars[star].iap].wt;
      locskybuf[nmeasexist+pt] = lcbuf[pt].sky;
      peakbuf[nmeasexist+pt] = lcbuf[pt].peak;

      flagbuf[nmeasexist+pt] = flags;

#ifdef HJD
      /* Calculate HJD (as UTC) */
      hjdbuf[nmeasexist+pt] = mefinfo->mjdref + mefinfo->frames[pt].mjd +
	                      hjdcorr(epos + 3*pt,
				      mefinfo->stars[star].ra,
				      mefinfo->stars[star].dec);
#endif
      bjdbuf[nmeasexist+pt] = lcbuf[pt].bjd;
    }

    /* Do other apertures */
    for(ap = ap1; ap < ap2; ap++) {
      soff = (ap-ap1+1)*nmeasout + nmeasexist;
      
      /* Fill in buffer */
      for(pt = 0; pt < mefinfo->nf; pt++) {
	if(lcbuf[pt].aper[ap].flux != 0.0) {
	  fluxbuf[soff+pt] = mefinfo->zp - lcbuf[pt].aper[ap].flux;
	  if(lcbuf[pt].aper[ap].fluxvarcom > 0.0)
	    fluxerrbuf[soff+pt] = sqrtf(lcbuf[pt].aper[ap].fluxvarcom);
	  else
	    fluxerrbuf[soff+pt] = -999.0;
	}
        else {
          fluxbuf[soff+pt] = -999.0;
          fluxerrbuf[soff+pt] = -999.0;
        }

	wtbuf[soff+pt] = lcbuf[pt].aper[ap].wt;
      }
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
    for(ap = 0; ap < napcol; ap++) {
      soff = ap*nmeasout;

      nmed = 0;
      
      for(pt = 0; pt < nmeasout; pt++) {
	if(fluxbuf[soff+pt] != -999.0 && fluxerrbuf[soff+pt] != -999.0) {
	  medlist[nmed] = fluxbuf[soff+pt];
	  nmed++;
	}
      }
      
      if(nmed > 0) {
	fmedsig(medlist, nmed, &medflux, &sigflux);
	
	if(nmed == 1)
	  sigflux = -999.0;

	medfluxbuf[ap] = medflux;
	rmsbuf[ap] = sigflux;
	
	if(ap == 0) {
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
      }
      else {
	medfluxbuf[ap] = -999.0;
	rmsbuf[ap] = -999.0;

	if(ap == 0) {
	  chisq = -999.0;
	  nchisq = 0;
	}
      }
    }

    /* Write in modified data */
    icol = 0;

    ffpclj(fits, gcols[icol++], starout+1, 1, 1, &pointer, &status);
    ffpcnd(fits, gcols[icol++], starout+1, 1, nmeasout, bjdbuf, -999.0, &status);
#ifdef HJD
    ffpcnd(fits, gcols[icol++], starout+1, 1, nmeasout, hjdbuf, -999.0, &status);
#endif

    for(ap = 0; ap < napcol; ap++)
      ffpcne(fits, gcols[icol++], starout+1, 1, nmeasout, fluxbuf+ap*nmeasout, -999.0, &status);

    for(ap = 0; ap < napcol; ap++)
      ffpcne(fits, gcols[icol++], starout+1, 1, nmeasout, fluxerrbuf+ap*nmeasout, -999.0, &status);

    ffpcnd(fits, gcols[icol++], starout+1, 1, nmeasout, xlcbuf, -999.0, &status);
    ffpcnd(fits, gcols[icol++], starout+1, 1, nmeasout, ylcbuf, -999.0, &status);
    ffpcne(fits, gcols[icol++], starout+1, 1, nmeasout, airbuf, -999.0, &status);
    ffpcne(fits, gcols[icol++], starout+1, 1, nmeasout, habuf, -999.0, &status);

    for(ap = 0; ap < napcol; ap++)
      ffpcne(fits, gcols[icol++], starout+1, 1, nmeasout, wtbuf+ap*nmeasout, -999.0, &status);

    ffpcne(fits, gcols[icol++], starout+1, 1, nmeasout, locskybuf, -999.0, &status);
    ffpcne(fits, gcols[icol++], starout+1, 1, nmeasout, peakbuf, -999.0, &status);
    ffpclb(fits, gcols[icol++], starout+1, 1, nmeasout, flagbuf, &status);

    for(ap = 0; ap < napcol; ap++)
      ffpcne(fits, gcols[icol++], starout+1, 1, 1, medfluxbuf+ap, -999.0, &status);

    for(ap = 0; ap < napcol; ap++)
      ffpcne(fits, gcols[icol++], starout+1, 1, 1, rmsbuf+ap, -999.0, &status);

    ffpcne(fits, gcols[icol++], starout+1, 1, 1, &chisq, -999.0, &status);
    ffpcnj(fits, gcols[icol++], starout+1, 1, 1, &nchisq, -999, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpcl");
      goto error;
    }

    starout++;
  }

  free((void *) gcols);
  gcols = (int *) NULL;
  free((void *) lcbuf);
  lcbuf = (struct lc_point *) NULL;
#ifdef HJD
  free((void *) epos);
  epos = (double *) NULL;
#endif
  free((void *) copysect);
  copysect = (struct table_section *) NULL;
  free((void *) medfluxbuf);
  medfluxbuf = (float *) NULL;
  free((void *) fluxbuf);
  fluxbuf = (float *) NULL;
  free((void *) airbuf);
  airbuf = (float *) NULL;
  free((void *) bjdbuf);
  bjdbuf = (double *) NULL;
  free((void *) flagbuf);
  flagbuf = (unsigned char *) NULL;
  free((void *) rawbuf);
  rawbuf = (unsigned char *) NULL;

  return(0);

 error:
  if(gcols)
    free((void *) gcols);
  if(lcbuf)
    free((void *) lcbuf);
#ifdef HJD
  if(epos)
    free((void *) epos);
#endif
  if(copysect)
    free((void *) copysect);
  if(medfluxbuf)
    free((void *) medfluxbuf);
  if(fluxbuf)
    free((void *) fluxbuf);
  if(airbuf)
    free((void *) airbuf);
  if(bjdbuf)
    free((void *) bjdbuf);
  if(flagbuf)
    free((void *) flagbuf);
  if(rawbuf)
    free((void *) rawbuf);

  return(1);
}

