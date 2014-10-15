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

static int write_lc (fitsfile *reff, fitsfile *fits,
		     struct buffer_info *buf, struct lc_mef *mefinfo,
		     int outcls, int wantoutcls,
		     char *errstr);
static int write_goodlist (char *outfile, struct lc_mef *meflist, int nmefs,
			   char **fnlist, char *errstr);

/* Getopt stuff */
extern char *optarg;
extern int optind;

/* Verbose flag */
int verbose = 1;

static void usage (char *av) {
  fprintf(stderr, "Usage:\t%s [options] reffile file [...]\n\n", av);
  fprintf(stderr,
	  "Lightcurve processing:\n"
	  "         -a aper   Specify aperture or 'sel' for only auto (default both).\n"
	  "         -d        Enables difference imaging mode.\n"
	  "         -f degree Apply polynomial of 'degree' for systematics removal.\n"
	  "         -i file   Apply intrapixel correction from 'file'.\n"
	  "         -m        Allow for meridian offset in rms and weighting.\n"
	  "         -mm       Same, only also removes it.\n"
	  "         -n        Do not renormalise median to reference magnitude.\n"
          "         -S        Use theoretical sky noise rather than empirical.\n"
	  "         -s level  Override saturation level to 'level'.\n"
	  "         -u mag    Set upper(,lower) mag limit for systematics correction.\n"
	  "         -V file   Use 'instrument version' table from 'file' (MEarth only).\n\n"
	  "Output:\n"
	  "         -c cls    Write out only comparison stars and class==cls.\n"
	  "         -g file   Writes good frames list to 'file'.  Cannot be used with -c.\n"
	  "         -o file   Writes lightcurves to 'file'.\n"
	  "         -p        Disables plots.\n"
	  "         -q        Decreases the verbosity level of the program.\n"
  	  "         -v        Increases the verbosity level of the program.\n");
  exit(1);
}  

int main (int argc, char *argv[]) {
  char *pn = (char *) NULL, *avzero, *p, *ep;
  int c;

  char errstr[ERRSTR_LEN];

  char *refname, **fnlist = (char **) NULL;
  struct lc_mef *meflist = (struct lc_mef *) NULL;
  struct intra *intralist = (struct intra *) NULL;
  struct buffer_info buf;

  long f, nf = 0;

  fitsfile *inf, *outf;
  int status = 0, ext, mef, nmefs;

  char outfile[FLEN_FILENAME-1], fnbuf[FLEN_FILENAME];
  int dooutput = 0;

  char goodfile[FLEN_FILENAME];
  int dogood = 0;

  char intrafile[FLEN_FILENAME];
  int dointra = 0;

  char instversfile[FLEN_FILENAME];
  int doinstvers = 0;
  struct instvers *instverslist = (struct instvers *) NULL;
  int ninstvers = 0;

  int aperture = 0;
  int apselmode = APSEL_SEL | APSEL_ALL;

  int domerid = 0;
  int theosky = 0;
  int norenorm = 0;
  int polydeg = -1;

  int outcls = 0;
  int wantoutcls = 0;

  int rv;
  struct dtai_table dtab, *dtptr = NULL;
  struct iers_table itab, *itptr = NULL;
  struct jpleph_table jtab, ttab, *jtptr = NULL, *ttptr = NULL;

  int len, maxflen, fspc;
  float *medbuf1 = (float *) NULL, *medbuf2, medsat, medlim;
  long nmedsat, nmedlim, nstartot;

  long star;
  int cflagmin;

  float satlev = -1.0;
  float sysulim = -1.0, sysllim = -1.0;

  int noplots = 0;
  int diffmode = 0;

  /* Set the program name for error reporting */
  if(argv[0])
    pn = basename(argv[0]);
  
  if(!pn)
    pn = "lightcurves";

  setprogname(pn);

  avzero = argv[0];

  /* Extract command-line arguments */
  while((c = getopt(argc, argv, "a:c:df:g:i:mno:pqSs:u:vV:")) != -1)
    switch(c) {
    case 'a':
      if(!strncasecmp(optarg, "sel", 3))
	apselmode = APSEL_SEL;  /* just selected aperture */
      else {
	aperture = (int) strtol(optarg, &ep, 0);
	if(*ep != '\0' || aperture < 0)
	  fatal(1, "invalid aperture: %s", optarg);
	apselmode = APSEL_SEL;
      }
      break;
    case 'c':
      outcls = (int) strtol(optarg, &ep, 0);
      if(*ep != '\0')
	fatal(1, "invalid class flag: %s", optarg);
      wantoutcls = 1;
      break;
    case 'd':
      diffmode++;
      break;
    case 'f':
      polydeg = (int) strtol(optarg, &ep, 10);
      if(*ep != '\0' || polydeg < 0)
	fatal(1, "invalid degree value: %s", optarg);
      break;
    case 'g':
      strncpy(goodfile, optarg, sizeof(goodfile)-1);
      goodfile[sizeof(goodfile)-1] = '\0';
      dogood = 1;
      break;
    case 'i':
      strncpy(intrafile, optarg, sizeof(intrafile)-1);
      intrafile[sizeof(intrafile)-1] = '\0';
      dointra = 1;
      break;
    case 'm':
      domerid++;
      break;
    case 'n':
      norenorm = 1;
      break;
    case 'o':
      strncpy(outfile, optarg, sizeof(outfile)-1);
      outfile[sizeof(outfile)-1] = '\0';
      dooutput = 1;
      break;
    case 'p':
      noplots++;
      break;
    case 'q':
      verbose--;
      break;
    case 'S':
      theosky++;
      break;
    case 's':
      satlev = (float) strtod(optarg, &ep);
      if(*ep != '\0' || satlev < 0)
	fatal(1, "invalid satlev value: %s", optarg);
      break;
    case 'u':
      sysulim = (float) strtod(optarg, &ep);
      if(ep == optarg || sysulim < 0)
	fatal(1, "invalid syslim value: %s", optarg);

      while(*ep && isspace((unsigned char) *ep))
	ep++;

      p = ep;
      if(*p == ',') {
	p++;
	
	sysllim = (float) strtod(p, &ep);
	if(ep == p || sysllim <= sysulim)
	  fatal(1, "invalid syslim value: %s", optarg);
      }

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

  fnlist = read_file_list(argc, argv, &nf, errstr);
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

  /* Create output file */
  if(dooutput) {
    /* Form output name */
    fnbuf[0] = '!';
    strncpy(&(fnbuf[1]), outfile, sizeof(fnbuf)-1);
    fnbuf[sizeof(fnbuf)-1] = '\0';

    /* Create file */
    ffinit(&outf, fnbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffinit: %s", outfile);
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
    meflist[mef].theosky = theosky;

    meflist[mef].degree = polydeg;
    meflist[mef].aperture = aperture;
    meflist[mef].apselmode = apselmode;
    meflist[mef].domerid = domerid;

    meflist[mef].avsigma = 0.0;
    meflist[mef].avskyfit = 0.0;
    meflist[mef].avapcor = 0.0;
    meflist[mef].avscint = 0.0;

    meflist[mef].frames = (struct lc_frame *) NULL;
    meflist[mef].nf = nf;

    if(ext == -99) {
      /* Move there */
      ffmahd(inf, mef+2, (int *) NULL, &status);
      if(status) {
	fitsio_err(errstr, status, "ffmahd: %s: HDU %d", refname, mef+2);
	fatal(1, "%s", errstr);
      }
    }

    /* Read it in */
    if(read_ref(inf, &(meflist[mef]), diffmode, satlev,
		sysllim, sysulim, outcls, wantoutcls, errstr))
      fatal(1, "read_ref: HDU %d: %s", mef+2, errstr);

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
	printf("\r Reading %*s (%*ld of %*ld)", maxflen, fnlist[f], fspc, f+1, fspc, nf);

      if(read_cat(fnlist[f], f, mef, &(meflist[mef]), &buf,
		  dtptr, itptr, jtptr, ttptr,
		  dointra, &(intralist[mef]),
		  doinstvers, instverslist, ninstvers,
		  diffmode, satlev, errstr))
	fatal(1, "read_cat: %s: %s", fnlist[f], errstr);
    }

    if(verbose && isatty(1))
      printf("\n");

    if(meflist[mef].nstars > 0) {
      /* Sort out averages */
      meflist[mef].avsigma = sqrtf(meflist[mef].avsigma / nf);
      meflist[mef].avskyfit = sqrtf(meflist[mef].avskyfit / nf);
      meflist[mef].avapcor /= nf;
      meflist[mef].avscint /= nf;
      
      /* Fix the cflag column - sometimes in difference imaging there
       * are frames with zero confidence all-over, so subtract off
       * the minimum cflag value.
       */
      cflagmin = meflist[mef].stars[0].cflag;

      for(star = 1; star < meflist[mef].nstars; star++)
        if(meflist[mef].stars[star].cflag < cflagmin)
          cflagmin = meflist[mef].stars[star].cflag;
      
      for(star = 0; star < meflist[mef].nstars; star++)
        meflist[mef].stars[star].cflag -= cflagmin;
    }
    else
      cflagmin = 0;

    /* Change MJD to be relative to the first frame */
    meflist[mef].mjdref = floor(meflist[mef].frames[0].mjd);

    for(f = 0; f < nf; f++)
      meflist[mef].frames[f].mjd -= meflist[mef].mjdref;

    if(verbose) {
      printf("  Number of objects:    %ld\n", meflist[mef].nstars);

      if(meflist[mef].satmag != -999.0)
	printf("  Saturation level:     %.1f\n", meflist[mef].satmag);
      else
	printf("  Saturation level:     undetermined\n");

      printf("  5-sigma limit:        %.1f\n"
	     "  Minimum cflag:        %d\n",
	     meflist[mef].refflim, cflagmin);
    }

    if(meflist[mef].satmag != -999.0) {
      medbuf1[nmedsat] = meflist[mef].satmag;
      nmedsat++;
    }

    medbuf2[nmedlim] = meflist[mef].refflim;
    nmedlim++;

    nstartot += meflist[mef].nstars;

    /* Call into the main part of the program */
    if(lightcurves(&buf, &(meflist[mef]), norenorm, wantoutcls, errstr))
      fatal(1, "%s", errstr);

    /* Calculate average extinction */
    if(meflist[mef].degree >= 0) {
      meflist[mef].avextinc = 0.0;
      meflist[mef].avsigm = 0.0;

      for(f = 0; f < nf; f++) {
	meflist[mef].avextinc += powf(10.0, -0.4 * meflist[mef].frames[f].extinc);
	meflist[mef].avsigm += meflist[mef].frames[f].sigm;
      }      

      meflist[mef].avextinc /= nf;
      meflist[mef].avsigm /= nf;
    }
    else {
      meflist[mef].avextinc = 1.0;
      meflist[mef].avsigm = 0.0;
    }

    /* Write out lightcurves for this MEF if requested */
    if(dooutput) {
      if(verbose)
	printf(" Writing %s\n", outfile);

      if(write_lc(inf, outf, &buf, &(meflist[mef]), outcls, wantoutcls, errstr))
	fatal(1, "write_lc: %s", errstr);
    }

    /* Flush */
    if(dooutput) {
      ffflus(outf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffflus");
	fatal(1, "%s", errstr);
      }
    }
  }

#ifdef DEBUG
  cpgclos();
#endif

  /* Release disk buffer */
  buffer_close(&buf);

  /* Close reference */
  ffclos(inf, &status);
  if(status) {
    fitsio_err(errstr, status, "ffclos");
    fatal(1, "%s", errstr);
  }

  /* Close output file */
  if(dooutput) {
    ffclos(outf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffclos");
      fatal(1, "%s", errstr);
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

  if(dogood) {
    if(write_goodlist(goodfile, meflist, nmefs, fnlist, errstr))
      fatal(1, "do_goodlist: %s");
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

static int write_lc (fitsfile *reff, fitsfile *fits,
		     struct buffer_info *buf, struct lc_mef *mefinfo,
		     int outcls, int wantoutcls,
		     char *errstr) {
  int status = 0, col, icol, ocol, ntmpl, ncols;

  char tsbuf[FLEN_VALUE], tfbuf[FLEN_VALUE], tdbuf[FLEN_VALUE], tbbuf[FLEN_VALUE];

  /* Column templates */
  char *tmpl_ttype[] = { "x", "y", "medflux", "rms", "chisq", "nchisq",
			 "class", "bflag", "cflag", "sflag", "pointer",
			 "offsets", "apnum", "apradius", "compok",
			 "bjd",
#ifdef HJD
                         "hjd",
#endif
                         "flux", "fluxerr", "xlc", "ylc", "airmass", "ha",
			 "weight", "sky", "peak", "flags",
			 "ra", "dec", "pmra", "pmdec", "refmag" };
  char *tmpl_tform[] = { "1D", "1D", "1E", "1E", "1E", "1J",
			 "1I", "1I", "1J", "1J", "1J",
			 tsbuf, "1I", "1E", "1B",
			 tdbuf,
#ifdef HJD
			 tdbuf,
#endif
                         tfbuf, tfbuf, tdbuf, tdbuf, tfbuf, tfbuf,
			 tfbuf, tfbuf, tfbuf, tbbuf,
			 "1D", "1D", "1E", "1E", "1E" };
  char *tmpl_tunit[] = { "pixels", "pixels", "mag", "mag", "", "",
			 "", "", "", "", "",
			 "mag", "", "pixels", "",
			 "days",
#ifdef HJD
			 "days",
#endif
                         "mag", "mag", "pixels", "pixels", "", "radians",
			 "", "counts", "counts", "",
			 "radians", "radians", "arcsec/yr", "arcsec/yr", "mag" };
  char *tmpl_tdisp[] = { "F8.2", "F8.2", "F7.4", "F7.4", "F10.1", "I4",
			 "I2", "I2", "I4", "I8", "I8",
			 "F7.4", "", "F4.2", "I1",
			 "F14.6",
#ifdef HJD
			 "F14.6",
#endif
                         "F7.4", "F7.4", "F8.2", "F8.2", "F6.4", "F9.6",
			 "F9.0", "F8.2", "F5.0", "I3",
			 "F9.6", "F9.6", "F7.3", "F7.3", "F7.4" };

  /* Repeat column for every aperture?  1 = always, 2 = only if writing all requested */
  unsigned char tpap[] = { 0, 0, 1, 1, 0, 0,
			   0, 0, 0, 0, 0,
			   1, 0, 0, 0,
			   0,
#ifdef HJD
			   0,
#endif
                           2, 2, 0, 0, 0, 0,
			   2, 0, 0, 0,
			   0, 0, 0, 0, 0 };

  /* Real arrays, we build these later */
  char **ttype = (char **) NULL, **tform, **tunit, **tdisp;

  char kbuf[FLEN_KEYWORD];
  char vbuf[FLEN_VALUE];
  char cbuf[FLEN_COMMENT];

  long pt, star, nstarout;

  struct lc_point *lcbuf = (struct lc_point *) NULL;

#ifdef HJD
  double *epos = (double *) NULL;
#endif

  double *xbuf = (double *) NULL, *ybuf, *rabuf, *decbuf;
  float *medbuf = (float *) NULL, *rmsbuf;
  float *apbuf = (float *) NULL, *chibuf, *refmagbuf, *pmabuf, *pmdbuf;
  long *ptrbuf = (long *) NULL, *cfbuf, *sfbuf, *nchibuf;
  short *clsbuf = (short *) NULL, *bfbuf, *apnumbuf;
  float *offbuf = (float *) NULL;
  unsigned char *compokbuf = (unsigned char *) NULL;
  float *fluxbuf = (float *) NULL, *fluxerrbuf, *wtbuf;
  float *airbuf = (float *) NULL, *habuf, *locskybuf, *peakbuf;
  double *bjdbuf = (double *) NULL, *xlcbuf, *ylcbuf;
#ifdef HJD
  double *hjdbuf;
#endif
  unsigned char *flagbuf = (unsigned char *) NULL;

  int ikey, nkeys, kclass;
  char card[FLEN_CARD];
  long r, frow, rblksz, soff;

  long satflag;
  unsigned char flags;

  int ap, ap1, ap2, lap1, lap2, napcol, nlapcol;
  int allast;

  int iseg;

  /* Figure out apertures */
  if(mefinfo->aperture) {
    ap1 = mefinfo->aperture-1;
    ap2 = ap1;
    napcol = 1;
    lap1 = ap1;
    lap2 = ap2;
    nlapcol = 1;
  }
  else {
    ap1 = 0;
    ap2 = NFLUX;
    napcol = ap2-ap1 + 1;

    if(mefinfo->apselmode & APSEL_ALL) {
      lap1 = ap1;
      lap2 = ap2;
      nlapcol = napcol;
    }
    else {
      lap1 = 0;
      lap2 = 0;
      nlapcol = 1;
    }
  }

  /* Generate tform specifier for fluxes and errors */
  snprintf(tsbuf, sizeof(tsbuf), "%ldE", mefinfo->nseg);
  snprintf(tdbuf, sizeof(tdbuf), "%ldD", mefinfo->nf);
  snprintf(tfbuf, sizeof(tfbuf), "%ldE", mefinfo->nf);
  snprintf(tbbuf, sizeof(tbbuf), "%ldB", mefinfo->nf);

  /* Figure out real number of columns */
  ntmpl = sizeof(tmpl_ttype) / sizeof(tmpl_ttype[0]);

  ncols = 0;
  for(icol = 0; icol < ntmpl; icol++) {
    ncols++;

    if(tpap[icol] == 1 || (tpap[icol] == 2 && mefinfo->apselmode & APSEL_ALL))
      ncols += ap2-ap1;
  }
   
  /* Allocate workspace */
  ttype = (char **) malloc(4 * ncols * sizeof(char *));
  if(!ttype) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  tform = ttype + ncols;
  tunit = ttype + 2*ncols;
  tdisp = ttype + 3*ncols;

  /* Make sure populated with null pointers */
  ocol = 0;

  for(icol = 0; icol < ntmpl; icol++) {
    ocol++;

    if(tpap[icol] == 1 || (tpap[icol] == 2 && mefinfo->apselmode & APSEL_ALL))
      for(ap = ap1; ap < ap2; ap++) {
	ttype[ocol] = NULL;
	ocol++;
      }
  }

  /* Build real column descriptors */
  ocol = 0;

  for(icol = 0; icol < ntmpl; icol++) {
    ttype[ocol] = tmpl_ttype[icol];
    tform[ocol] = tmpl_tform[icol];
    tunit[ocol] = tmpl_tunit[icol];
    tdisp[ocol] = tmpl_tdisp[icol];
    
    ocol++;

    if(tpap[icol] == 1 || (tpap[icol] == 2 && mefinfo->apselmode & APSEL_ALL))
      for(ap = ap1; ap < ap2; ap++) {
	snprintf(vbuf, sizeof(vbuf), "%s%d", tmpl_ttype[icol], ap+1);
	ttype[ocol] = strdup(vbuf);
	if(!ttype[ocol]) {
	  report_syserr(errstr, "malloc");
	  goto error;
	}

	tform[ocol] = tmpl_tform[icol];
	tunit[ocol] = tmpl_tunit[icol];
	tdisp[ocol] = tmpl_tdisp[icol];

	ocol++;
      }
  }

  /* Figure out how large the table is */
  if(wantoutcls) {
    nstarout = 0;
    for(star = 0; star < mefinfo->nstars; star++) {
      if(mefinfo->stars[star].compok ||
	 mefinfo->stars[star].cls == outcls)
	nstarout++;
    }
  }
  else
    nstarout = mefinfo->nstars;

  /* Create table */
  ffcrtb(fits, BINARY_TBL, nstarout, ncols, ttype, tform, tunit, "", &status);
  if(status) {
    for(col = 0; col < ncols; col++)
      fprintf(stderr, "%d %s %s %s\n", col+1, ttype[col], tform[col], tunit[col]);

    fitsio_err(errstr, status, "ffcrtb");
    goto error;
  }

  /* Write TDISP keywords */
  for(col = 0; col < ncols; col++)
    if(tdisp[col] && *(tdisp[col]) != '\0') {
      snprintf(kbuf, sizeof(kbuf), "TDISP%d", col+1);
      snprintf(cbuf, sizeof(cbuf), "Display format for column %d", col+1);
      ffpkys(fits, kbuf, tdisp[col], cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkys: %s", kbuf);
	goto error;
      }
    }

  /* Dispense with workspace */
  ocol = 0;

  for(icol = 0; icol < ntmpl; icol++) {
    ocol++;

    if(tpap[icol] == 1 || (tpap[icol] == 2 && mefinfo->apselmode & APSEL_ALL))
      for(ap = ap1; ap < ap2; ap++) {
	free((void *) ttype[ocol]);
	ocol++;
      }
  }

  free((void *) ttype);
  ttype = (char **) NULL;

#ifdef HJD
  /* Allocate buffer for earth positions */
  epos = (double *) malloc(3 * mefinfo->nf * sizeof(double));
  if(!epos) {
    report_syserr(errstr, "malloc");
    goto error;
  }
#endif

  /* Write out frame information */
  ffpkyj(fits, "NMEAS", mefinfo->nf,
	 "Number of points in each lightcurve", &status);
  ffpkyj(fits, "NROWMAST", mefinfo->nrows,
	 "Number of rows in master catalogue", &status);
  ffpkyg(fits, "MJDBASE", mefinfo->mjdref, 7,
	 "Base MJD for time axis", &status);
  ffpkyf(fits, "SATMAG", mefinfo->satmag, 4,
	 "Approximate saturation magnitude", &status);
  ffpkyf(fits, "FLIM", mefinfo->refflim, 4,
	 "Flux limit of reference catalogue", &status);
  ffpkyf(fits, "ZP", mefinfo->zp, 4,
	 "Zeropoint for magnitudes", &status);
  ffpkyf(fits, "UMLIM", mefinfo->zp - mefinfo->sysulim, 4,
	 "Upper mag limit for fit", &status);
  ffpkyf(fits, "LMLIM", mefinfo->zp - mefinfo->sysllim, 4,
	 "Lower mag limit for fit", &status);
  ffpkyl(fits, "THEOSKY", mefinfo->theosky,
         "T theoretical sky noise, F empirical", &status);
  ffpkyj(fits, "POLYDEG", mefinfo->degree,
	 "Polynomial degree in fit", &status);
  ffpkyj(fits, "APSEL", mefinfo->aperture,
	 "Aperture used (0 = automatic)", &status);
  ffpkyj(fits, "APMODE", mefinfo->apselmode,
	 "Aperture output mode", &status);
  ffpkyj(fits, "DOMERID", mefinfo->domerid,
	 "Meridian flip removal?", &status);
  ffpkyf(fits, "REFFANG", mefinfo->reffang, 6,
	 "Reference file field angle", &status);
  ffpkyj(fits, "NSEGME", mefinfo->nseg,
	 "Number of segments", &status);
  if(status) {
    fitsio_err(errstr, status, "ffkpy: frame info");
    goto error;
  }

  for(iseg = 0; iseg < mefinfo->nseg; iseg++) {
    snprintf(kbuf, sizeof(kbuf), "SEGV%d", iseg+1);
    snprintf(cbuf, sizeof(cbuf), "Segment %d instrument version number", iseg+1);
    ffpkyj(fits, kbuf,
	   mefinfo->segs[iseg].instvers ?
	   mefinfo->segs[iseg].instvers->iver : -1, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "SEGD%d", iseg+1);
    snprintf(cbuf, sizeof(cbuf), "Segment %d instrument change date", iseg+1);
    ffpkyj(fits, kbuf,
	   mefinfo->segs[iseg].instvers ?
	   mefinfo->segs[iseg].instvers->date : -1, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "SEGA%d", iseg+1);
    snprintf(cbuf, sizeof(cbuf), "Segment %d angle", iseg+1);
    ffpkyj(fits, kbuf, mefinfo->segs[iseg].iang, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
      goto error;
    }
  }

  allast = 1;

  for(pt = 0; pt < mefinfo->nf; pt++) {
    snprintf(kbuf, sizeof(kbuf), "TV%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Time value for datapoint %ld", pt+1);
    ffpkyg(fits, kbuf, mefinfo->frames[pt].mjd, 14, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyg: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "TEXP%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Exposure time for datapoint %ld", pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].exptime, 3, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyg: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "OFF%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Frame offset for datapoint %ld", pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].offset, 4, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "RMS%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Frame RMS for datapoint %ld", pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].rms, 4, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "EXTC%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Extinction correction for datapoint %ld", pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].extinc, 4, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "SEE%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Seeing for datapoint %ld", pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].seeing, 3, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "ELL%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Ellipticity for datapoint %ld", pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].ellipt, 3, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "SKY%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Sky level for datapoint %ld", pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].skylev, 2, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "NOIS%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Sky noise for datapoint %ld", pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].skynoise, 2, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "FANG%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Field angle for datapoint %ld", pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].fang, 6, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    snprintf(kbuf, sizeof(kbuf), "IANG%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Field angle modulo pi for datapoint %ld", pt+1);
    ffpkyj(fits, kbuf, mefinfo->frames[pt].iang, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    if(mefinfo->frames[pt].tamb != -999) {
      snprintf(kbuf, sizeof(kbuf), "TAMB%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "[degC] Ambient temp for datapoint %ld", pt+1);
      ffpkyf(fits, kbuf, mefinfo->frames[pt].tamb, 1, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].humid != -999) {
      snprintf(kbuf, sizeof(kbuf), "HUM%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "[%%] Humidity for datapoint %ld", pt+1);
      ffpkyf(fits, kbuf, mefinfo->frames[pt].humid, 1, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].press != -999) {
      snprintf(kbuf, sizeof(kbuf), "PRES%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "[hPa] Pressure for datapoint %ld", pt+1);
      ffpkyf(fits, kbuf, mefinfo->frames[pt].press, 1, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].skytemp != -999) {
      snprintf(kbuf, sizeof(kbuf), "TSKY%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "[degC] Sky temp for datapoint %ld", pt+1);
      ffpkyf(fits, kbuf, mefinfo->frames[pt].skytemp, 1, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].split_nexp >= 0) {
      snprintf(kbuf, sizeof(kbuf), "IEXP%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Exposure index for datapoint %ld", pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].split_iexp, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }

      snprintf(kbuf, sizeof(kbuf), "NEXP%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Exposure count for datapoint %ld", pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].split_nexp, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].rtstat >= 0) {
      snprintf(kbuf, sizeof(kbuf), "RTST%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Trigger status for datapoint %ld", pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].rtstat, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].schpri >= 0) {
      snprintf(kbuf, sizeof(kbuf), "SPRI%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Scheduling priority for datapoint %ld", pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].schpri, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].schcad >= 0) {
      snprintf(kbuf, sizeof(kbuf), "SCAD%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Scheduling cadence for datapoint %ld", pt+1);
      ffpkye(fits, kbuf, mefinfo->frames[pt].schcad, 6, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkye: %s", kbuf);
	goto error;
      }
    }
    
    if(mefinfo->frames[pt].schtype[0]) {
      snprintf(kbuf, sizeof(kbuf), "STYP%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Observation type code for datapoint %ld", pt+1);
      ffpkys(fits, kbuf, mefinfo->frames[pt].schtype, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkys: %s", kbuf);
	goto error;
      }
    }

    if(mefinfo->frames[pt].cadencenum != -999) {
      snprintf(kbuf, sizeof(kbuf), "CADN%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Cadence number for datapoint %ld", pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].cadencenum, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }
    }

    snprintf(kbuf, sizeof(kbuf), "ISEG%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Segment number for datapoint %ld", pt+1);
    ffpkyj(fits, kbuf, mefinfo->frames[pt].iseg+1, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
      goto error;
    }

    if(mefinfo->frames[pt].instvers) {
      snprintf(kbuf, sizeof(kbuf), "IVER%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Instrument version for datapoint %ld", pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].instvers->iver, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }

      snprintf(kbuf, sizeof(kbuf), "IDAT%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Inst last change date for datapoint %ld", pt+1);
      ffpkyj(fits, kbuf, mefinfo->frames[pt].instvers->date, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
	goto error;
      }
    }

    if(strcmp(mefinfo->frames[pt].schtype, "a"))
      allast = 0;

    snprintf(kbuf, sizeof(kbuf), "IUPD%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Update number when datapoint %ld was added", pt+1);
    ffpkyj(fits, kbuf, 0, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyj: %s", kbuf);
      goto error;
    }

    /* Suppress output of frame transformation if we only had the
       photometric comparison stars.  It really needs everything. */
    if(!wantoutcls) {
      snprintf(kbuf, sizeof(kbuf), "LXX%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Linear transformation to ref. for frame %ld", pt+1);
      ffpkyd(fits, kbuf, mefinfo->frames[pt].tr[0], 12, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyd: %s", kbuf);
	goto error;
      }
      
      snprintf(kbuf, sizeof(kbuf), "LXY%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Linear transformation to ref. for frame %ld", pt+1);
      ffpkyd(fits, kbuf, mefinfo->frames[pt].tr[1], 12, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyd: %s", kbuf);
	goto error;
      }
      
      snprintf(kbuf, sizeof(kbuf), "LXD%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Linear transformation to ref. for frame %ld", pt+1);
      ffpkyd(fits, kbuf, mefinfo->frames[pt].tr[2], 12, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyd: %s", kbuf);
	goto error;
      }
      
      snprintf(kbuf, sizeof(kbuf), "LYY%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Linear transformation to ref. for frame %ld", pt+1);
      ffpkyd(fits, kbuf, mefinfo->frames[pt].tr[3], 12, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyd: %s", kbuf);
	goto error;
      }
      
      snprintf(kbuf, sizeof(kbuf), "LYX%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Linear transformation to ref. for frame %ld", pt+1);
      ffpkyd(fits, kbuf, mefinfo->frames[pt].tr[4], 12, cbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpkyd: %s", kbuf);
	goto error;
      }
      
      snprintf(kbuf, sizeof(kbuf), "LYD%ld", pt+1);
      snprintf(cbuf, sizeof(cbuf), "Linear transformation to ref. for frame %ld", pt+1);
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

  /* Write out astrometry flag for web page */
  ffpkyl(fits, "ASTONLY", allast, "Are all observations for astrometry?", &status);
  if(status) {
    fitsio_err(errstr, status, "ffpkyl: ASTONLY");
    goto error;
  }

  /* Copy in reference header keywords */
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

    kclass = ffgkcl(card);
    switch(kclass) {
    case TYP_STRUC_KEY:
    case TYP_CMPRS_KEY:
    case TYP_SCAL_KEY:
    case TYP_NULL_KEY:
    case TYP_DIM_KEY:
    case TYP_RANG_KEY:
    case TYP_UNIT_KEY:
    case TYP_DISP_KEY:
    case TYP_HDUID_KEY:
    case TYP_CKSUM_KEY:
      break;
    default:
      ffprec(fits, card, &status);
      if(status) {
	fitsio_err(errstr, status, "ffprec: card %d", ikey+1);
	goto error;
      }
    }
  }

  /* Get optimal block size for writing */
  ffgrsz(fits, &rblksz, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgrsz");
    goto error;
  }

  /* Allocate input buffer */
  lcbuf = (struct lc_point *) malloc(mefinfo->nf * sizeof(struct lc_point));
  if(!lcbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Allocate output buffers */
  xbuf = (double *) malloc(4 * rblksz * sizeof(double));
  medbuf = (float *) malloc(2 * rblksz * napcol * sizeof(float));
  apbuf = (float *) malloc(5 * rblksz * sizeof(float));
  ptrbuf = (long *) malloc(4 * rblksz * sizeof(long));
  clsbuf = (short *) malloc(3 * rblksz * sizeof(short));
  offbuf = (float *) malloc(mefinfo->nseg * napcol * rblksz * sizeof(float));
  compokbuf = (unsigned char *) malloc(rblksz * sizeof(unsigned char));
#ifdef HJD
  bjdbuf = (double *) malloc(4 * rblksz * mefinfo->nf * sizeof(double));
#else
  bjdbuf = (double *) malloc(3 * rblksz * mefinfo->nf * sizeof(double));
#endif
  fluxbuf = (float *) malloc(3 * rblksz * mefinfo->nf * nlapcol * sizeof(float));
  airbuf = (float *) malloc(4 * rblksz * mefinfo->nf * sizeof(float));
  flagbuf = (unsigned char *) malloc(rblksz * mefinfo->nf * sizeof(unsigned char));
  if(!xbuf || !medbuf || !apbuf || !ptrbuf || !clsbuf || !offbuf || !compokbuf || !bjdbuf || !fluxbuf || !airbuf || !flagbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  rmsbuf = medbuf + napcol * rblksz;

  chibuf = apbuf + rblksz;
  refmagbuf = apbuf + 2 * rblksz;
  pmabuf = apbuf + 3 * rblksz;
  pmdbuf = apbuf + 4 * rblksz;

  ybuf = xbuf + rblksz;
  rabuf = xbuf + 2 * rblksz;
  decbuf = xbuf + 3 * rblksz;

  bfbuf = clsbuf + rblksz;
  apnumbuf = clsbuf + 2*rblksz;

  cfbuf = ptrbuf + rblksz;
  sfbuf = ptrbuf + 2 * rblksz;
  nchibuf = ptrbuf + 3 * rblksz;

  fluxerrbuf = fluxbuf + nlapcol * rblksz * mefinfo->nf;
  wtbuf = fluxbuf + 2 * nlapcol * rblksz * mefinfo->nf;  

  habuf = airbuf + rblksz * mefinfo->nf;
  locskybuf = airbuf + 2 * rblksz * mefinfo->nf;
  peakbuf = airbuf + 3 * rblksz * mefinfo->nf;

  xlcbuf = bjdbuf + rblksz * mefinfo->nf;
  ylcbuf = bjdbuf + 2 * rblksz * mefinfo->nf;
#ifdef HJD
  hjdbuf = bjdbuf + 3 * rblksz * mefinfo->nf;
#endif

  /* Loop through all stars - read in the lightcurve points one
   * at a time and write out blocks of 'rblksz' objects.
   */
  r = 0;  /* number of rows buffered */
  frow = 1;

#define FORAP for(ap = 0; ap < napcol; ap++)
#define FORLAP for(ap = 0; ap < nlapcol; ap++)
#define WRITE_COL(func, buf, len)			\
  func(fits, ++ocol, frow, 1, r*(len), buf, &status)
#define WRITE_COL_NULL(func, buf, len, nullval)		\
  func(fits, ++ocol, frow, 1, r*(len), buf, nullval, &status)
#ifdef HJD
#define WRITE_HJD WRITE_COL_NULL(ffpcnd, hjdbuf, mefinfo->nf, -999.0);
#else
#define WRITE_HJD
#endif

#define TABLE_FLUSH() {							\
  ocol = 0;								\
  WRITE_COL(ffpcld, xbuf, 1);						\
  WRITE_COL(ffpcld, ybuf, 1);						\
  FORAP WRITE_COL_NULL(ffpcne, medbuf+ap*rblksz, 1, -999.0);		\
  FORAP WRITE_COL_NULL(ffpcne, rmsbuf+ap*rblksz, 1, -999.0);		\
  WRITE_COL(ffpcle, chibuf, 1);						\
  WRITE_COL(ffpclj, nchibuf, 1);					\
  WRITE_COL(ffpcli, clsbuf, 1);						\
  WRITE_COL(ffpcli, bfbuf, 1);						\
  WRITE_COL(ffpclj, cfbuf, 1);						\
  WRITE_COL(ffpclj, sfbuf, 1);						\
  WRITE_COL(ffpclj, ptrbuf, 1);						\
  FORAP WRITE_COL_NULL(ffpcne, offbuf+ap*rblksz*mefinfo->nseg, mefinfo->nseg, -999.0);	\
  WRITE_COL(ffpcli, apnumbuf, 1);					\
  WRITE_COL(ffpcle, apbuf, 1);						\
  WRITE_COL(ffpclb, compokbuf, 1);		       			\
  WRITE_COL_NULL(ffpcnd, bjdbuf, mefinfo->nf, -999.0);			\
  WRITE_HJD                                           			\
  FORLAP WRITE_COL_NULL(ffpcne, fluxbuf+ap*rblksz*mefinfo->nf, mefinfo->nf, -999.0); \
  FORLAP WRITE_COL_NULL(ffpcne, fluxerrbuf+ap*rblksz*mefinfo->nf, mefinfo->nf, -999.0); \
  WRITE_COL_NULL(ffpcnd, xlcbuf, mefinfo->nf, -999.0);			\
  WRITE_COL_NULL(ffpcnd, ylcbuf, mefinfo->nf, -999.0);			\
  WRITE_COL_NULL(ffpcne, airbuf, mefinfo->nf, -999.0);			\
  WRITE_COL_NULL(ffpcne, habuf, mefinfo->nf, -999.0);			\
  FORLAP WRITE_COL_NULL(ffpcne, wtbuf+ap*rblksz*mefinfo->nf, mefinfo->nf, -999.0); \
  WRITE_COL_NULL(ffpcne, locskybuf, mefinfo->nf, -999.0);		\
  WRITE_COL_NULL(ffpcne, peakbuf, mefinfo->nf, -999.0);			\
  WRITE_COL(ffpclb, flagbuf, mefinfo->nf);				\
  WRITE_COL(ffpcld, rabuf, 1);						\
  WRITE_COL(ffpcld, decbuf, 1);						\
  WRITE_COL_NULL(ffpcne, pmabuf, 1, -999.0);                            \
  WRITE_COL_NULL(ffpcne, pmdbuf, 1, -999.0);                            \
  WRITE_COL_NULL(ffpcne, refmagbuf, 1, -999.0);	       			\
  if(status) {								\
    fitsio_err(errstr, status, "ffpcl");				\
    goto error;								\
  }									\
									\
  frow += r;								\
  r = 0;								\
}

  for(star = 0; star < mefinfo->nstars; star++) {
    if(wantoutcls &&
       !mefinfo->stars[star].compok &&
       mefinfo->stars[star].cls != outcls)
      continue;  /* skip star */

    /* Fill in buffers */
    xbuf[r] = mefinfo->stars[star].x;
    ybuf[r] = mefinfo->stars[star].y;
    medbuf[r] = (mefinfo->stars[star].med > 0.0 ?
		 mefinfo->zp - mefinfo->stars[star].med : -999.0);
    rmsbuf[r] = mefinfo->stars[star].rms;

    for(iseg = 0; iseg < mefinfo->nseg; iseg++)
      offbuf[r*mefinfo->nseg + iseg]
	= mefinfo->stars[star].segs[iseg].corr[mefinfo->stars[star].iap];
    

    for(ap = ap1; ap < ap2; ap++) {
      medbuf[(ap-ap1+1)*rblksz+r] = (mefinfo->stars[star].medflux[ap] > 0.0 ?
				     mefinfo->zp - mefinfo->stars[star].medflux[ap] :
				     -999.0);
      rmsbuf[(ap-ap1+1)*rblksz+r] = mefinfo->stars[star].sigflux[ap];

      for(iseg = 0; iseg < mefinfo->nseg; iseg++)
	offbuf[((ap-ap1+1)*rblksz + r) * mefinfo->nseg + iseg]
	  = mefinfo->stars[star].segs[iseg].corr[ap];
    }

    chibuf[r] = mefinfo->stars[star].chisq;
    nchibuf[r] = mefinfo->stars[star].nchisq;

    clsbuf[r] = mefinfo->stars[star].cls;
    bfbuf[r] = mefinfo->stars[star].bflag;
    cfbuf[r] = mefinfo->stars[star].cflag;
    ptrbuf[r] = mefinfo->stars[star].ptr;

    apnumbuf[r] = mefinfo->stars[star].iap+1;
    apbuf[r] = mefinfo->stars[star].apradius;
    compokbuf[r] = mefinfo->stars[star].compok;
    rabuf[r] = mefinfo->stars[star].ra;
    decbuf[r] = mefinfo->stars[star].dec;

    if(mefinfo->stars[star].havepm) {
      pmabuf[r] = mefinfo->stars[star].pmra * RAD_TO_AS;
      pmdbuf[r] = mefinfo->stars[star].pmdec * RAD_TO_AS;
    }
    else {
      pmabuf[r] = -999.0;
      pmdbuf[r] = -999.0;
    }

    refmagbuf[r] = (mefinfo->stars[star].refmag > 0.0 ?
		    mefinfo->zp - mefinfo->stars[star].refmag : -999.0);

    /* Get lightcurve */
    if(buffer_fetch_object(buf, lcbuf, 0, mefinfo->nf, star, errstr))
      goto error;

    /* Fill in buffer */
    soff = r * mefinfo->nf;

    satflag = 0;
    for(pt = 0; pt < mefinfo->nf; pt++) {
      flags = 0;

      if(lcbuf[pt].aper[mefinfo->stars[star].iap].flux != 0.0) {
	fluxbuf[soff+pt] = mefinfo->zp - lcbuf[pt].aper[mefinfo->stars[star].iap].flux;

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
	  fluxerrbuf[soff+pt] = sqrtf(lcbuf[pt].aper[mefinfo->stars[star].iap].fluxvarcom);
	else
	  fluxerrbuf[soff+pt] = -999.0;
      }
      else {
	fluxbuf[soff+pt] = -999.0;
	fluxerrbuf[soff+pt] = -999.0;
	flags |= FLAG_NODP;
      }

      xlcbuf[soff+pt] = lcbuf[pt].x;
      ylcbuf[soff+pt] = lcbuf[pt].y;
      airbuf[soff+pt] = lcbuf[pt].airmass;
      habuf[soff+pt] = lcbuf[pt].ha;
      wtbuf[soff+pt] = lcbuf[pt].aper[mefinfo->stars[star].iap].wt;
      locskybuf[soff+pt] = lcbuf[pt].sky;
      peakbuf[soff+pt] = lcbuf[pt].peak;

      flagbuf[soff+pt] = flags;

#ifdef HJD
      /* Calculate HJD (as UTC) */
      hjdbuf[soff+pt] = mefinfo->mjdref + mefinfo->frames[pt].mjd +
	                hjdcorr(epos + 3*pt,
				mefinfo->stars[star].ra,
				mefinfo->stars[star].dec);
#endif
      bjdbuf[soff+pt] = lcbuf[pt].bjd;
    }

    sfbuf[r] = satflag;

    /* Do other apertures */
    for(ap = lap1; ap < lap2; ap++) {
      soff = ((ap-lap1+1)*rblksz + r) * mefinfo->nf;
      
      /* Fill in buffer */
      for(pt = 0; pt < mefinfo->nf; pt++) {
	if(lcbuf[pt].aper[ap].flux != 0.0) {
	  fluxbuf[soff+pt] = mefinfo->zp - lcbuf[pt].aper[ap].flux;
	  if(lcbuf[pt].aper[ap].fluxvarcom > 0.0)
	    fluxerrbuf[soff+pt] = sqrtf(lcbuf[pt].aper[ap].fluxvarcom);
	  else
	    fluxerrbuf[soff+pt] = -999.0;
	}

	wtbuf[soff+pt] = lcbuf[pt].aper[ap].wt;
      }
    }

    r++;

    if(r >= rblksz) {
      /* Flush */
      TABLE_FLUSH();
    }
  }

  /* Flush out buffers */
  if(r > 0) {
    TABLE_FLUSH();
  }    

  free((void *) lcbuf);
  lcbuf = (struct lc_point *) NULL;
#ifdef HJD
  free((void *) epos);
  epos = (double *) NULL;
#endif
  free((void *) medbuf);
  medbuf = (float *) NULL;
  free((void *) apbuf);
  apbuf = (float *) NULL;
  free((void *) xbuf);
  xbuf = (double *) NULL;
  free((void *) ptrbuf);
  ptrbuf = (long *) NULL;
  free((void *) clsbuf);
  clsbuf = (short *) NULL;
  free((void *) offbuf);
  offbuf = (float *) NULL;
  free((void *) compokbuf);
  compokbuf = (unsigned char *) NULL;
  free((void *) fluxbuf);
  fluxbuf = (float *) NULL;
  free((void *) airbuf);
  airbuf = (float *) NULL;
  free((void *) bjdbuf);
  bjdbuf = (double *) NULL;
  free((void *) flagbuf);
  flagbuf = (unsigned char *) NULL;

  return(0);

 error:
  if(ttype) {
    ocol = 0;
    
    for(icol = 0; icol < ntmpl; icol++) {
      ocol++;
      
      if(tpap[icol] == 1 || (tpap[icol] == 2 && mefinfo->apselmode & APSEL_ALL))
	for(ap = ap1; ap < ap2; ap++) {
	  if(ttype[ocol])
	    free((void *) ttype[ocol]);

	  ocol++;
	}
    }
    
    free((void *) ttype);
  }
  if(lcbuf)
    free((void *) lcbuf);
#ifdef HJD
  if(epos)
    free((void *) epos);
#endif
  if(medbuf)
    free((void *) medbuf);
  if(apbuf)
    free((void *) apbuf);
  if(xbuf)
    free((void *) xbuf);
  if(ptrbuf)
    free((void *) ptrbuf);
  if(clsbuf)
    free((void *) clsbuf);
  if(offbuf)
    free((void *) offbuf);
  if(compokbuf)
    free((void *) compokbuf);
  if(fluxbuf)
    free((void *) fluxbuf);
  if(airbuf)
    free((void *) airbuf);
  if(bjdbuf)
    free((void *) bjdbuf);
  if(flagbuf)
    free((void *) flagbuf);

  return(1);
}

static int write_goodlist (char *outfile, struct lc_mef *meflist, int nmefs,
			   char **fnlist, char *errstr) {
  FILE *fp;
  long f, nf;
  int rv, mef, isok, iseg;

  long nok;

  /* Create output file */
  fp = fopen(outfile, "w");
  if(!fp) {
    report_syserr(errstr, "open: %s", outfile);
    goto error;
  }

  /* Loop through files, decide if we should include each one */
  nf = meflist[0].nf;

  nok = 0;
  for(f = 0; f < nf; f++) {
    /* Check each MEF */
    isok = 1;

    for(mef = 0; mef < nmefs; mef++) {
#if 1
      /* Old criterion: residual offset or fit rms > 0.05 */
      if(fabsf(meflist[mef].frames[f].offset) > 0.05 ||
	 fabsf(meflist[mef].frames[f].rms) > 0.05)
	isok = 0;
#endif

#if 1
      iseg = meflist[mef].frames[f].iseg;

      /* MEarth criterion: delta mag > 0.6 and position within 10 sigma */
      if(meflist[mef].frames[f].extinc < -0.6 ||
	 fabsf(meflist[mef].frames[f].xoff-meflist[mef].segs[iseg].medxoff) > 10*meflist[mef].segs[iseg].sigxoff ||
	 fabsf(meflist[mef].frames[f].yoff-meflist[mef].segs[iseg].medyoff) > 10*meflist[mef].segs[iseg].sigyoff) {
	printf("%ld %ld %f %f %f %f %f %f\n", f+1, meflist[mef].frames[f].iang,
	       meflist[mef].frames[f].xoff, meflist[mef].segs[iseg].medxoff, meflist[mef].segs[iseg].sigxoff,
	       meflist[mef].frames[f].yoff, meflist[mef].segs[iseg].medyoff, meflist[mef].segs[iseg].sigyoff);
	isok = 0;
      }
#endif
    }

    if(isok) {
      nok++;

      /* It's OK, write out to the list */
      rv = fprintf(fp, "%s\n", fnlist[f]);
      if(rv <= 0) {
	report_syserr(errstr, "write");
	goto error;
      }
    }
  }

  /* Done */
  rv = fclose(fp);
  if(rv == EOF) {
    report_syserr(errstr, "close");
    goto error;
  }

  if(verbose)
    printf("Kept %ld good frames out of %ld, discarded %ld\n", nok, nf, nf - nok);

  return(0);

 error:
  return(1);
}
