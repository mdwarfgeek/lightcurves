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

#include "cvtunit.h"
#include "floatmath.h"
#include "util.h"

static int read_ref (fitsfile *fits, struct lc_mef *mefinfo, char *errstr);
static int read_cat (char *catfile, int iframe, int mef, struct lc_mef *mefinfo,
		     struct buffer_info *buf,
		     int dointra, struct intra *icorr, char *errstr);

static int write_lc (fitsfile *reff, fitsfile *fits,
		     struct buffer_info *buf, struct lc_mef *mefinfo, char *errstr);
static int write_goodlist (char *outfile, struct lc_mef *meflist, int nmefs,
			   char **fnlist, char *errstr);

static char *sstrip (char *str);

/* Flux keywords and aperture sizes in terms of rcore */
static char *flux_keys[NFLUX] = { "Core_flux",
				  "Core2_flux",
				  "Core3_flux" };
float flux_apers[NFLUX] = { 1.0,
			    M_SQRT2,
			    2.0 };
static char *apcor_keys[NFLUX] = { "APCOR",
				   "APCOR2",
				   "APCOR3" };

/* Getopt stuff */
extern char *optarg;
extern int optind;

/* Verbose flag */
int verbose = 1;

static int noplots = 0;
static int diffmode = 0;
static float sysbodge = 0.0;

static void usage (char *av) {
  fprintf(stderr, "Usage:\t%s [options] reffile file [...]\n\n", av);
  fprintf(stderr,
	  "Lightcurve processing:\n"
	  "         -a        Disables aperture selection code.\n"
	  "         -b bfac   Add 'bfac' mmag in quadrature with errors.\n"
	  "         -d        Enables difference imaging mode.\n"
	  "         -f degree Apply polynomial of 'degree' for systematics removal.\n"
	  "         -i file   Apply intrapixel correction from 'file'.\n"
	  "         -u mag    Set upper mag limit for systematics correction.\n\n"
	  "Output:\n"
	  "         -g file   Writes good frames list to 'file'.\n"
	  "         -o file   Writes lightcurves to 'file'.\n"
	  "         -p        Disables plots.\n"
	  "         -q        Decreases the verbosity level of the program.\n"
  	  "         -v        Increases the verbosity level of the program.\n");
  exit(1);
}  

int main (int argc, char *argv[]) {
  char *pn = (char *) NULL, *avzero, *ep;
  int c;

  char errstr[ERRSTR_LEN], line[16384];

  char *refname, **fnlist = (char **) NULL, **fnp;
  struct lc_mef *meflist = (struct lc_mef *) NULL;
  struct intra *intralist = (struct intra *) NULL;
  struct buffer_info buf;

  long a, f, nf, op;
  char *p;

  FILE *fp;

  fitsfile *inf, *outf;
  int status = 0, ext, mef, nmefs;

  char outfile[FLEN_FILENAME-1], fnbuf[FLEN_FILENAME];
  int dooutput = 0;

  char goodfile[FLEN_FILENAME];
  int dogood = 0;

  char intrafile[FLEN_FILENAME];
  int dointra = 0;

  int noapsel = 0;
  int polydeg = -1;
  int pcafit = 0;

  int len, maxflen, fspc;
  float *medbuf1 = (float *) NULL, *medbuf2, medsat, medlim;
  long nmedsat, nmedlim, nstartot;

  long cflagmed, star;
  long *tmpmed = (long *) NULL;

  float syslim = -1.0;

  /* Set the program name for error reporting */
  if(argv[0])
    pn = basename(argv[0]);
  
  if(!pn)
    pn = "lightcurves";

  setprogname(pn);

  avzero = argv[0];

  /* Extract command-line arguments */
  while((c = getopt(argc, argv, "ab:df:g:i:o:pqtu:v")) != -1)
    switch(c) {
    case 'a':
      noapsel++;
      break;
    case 'b':
      sysbodge = (float) strtod(optarg, &ep);
      if(*ep != '\0' || sysbodge < 0)
	fatal(1, "invalid systematics bodge value: %s", optarg);

      sysbodge /= 1000;  /* convert to mag */

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
      outfile[sizeof(goodfile)-1] = '\0';
      dogood = 1;
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
    case 'p':
      noplots++;
      break;
    case 'q':
      verbose--;
      break;
    case 't':
      pcafit++;
      break;
    case 'u':
      syslim = (float) strtod(optarg, &ep);
      if(*ep != '\0' || syslim < 0)
	fatal(1, "invalid syslim value: %s", optarg);
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

  op = 0;
  nf = 0;
  for(a = 1; a < argc; a++) {
    if(*(argv[a]) == '@') {
      /* @list form */
      fp = fopen(argv[a] + 1, "r");
      if(!fp)
	error(1, "open: %s", argv[a] + 1);
      
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

    meflist[mef].degree = polydeg;

    meflist[mef].avsigma = 0.0;
    meflist[mef].avapcor = 0.0;

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
    if(read_ref(inf, &(meflist[mef]), errstr))
      fatal(1, "read_ref: HDU %d: %s", mef+2, errstr);

    /* Get disk buffer */
    if(buffer_alloc(&buf, meflist[mef].nstars, nf, NFLUX, errstr))
      fatal(1, "buffer_alloc: %s", errstr);

    /* Allocate buffer for frame info */
    meflist[mef].frames = (struct lc_frame *) malloc(nf * sizeof(struct lc_frame));
    if(!meflist[mef].frames)
      error(1, "malloc");

    /* Read frames */
    for(f = 0; f < nf; f++) {
      if(verbose)
	printf("\r Reading %*s (%*ld of %*ld)", maxflen, fnlist[f], fspc, f+1, fspc, nf);

      if(read_cat(fnlist[f], f, mef, &(meflist[mef]), &buf,
		  dointra, &(intralist[mef]), errstr))
	fatal(1, "read_cat: %s: %s", fnlist[f], errstr);
    }

    if(verbose)
      printf("\n");

    /* Sort out averages */
    meflist[mef].avsigma /= nf;
    meflist[mef].avapcor /= nf;

    /* Fix the cflag column - sometimes in difference imaging there
     * are frames with zero confidence all-over, so subtract off
     * the median cflag value.
     */
    tmpmed = (long *) malloc(meflist[mef].nstars * sizeof(long));
    if(!tmpmed)
      error(1, "malloc");

    for(star = 0; star < meflist[mef].nstars; star++)
      tmpmed[star] = meflist[mef].stars[star].cflag;

    sortlong(tmpmed, meflist[mef].nstars);
    cflagmed = meflist[mef].nstars % 2 ?
               tmpmed[meflist[mef].nstars/2] :
               (tmpmed[meflist[mef].nstars/2-1] +  tmpmed[meflist[mef].nstars/2]) / 2;

    free((void *) tmpmed);
    tmpmed = (long *) NULL;

    for(star = 0; star < meflist[mef].nstars; star++) {
      meflist[mef].stars[star].cflag -= cflagmed;

      if(meflist[mef].stars[star].cflag < 0)
	meflist[mef].stars[star].cflag = 0;
    }

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
	     "  Median cflag:         %ld\n",
	     meflist[mef].refflim, cflagmed);

    }

    if(meflist[mef].satmag != -999.0) {
      medbuf1[nmedsat] = meflist[mef].satmag;
      nmedsat++;
    }

    medbuf2[nmedlim] = meflist[mef].refflim;
    nmedlim++;

    nstartot += meflist[mef].nstars;

    /* Calculate syslim for this frame */
    if(syslim < 0)
      meflist[mef].syslim = -1.0;  /* calculate it later */
    else 
      meflist[mef].syslim = meflist[mef].zp - syslim;  /* user-supplied */

    /* Call into the main part of the program */
    if(lightcurves(&buf, &(meflist[mef]), noapsel, pcafit, errstr))
      fatal(1, "%s", errstr);

    /* Write out lightcurves for this MEF if requested */
    if(dooutput) {
      if(verbose)
	printf(" Writing %s\n", outfile);

      if(write_lc(inf, outf, &buf, &(meflist[mef]), errstr))
	fatal(1, "write_lc: %s", errstr);
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
    if(do_plots(meflist, nmefs, medsat, medlim, sysbodge, errstr))
      fatal(1, "do_plots: %s");
  }

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

#define RADECZP(x1, y1, ra, dec) {						\
  float x, y, xi, xn, rv, rfac;							\
  float aa, alpha, delta;							\
										\
  x = x1 - c;									\
  y = y1 - f;									\
										\
  xi = a * x + b * y;								\
  xn = d * x + e * y;								\
  rv = sqrtf(xi * xi + xn * xn);						\
										\
  rfac = projp1 + projp3 * rv * rv;   /* NB this is only a 1st order approx */	\
  rv /= rfac;									\
										\
  rfac = projp1 + projp3 * rv * rv;   /* now 2nd order correction		\
                                     * accurate to few 100ths of pixel */	\
  xi /= rfac;									\
  xn /= rfac;									\
										\
  aa    = atanf(xi * secd / (1.0 - xn * tand));					\
  alpha = aa + tpa;								\
										\
  if(xi != 0.0)									\
    delta = atanf((xn + tand) * sinf(aa) / (xi * secd));			\
  else										\
    delta = atanf((xn + tand) / (1.0 - xn * tand));				\
										\
  if(alpha > tpi)								\
    alpha -= tpi;								\
										\
  if(alpha < 0.0)								\
    alpha += tpi;								\
										\
  (ra)  = alpha;								\
  (dec) = delta;								\
}

static int read_ref (fitsfile *fits, struct lc_mef *mefinfo, char *errstr) {
  int status = 0;

  char *colnames[5] = { "X_coordinate", "Y_coordinate", "Peak_height",
			"Classification", "Areal_7_profile" };
  int gcols[5+NFLUX], col, collim;

  struct lc_star *stars = (struct lc_star *) NULL;

  float *xbuf = (float *) NULL, *ybuf, *fluxbuf, *pkhtbuf, *clsbuf, *a7buf;
  float *sattmp = (float *) NULL;
  long nsattmp;

  float tpa, tpd, a, b, c, d, e, f, projp1, projp3, secd, tand;
  float skylev, skynoise, satlev, exptime, rcore, gain, magzpt, percorr;
  float skyvar, tpi;
  float apcor[NFLUX];

  long nrows, rblksz, roff, rout, remain, rread, r;
  float satflux;

  char filter[FLEN_VALUE];
  float airmass = 1.0, extinct = 0.0;
  int l1, l2, i, ilim;

  int noexp = 0;

  struct {
    char *filt;
    float extinct;
  } default_extinct_tab[] = {
    { "r ",          0.09 },
    { "g ",          0.19 },
    { "U ",          0.46 },
    { "i ",          0.05 },
    { "z ",          0.05 },
    { "B ",          0.22 },
    { "V ",          0.12 },
    { "R ",          0.08 },
    { "I ",          0.04 },
    { "stromgren u", 0.51 },
    { "stromgren v", 0.26 },
    { "stromgren b", 0.15 },
    { "stromgren y", 0.10 }
  };

  /* Read number of rows */
  ffgnrw(fits, &nrows, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgnrw");
    goto error;
  }
  
  /* Get column numbers */
  collim = sizeof(colnames) / sizeof(colnames[0]);
  for(col = 0; col < collim; col++) {
    ffgcno(fits, CASEINSEN, colnames[col], &(gcols[col]), &status);
    if(status == COL_NOT_UNIQUE)
      status = 0;  /* ignore */
    else if(status) {
      fitsio_err(errstr, status, "ffgcno: %s", colnames[col]);
      goto error;
    }
  }

  for(col = 0; col < NFLUX; col++) {
    ffgcno(fits, CASEINSEN, flux_keys[col], &(gcols[collim+col]), &status);
    if(status == COL_NOT_UNIQUE)
      status = 0;  /* ignore */
    else if(status) {
      fitsio_err(errstr, status, "ffgcno: %s", flux_keys[col]);
      goto error;
    }
  }

  /* Read WCS info */
  ffgkye(fits, "CRVAL1", &tpa, (char *) NULL, &status);
  ffgkye(fits, "CRVAL2", &tpd, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkye: CRVAL[12]");
    goto error;
  }

  ffgkye(fits, "CRPIX1", &c, (char *) NULL, &status);
  ffgkye(fits, "CRPIX2", &f, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkye: CRPIX[12]");
    goto error;
  }

  ffgkye(fits, "CD1_1", &a, (char *) NULL, &status);
  ffgkye(fits, "CD1_2", &b, (char *) NULL, &status);
  ffgkye(fits, "CD2_1", &d, (char *) NULL, &status);
  ffgkye(fits, "CD2_2", &e, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkye: CD[12]_[12]");
    goto error;
  }

  ffgkye(fits, "PROJP1", &projp1, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    projp1 = 1.0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: PROJP1");
    goto error;
  }

  ffgkye(fits, "PROJP3", &projp3, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    projp3 = 220.0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: PROJP3");
    goto error;
  }

  tpa *= DEG_TO_RAD;
  tpd *= DEG_TO_RAD;

  a *= DEG_TO_RAD;
  b *= DEG_TO_RAD;
  d *= DEG_TO_RAD;
  e *= DEG_TO_RAD;

  secd = 1.0 / cosf(tpd);
  tand = tanf(tpd);

  /* Get saturation level */
  ffgkye(fits, "SATURATE", &satlev, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    satlev = 65535;

    if(verbose)
      printf("Warning: using default satlev = %.1f\n", satlev);
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: SATURATE");
    goto error;
  }

  /* Read keywords for photometry */
  ffgkye(fits, "EXPTIME", &exptime, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkye(fits, "EXPOSED", &exptime, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkye(fits, "EXP_TIME", &exptime, (char *) NULL, &status);
      if(status) {
	fitsio_err(errstr, status, "ffgkye: EXP_TIME");
	goto error;
      }
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: EXPOSED");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: EXPTIME");
    goto error;
  }

  exptime = fabsf(exptime);
  if(exptime < 1.0)
    exptime = 1.0;

  ffgkye(fits, "SKYLEVEL", &skylev, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkye: SKYLEVEL");
    goto error;
  }

  ffgkye(fits, "SKYNOISE", &skynoise, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkye: SKYNOISE");
    goto error;
  }

  ffgkye(fits, "RCORE", &rcore, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkye: RCORE");
    goto error;
  }

  ffgkye(fits, "GAIN", &gain, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    gain = 1.0;  /* !!! */
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: GAIN");
    goto error;
  }

  ffgkye(fits, "MAGZPT", &magzpt, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkye(fits, "ZMAG", &magzpt, (char *) NULL, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgkye: ZMAG");
      goto error;
    }
    else {
      noexp = 1;  /* don't add in 2.5log10(exptime) */
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: MAGZPT");
    goto error;
  }

  for(col = 0; col < NFLUX; col++) {
    ffgkye(fits, apcor_keys[col], &(apcor[col]), (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      apcor[col] = 1.0;
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: %s", apcor_keys[col]);
      goto error;
    }
    else {
      apcor[col] = powf(10.0, 0.4 * apcor[col]);
    }
  }

  ffgkye(fits, "PERCORR", &percorr, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    percorr = 1.0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: PERCORR");
    goto error;
  }
  else {
    percorr = powf(10.0, 0.4 * percorr);
  }

  /* Get airmass */
  ffgkye(fits, "AIRMASS", &airmass, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkye(fits, "AMSTART", &airmass, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      airmass = 1.0;
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: AMSTART");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: AIRMASS");
    goto error;
  }

  /* Get filter name in case we can't get extinction */
  ffgkys(fits, "WFFBAND", filter, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkys(fits, "FILTER", filter, (char *) NULL, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgkye: FILTER");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: WFFBAND");
    goto error;
  }

  /* Append a space */
  l1 = strlen(filter);
  if(l1+1 < sizeof(filter)) {
    filter[l1] = ' ';
    filter[l1+1] = '\0';
    l1++;
  }

  /* Attempt to get extinction */
  ffgkye(fits, "EXTINCT", &extinct, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    extinct = 0.0;

    /* Attempt to find it in the table of defaults */
    ilim = sizeof(default_extinct_tab) / sizeof(default_extinct_tab[0]);

    for(i = 0; i < ilim; i++) {
      l2 = strlen(default_extinct_tab[i].filt);

      if(l1 >= l2 && !strncmp(filter, default_extinct_tab[i].filt, l2)) {
	/* Found it */
	extinct = default_extinct_tab[i].extinct;
	break;
      }
    }
  }

  if(noexp)
    mefinfo->zp = magzpt;
  else
    mefinfo->zp = magzpt + 2.5 * log10f(exptime) - (airmass - 1.0)*extinct;

  skyvar = M_PI * rcore * rcore * skynoise * skynoise;
  tpi = 2.0 * M_PI;

  /* Get block size for row I/O */
  ffgrsz(fits, &rblksz, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgrsz");
    goto error;
  }
  
  /* Allocate column buffers */
  xbuf = (float *) malloc(6 * rblksz * sizeof(float));
  if(!xbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }
  
  ybuf = xbuf + rblksz;
  fluxbuf = xbuf + 2 * rblksz;
  pkhtbuf = xbuf + 3 * rblksz;
  clsbuf = xbuf + 4 * rblksz;
  a7buf = xbuf + 5 * rblksz;
  
  /* Allocate memory for catalogue stars */
  stars = (struct lc_star *) malloc(nrows * sizeof(struct lc_star));
  sattmp = (float *) malloc(nrows * sizeof(float));
  if(!stars || !sattmp) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Read catalogue */
  roff = 0L;
  remain = nrows;
  rout = 0L;

  satflux = 0.0;

  nsattmp = 0;
  while(remain > 0) {
    rread = (remain > rblksz ? rblksz : remain);
    
    ffgcve(fits, gcols[0], roff + 1L, 1L, rread, 0.0, xbuf, (int *) NULL, &status);
    ffgcve(fits, gcols[1], roff + 1L, 1L, rread, 0.0, ybuf, (int *) NULL, &status);
    ffgcve(fits, gcols[2], roff + 1L, 1L, rread, 0.0, pkhtbuf, (int *) NULL, &status);
    ffgcve(fits, gcols[3], roff + 1L, 1L, rread, 0.0, clsbuf, (int *) NULL, &status);
    ffgcve(fits, gcols[4], roff + 1L, 1L, rread, 0.0, a7buf, (int *) NULL, &status);

    for(col = 0; col < NFLUX; col++) {
      ffgcve(fits, gcols[5+col], roff + 1L, 1L, rread, 0.0, fluxbuf,
	     (int *) NULL, &status);

      for(r = 0; r < rread; r++) {
	rout = roff + r;

	stars[rout].ref[col].flux = fluxbuf[r] * apcor[col] * percorr;
	stars[rout].ref[col].fluxerr = fabsf(fluxbuf[r]) * apcor[col] / gain;
	/* +skyvar * flux_apers[col] * flux_apers[col] ? */

	if(pkhtbuf[r] > satlev)
	  stars[rout].ref[col].satur = 1;
	else
	  stars[rout].ref[col].satur = 0;
      }
    }

    if(status) {
      fitsio_err(errstr, status, "ffgcve");
      goto error;
    }
    
    for(r = 0; r < rread; r++) {
      rout = roff + r;

      stars[rout].ptr = rout + 1;
      stars[rout].x = xbuf[r];
      stars[rout].y = ybuf[r];

      RADECZP(xbuf[r], ybuf[r], stars[rout].ra, stars[rout].dec);

      stars[rout].cls = NINT(clsbuf[r]);
      stars[rout].bflag = (a7buf[r] < 0.0 ? 1 : 0);
      stars[rout].cflag = 0;

      if(pkhtbuf[r] > satlev) {
	sattmp[nsattmp] = stars[rout].ref[0].flux;
	nsattmp++;
      }
    }
    
    roff += rread;
    remain -= rread;
  }

  /* Free workspace */
  free((void *) xbuf);
  xbuf = (float *) NULL;

  /* Determine saturation level robustly - 10%ile */
  if(nsattmp > 0) {
    sortfloat(sattmp, nsattmp);
    satflux = sattmp[nsattmp/4];
  }

  mefinfo->stars = stars;
  mefinfo->nstars = nrows;
  if(satflux > 0.0)
    mefinfo->satmag = mefinfo->zp - 2.5 * log10f(satflux);
  else
    mefinfo->satmag = -999.0;
  mefinfo->refexp = exptime;
  mefinfo->refsigma = skynoise;
  mefinfo->refflim = mefinfo->zp - 2.5 * log10f(5.0 * sqrtf(M_PI * rcore * rcore) *
						skynoise * apcor[0]);

  mefinfo->refgain = gain;
  mefinfo->refrcore = rcore;

  memcpy(&(mefinfo->apcor), apcor, sizeof(apcor));
  mefinfo->percorr = percorr;

  free((void *) sattmp);
  sattmp = (float *) NULL;

  return(0);

 error:
  if(stars)
    free((void *) stars);
  if(xbuf)
    free((void *) xbuf);
  if(sattmp)
    free((void *) sattmp);

  return(1);
}

static int read_cat (char *catfile, int iframe, int mef, struct lc_mef *mefinfo,
		     struct buffer_info *buf,
		     int dointra, struct intra *icorr, char *errstr) {
  fitsfile *fits;
  int status = 0;

  char *colnames[4] = { "X_coordinate", "Y_coordinate", "Peak_height", "Skyrms" };
  char *optcolnames[1] = { "Bad_pixels" };
  int gcols[5+NFLUX], col, collim, optcollim;

  struct lc_point *points = (struct lc_point *) NULL;

  float *xbuf = (float *) NULL, *ybuf, *fluxbuf, *pkhtbuf, *skyrmsbuf, *badpixbuf;

  float skylev, skynoise, satlev, exptime, rcore, gain, percorr;
  float skyvar, tpi, tmp, expfac;
  double mjd;

  float apcor[NFLUX];

  long nrows, rblksz, roff, remain, rout, rread, r;
  float flux, fluxerr;

  /* Open catalogue */
  ffopen(&fits, catfile, READONLY, &status);
  if(status) {
    fitsio_err(errstr, status, "ffopen: %s", catfile);
    goto error;
  }

  /* Move to the right MEF */
  ffmahd(fits, mef+2, (int *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffmahd: %s: HDU %d", catfile, mef+2);
    goto error;
  }

  /* Read number of rows */
  ffgnrw(fits, &nrows, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgnrw");
    goto error;
  }
  
  /* Get column numbers */
  collim = sizeof(colnames) / sizeof(colnames[0]);
  for(col = 0; col < collim; col++) {
    ffgcno(fits, CASEINSEN, colnames[col], &(gcols[col]), &status);
    if(status == COL_NOT_UNIQUE)
      status = 0;  /* ignore */
    else if(status) {
      fitsio_err(errstr, status, "ffgcno: %s", colnames[col]);
      goto error;
    }
  }

  optcollim = sizeof(optcolnames) / sizeof(optcolnames[0]);
  for(col = 0; col < optcollim; col++) {
    ffgcno(fits, CASEINSEN, optcolnames[col], &(gcols[collim+col]), &status);
    if(status == COL_NOT_UNIQUE)
      status = 0;  /* ignore */
    else if(status == COL_NOT_FOUND) {
      status = 0;
      gcols[collim+col] = 0;
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgcno: %s", optcolnames[col]);
      goto error;
    }
  }

  collim += optcollim;

  for(col = 0; col < NFLUX; col++) {
    ffgcno(fits, CASEINSEN, flux_keys[col], &(gcols[collim+col]), &status);
    if(status == COL_NOT_UNIQUE)
      status = 0;  /* ignore */
    else if(status) {
      fitsio_err(errstr, status, "ffgcno: %s", flux_keys[col]);
      goto error;
    }
  }

  /* Get saturation level */
  ffgkye(fits, "SATURATE", &satlev, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    satlev = 65535;

    if(verbose)
      printf("Warning: using default satlev = %.1f\n", satlev);
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: SATURATE");
    goto error;
  }

  /* Read keywords for photometry */
  ffgkye(fits, "EXPTIME", &exptime, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkye(fits, "EXPOSED", &exptime, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkye(fits, "EXP_TIME", &exptime, (char *) NULL, &status);
      if(status) {
	fitsio_err(errstr, status, "ffgkye: EXP_TIME");
	goto error;
      }
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: EXPOSED");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: EXPTIME");
    goto error;
  }

  exptime = fabsf(exptime);
  if(exptime < 1.0)
    exptime = 1.0;

  expfac = mefinfo->refexp / exptime;

  ffgkye(fits, "SKYLEVEL", &skylev, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkye: SKYLEVEL");
    goto error;
  }

  ffgkye(fits, "SKYNOISE", &skynoise, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkye: SKYNOISE");
    goto error;
  }

  ffgkye(fits, "RCORE", &rcore, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkye: RCORE");
    goto error;
  }

  ffgkye(fits, "GAIN", &gain, (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkye: GAIN");
    goto error;
  }

  skyvar = M_PI * rcore * rcore * skynoise * skynoise;
  tpi = 2.0 * M_PI;

  for(col = 0; col < NFLUX; col++) {
    ffgkye(fits, apcor_keys[col], &(apcor[col]), (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      apcor[col] = mefinfo->apcor[col];  /* as a backup */
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: %s", apcor_keys[col]);
      goto error;
    }
    else {
      apcor[col] = powf(10.0, 0.4 * apcor[col]);
    }
  }

  ffgkye(fits, "PERCORR", &percorr, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    percorr = mefinfo->percorr;  /* as a backup */
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: PERCORR");
    goto error;
  }
  else {
    percorr = powf(10.0, 0.4 * percorr);
  }

  ffgkyd(fits, "MJD-OBS", &mjd, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkyd(fits, "MJD", &mjd, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkyd(fits, "JD", &mjd, (char *) NULL, &status);
      if(status) {
	fitsio_err(errstr, status, "ffgkyd: JD");
	goto error;
      }
      else
	mjd -= 2400000.5;
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkyd: MJD");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkyd: MJD-OBS");
    goto error;
  }

  /* Get block size for row I/O */
  ffgrsz(fits, &rblksz, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgrsz");
    goto error;
  }
  
  /* Allocate column buffers */
  xbuf = (float *) malloc(6 * rblksz * sizeof(float));
  if(!xbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }
  
  ybuf = xbuf + rblksz;
  fluxbuf = xbuf + 2 * rblksz;
  pkhtbuf = xbuf + 3 * rblksz;
  skyrmsbuf = xbuf + 4 * rblksz;
  badpixbuf = xbuf + 5 * rblksz;

  /* Allocate memory for lightcurve points */
  points = (struct lc_point *) malloc(rblksz * sizeof(struct lc_point));
  if(!points) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Read catalogue */
  roff = 0L;
  remain = nrows;
  rout = 0L;

  while(remain > 0) {
    rread = (remain > rblksz ? rblksz : remain);
    
    ffgcve(fits, gcols[0], roff + 1L, 1L, rread, 0.0, xbuf, (int *) NULL, &status);
    ffgcve(fits, gcols[1], roff + 1L, 1L, rread, 0.0, ybuf, (int *) NULL, &status);
    ffgcve(fits, gcols[2], roff + 1L, 1L, rread, 0.0, pkhtbuf, (int *) NULL, &status);
    ffgcve(fits, gcols[3], roff + 1L, 1L, rread, 0.0, skyrmsbuf, (int *) NULL, &status);

    if(gcols[4])
      ffgcve(fits, gcols[4], roff + 1L, 1L, rread, 0.0, badpixbuf, (int *) NULL, &status);

    for(col = 0; col < NFLUX; col++) {
      ffgcve(fits, gcols[5+col], roff + 1L, 1L, rread, 0.0, fluxbuf,
	     (int *) NULL, &status);

      for(r = 0; r < rread; r++) {
	rout = roff + r;

	points[r].x = xbuf[r];
	points[r].y = ybuf[r];

	if(diffmode) {
	  flux = fluxbuf[r] * mefinfo->apcor[col] * mefinfo->percorr;
	  fluxerr = fabsf(fluxbuf[r]) * mefinfo->apcor[col] / gain +
	    (skyvar + skyrmsbuf[r]*skyrmsbuf[r]) * flux_apers[col] * flux_apers[col];

	  if(flux == 0.0 || mefinfo->stars[rout].ref[col].flux == 0.0)
	    flux = 0.0;
	  else
	    flux += mefinfo->stars[rout].ref[col].flux;

	  fluxerr += mefinfo->stars[rout].ref[col].fluxerr;

	  points[r].satur = mefinfo->stars[rout].ref[col].satur;
	}
	else {
	  flux = fluxbuf[r] * apcor[col] * percorr;
	  fluxerr = fabsf(fluxbuf[r]) * apcor[col] / gain +
	    (skyvar + skyrmsbuf[r]*skyrmsbuf[r]) * flux_apers[col] * flux_apers[col];

	  if(pkhtbuf[r] > satlev || mefinfo->stars[rout].ref[col].satur)
	    points[r].satur = 1;
	  else
	    points[r].satur = 0;
	}

	points[r].flux = 2.5 * log10f(MAX(1.0, flux));

	if(flux > 0.0) {
	  points[r].fluxerr = 2.5 * log10f(1.0 + sqrtf(fluxerr) / flux);

	  if(sysbodge > 0.0)
	    points[r].fluxerr = sqrtf(points[r].fluxerr*points[r].fluxerr +
				      sysbodge*sysbodge);
	}
	else
	  points[r].fluxerr = 0.0;

	/* Apply intrapixel correction if requested */
	if(dointra)
	  points[r].flux += calc_intra(xbuf[r], ybuf[r], icorr);

	/* Accumulate counts of frames with bad pixels */
	if(gcols[4] && badpixbuf[r] > 0.0)
	  mefinfo->stars[rout].cflag++;
      }

      /* Write out those */
      if(buffer_put_frame(buf, points, roff, rread, iframe, col, errstr))
	goto error;
    }

    if(status) {
      fitsio_err(errstr, status, "ffgcve");
      goto error;
    }
    
    roff += rread;
    remain -= rread;
  }

  /* Free workspace */
  free((void *) xbuf);
  xbuf = (float *) NULL;

  free((void *) points);
  points = (struct lc_point *) NULL;

  /* Accumulate average sigma */
  tmp = skynoise * sqrtf(expfac);

  if(diffmode)
    mefinfo->avsigma += sqrtf(tmp * tmp + mefinfo->refsigma * mefinfo->refsigma);
  else
    mefinfo->avsigma += tmp;

  mefinfo->avapcor += apcor[0];

  /* Store this frame MJD */
  mefinfo->frames[iframe].mjd = mjd;

  /* Make sure rcore is correctly inserted into aperture radius */
  for(r = 0; r < nrows; r++)
    mefinfo->stars[r].apradius = rcore;

  /* Close file */
  ffclos(fits, &status);
  if(status) {
    fitsio_err(errstr, status, "ffclos");
    goto error;
  }

  return(0);

 error:
  if(xbuf)
    free((void *) xbuf);
  if(points)
    free((void *) points);

  return(1);
}

static int write_lc (fitsfile *reff, fitsfile *fits,
		     struct buffer_info *buf, struct lc_mef *mefinfo, char *errstr) {
  int status = 0, col, ncols;

  char *ttype[] = { "x", "y", "medflux", "rms", "chisq", "nchisq",
		    "class", "bflag", "cflag", "pointer",
		    "apradius",
		    "hjd", "flux", "fluxerr", "xlc", "ylc", "ra", "dec" };
  char *tform[] = { "1E", "1E", "1E", "1E", "1E", "1J",
		    "1I", "1I", "1J", "1J",
		    "1E",
		    "", "", "", "", "", "1E", "1E" };
  char *tunit[] = { "pixels", "pixels", "mag", "mag", "", "",
		    "", "", "", "",
		    "pixels",
		    "days", "mag", "mag", "pixels", "pixels", "radians", "radians" };
  char *tdisp[] = { "F8.2", "F8.2", "F7.4", "F7.4", "F10.1", "I4",
		    "I2", "I2", "I4", "I8",
		    "F4.2",
		    "F14.6", "F7.4", "F7.4", "F8.2", "F8.2", "F9.6", "F9.6" };
  char kbuf[FLEN_KEYWORD], tfbuf[FLEN_VALUE], tdbuf[FLEN_VALUE], cbuf[FLEN_COMMENT];

  long pt, star;

  struct lc_point *lcbuf = (struct lc_point *) NULL;

  double *epos = (double *) NULL;

  float *xbuf = (float *) NULL, *ybuf, *medbuf, *rmsbuf, *chibuf, *apbuf, *rabuf, *decbuf;
  long *nchibuf = (long *) NULL, *ptrbuf, *cfbuf;
  short *clsbuf = (short *) NULL, *bfbuf;
  float *fluxbuf = (float *) NULL, *fluxerrbuf, *xlcbuf, *ylcbuf;
  double *hjdbuf = (double *) NULL;

  int ikey, nkeys, kclass;
  char card[FLEN_CARD];
  long r, frow, rblksz, soff;

  /* Generate tform specifier for fluxes and errors */
  snprintf(tdbuf, sizeof(tdbuf), "%ldD", mefinfo->nf);
  snprintf(tfbuf, sizeof(tfbuf), "%ldE", mefinfo->nf);

  tform[11] = tdbuf;
  tform[12] = tfbuf;
  tform[13] = tfbuf;
  tform[14] = tfbuf;
  tform[15] = tfbuf;

  /* Create table */
  ncols = sizeof(ttype) / sizeof(ttype[0]);

  ffcrtb(fits, BINARY_TBL, mefinfo->nstars, ncols, ttype, tform, tunit, "", &status);
  if(status) {
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

  /* Allocate buffer for earth positions */
  epos = (double *) malloc(3 * mefinfo->nf * sizeof(double));
  if(!epos) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Write out frame information */
  ffpkyj(fits, "NMEAS", mefinfo->nf,
	 "Number of points in each lightcurve", &status);
  ffpkyg(fits, "MJDBASE", mefinfo->mjdref, 7,
	 "Base MJD for time axis", &status);
  ffpkyf(fits, "SATMAG", mefinfo->satmag, 4,
	 "Approximate saturation magnitude", &status);
  ffpkyf(fits, "FLIM", mefinfo->refflim, 4,
	 "Flux limit of reference catalogue", &status);
  ffpkyf(fits, "ZP", mefinfo->zp, 4,
	 "Zeropoint for magnitudes", &status);
  if(status) {
    fitsio_err(errstr, status, "ffkpy: frame info");
    goto error;
  }

  for(pt = 0; pt < mefinfo->nf; pt++) {
    snprintf(kbuf, sizeof(kbuf), "TV%ld", pt+1);
    snprintf(cbuf, sizeof(cbuf), "Time value for datapoint %ld", pt+1);
    ffpkyg(fits, kbuf, mefinfo->frames[pt].mjd, 7, cbuf, &status);
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
    snprintf(cbuf, sizeof(cbuf), "Frame extinction for datapoint %ld", pt+1);
    ffpkyf(fits, kbuf, mefinfo->frames[pt].extinc, 4, cbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpkyf: %s", kbuf);
      goto error;
    }

    /* Calculate Earth's heliocentric position at this MJD */
    getearth(mefinfo->mjdref + mefinfo->frames[pt].mjd, epos + 3*pt);
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
  xbuf = (float *) malloc(8 * rblksz * sizeof(float));
  nchibuf = (long *) malloc(3 * rblksz * sizeof(long));
  clsbuf = (short *) malloc(2 * rblksz * sizeof(short));
  fluxbuf = (float *) malloc(4 * rblksz * mefinfo->nf * sizeof(float));
  hjdbuf = (double *) malloc(rblksz * mefinfo->nf * sizeof(double));
  if(!xbuf || !nchibuf || !clsbuf || !fluxbuf || !hjdbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  ybuf = xbuf + rblksz;
  medbuf = xbuf + 2 * rblksz;
  rmsbuf = xbuf + 3 * rblksz;
  chibuf = xbuf + 4 * rblksz;
  rabuf = xbuf + 5 * rblksz;
  decbuf = xbuf + 6 * rblksz;
  apbuf = xbuf + 7 * rblksz;

  bfbuf = clsbuf + rblksz;

  ptrbuf = nchibuf + rblksz;
  cfbuf = nchibuf + 2 * rblksz;

  fluxerrbuf = fluxbuf + rblksz * mefinfo->nf;
  xlcbuf = fluxbuf + 2 * rblksz * mefinfo->nf;
  ylcbuf = fluxbuf + 3 * rblksz * mefinfo->nf;

  /* Loop through all stars - read in the lightcurve points one
   * at a time and write out blocks of 'rblksz' objects.
   */
  r = 0;
  frow = 1;

  for(star = 0; star < mefinfo->nstars; star++) {
    /* Fill in buffers */
    xbuf[r] = mefinfo->stars[star].x;
    ybuf[r] = mefinfo->stars[star].y;
    medbuf[r] = (mefinfo->stars[star].medflux[0] > 0.0 ?
		 mefinfo->zp - mefinfo->stars[star].medflux[0] : -999.0);
    rmsbuf[r] = mefinfo->stars[star].rms;
    chibuf[r] = mefinfo->stars[star].chisq;
    nchibuf[r] = mefinfo->stars[star].nchisq;
    clsbuf[r] = mefinfo->stars[star].cls;
    bfbuf[r] = mefinfo->stars[star].bflag;
    cfbuf[r] = mefinfo->stars[star].cflag;
    ptrbuf[r] = mefinfo->stars[star].ptr;
    apbuf[r] = mefinfo->stars[star].apradius;
    rabuf[r] = mefinfo->stars[star].ra;
    decbuf[r] = mefinfo->stars[star].dec;

    /* Get lightcurve */
    if(buffer_fetch_object(buf, lcbuf, 0, mefinfo->nf, star, 0, errstr))
      goto error;

    /* Fill in buffer */
    soff = r * mefinfo->nf;

    for(pt = 0; pt < mefinfo->nf; pt++) {
      if(lcbuf[pt].flux != 0.0 && !lcbuf[pt].satur) {
	fluxbuf[soff+pt] = mefinfo->zp - lcbuf[pt].flux;

	if(lcbuf[pt].fluxerr > 0.0)
	  fluxerrbuf[soff+pt] = lcbuf[pt].fluxerr;
	else
	  fluxerrbuf[soff+pt] = -999.0;
      }
      else
	fluxbuf[soff+pt] = -999.0;

      xlcbuf[soff+pt] = lcbuf[pt].x;
      ylcbuf[soff+pt] = lcbuf[pt].y;

      /* Calculate HJD */
      hjdbuf[soff+pt] = mefinfo->mjdref + mefinfo->frames[pt].mjd +
	                hjdcorr(epos + 3*pt,
				mefinfo->stars[star].ra,
				mefinfo->stars[star].dec);
    }

    r++;

    if(r >= rblksz) {
      /* Flush */
      ffpcle(fits, 1, frow, 1, r, xbuf, &status);
      ffpcle(fits, 2, frow, 1, r, ybuf, &status);
      ffpcne(fits, 3, frow, 1, r, medbuf, -999.0, &status);
      ffpcle(fits, 4, frow, 1, r, rmsbuf, &status);
      ffpcle(fits, 5, frow, 1, r, chibuf, &status);
      ffpclj(fits, 6, frow, 1, r, nchibuf, &status);
      ffpcli(fits, 7, frow, 1, r, clsbuf, &status);
      ffpcli(fits, 8, frow, 1, r, bfbuf, &status);
      ffpclj(fits, 9, frow, 1, r, cfbuf, &status);
      ffpclj(fits, 10, frow, 1, r, ptrbuf, &status);
      ffpcle(fits, 11, frow, 1, r, apbuf, &status);
      ffpcnd(fits, 12, frow, 1, r * mefinfo->nf, hjdbuf, -999.0, &status);
      ffpcne(fits, 13, frow, 1, r * mefinfo->nf, fluxbuf, -999.0, &status);
      ffpcne(fits, 14, frow, 1, r * mefinfo->nf, fluxerrbuf, -999.0, &status);
      ffpcne(fits, 15, frow, 1, r * mefinfo->nf, xlcbuf, -999.0, &status);
      ffpcne(fits, 16, frow, 1, r * mefinfo->nf, ylcbuf, -999.0, &status);
      ffpcle(fits, 17, frow, 1, r, rabuf, &status);
      ffpcle(fits, 18, frow, 1, r, decbuf, &status);
      if(status) {
	fitsio_err(errstr, status, "ffpcl");
	goto error;
      }

      frow += r;
      r = 0;
    }
  }

  /* Flush out buffers */
  if(r > 0) {
    /* Flush */
    ffpcle(fits, 1, frow, 1, r, xbuf, &status);
    ffpcle(fits, 2, frow, 1, r, ybuf, &status);
    ffpcne(fits, 3, frow, 1, r, medbuf, -999.0, &status);
    ffpcle(fits, 4, frow, 1, r, rmsbuf, &status);
    ffpcle(fits, 5, frow, 1, r, chibuf, &status);
    ffpclj(fits, 6, frow, 1, r, nchibuf, &status);
    ffpcli(fits, 7, frow, 1, r, clsbuf, &status);
    ffpcli(fits, 8, frow, 1, r, bfbuf, &status);
    ffpclj(fits, 9, frow, 1, r, cfbuf, &status);
    ffpclj(fits, 10, frow, 1, r, ptrbuf, &status);
    ffpcle(fits, 11, frow, 1, r, apbuf, &status);
    ffpcnd(fits, 12, frow, 1, r * mefinfo->nf, hjdbuf, -999.0, &status);
    ffpcne(fits, 13, frow, 1, r * mefinfo->nf, fluxbuf, -999.0, &status);
    ffpcne(fits, 14, frow, 1, r * mefinfo->nf, fluxerrbuf, -999.0, &status);
    ffpcne(fits, 15, frow, 1, r * mefinfo->nf, xlcbuf, -999.0, &status);
    ffpcne(fits, 16, frow, 1, r * mefinfo->nf, ylcbuf, -999.0, &status);
    ffpcle(fits, 17, frow, 1, r, rabuf, &status);
    ffpcle(fits, 18, frow, 1, r, decbuf, &status);
    if(status) {
      fitsio_err(errstr, status, "ffpcl");
      goto error;
    }
  }    

  free((void *) lcbuf);
  lcbuf = (struct lc_point *) NULL;
  free((void *) epos);
  epos = (double *) NULL;
  free((void *) xbuf);
  xbuf = (float *) NULL;
  free((void *) nchibuf);
  nchibuf = (long *) NULL;
  free((void *) clsbuf);
  clsbuf = (short *) NULL;
  free((void *) fluxbuf);
  fluxbuf = (float *) NULL;
  free((void *) hjdbuf);
  hjdbuf = (double *) NULL;

  return(0);

 error:
  if(lcbuf)
    free((void *) lcbuf);
  if(epos)
    free((void *) epos);
  if(xbuf)
    free((void *) xbuf);
  if(nchibuf)
    free((void *) nchibuf);
  if(clsbuf)
    free((void *) clsbuf);
  if(fluxbuf)
    free((void *) fluxbuf);
  if(hjdbuf)
    free((void *) hjdbuf);

  return(1);
}

static int write_goodlist (char *outfile, struct lc_mef *meflist, int nmefs,
			   char **fnlist, char *errstr) {
  FILE *fp;
  long f, nf;
  int rv, mef, isok;

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
    /* Check each MEF for frame offset amplitude or rms > 0.05 */
    isok = 1;

    for(mef = 0; mef < nmefs; mef++) {
      if(fabsf(meflist[mef].frames[f].offset) > 0.05 ||
	 fabsf(meflist[mef].frames[f].rms) > 0.05)
	isok = 0;
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

static char *sstrip (char *str) {
  char *p;

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

