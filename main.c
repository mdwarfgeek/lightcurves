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

static int write_lc (fitsfile *reff, fitsfile *fits,
		     struct buffer_info *buf, struct lc_mef *mefinfo,
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
	  "         -a        Disables aperture selection code.\n"
	  "         -d        Enables difference imaging mode.\n"
	  "         -f degree Apply polynomial of 'degree' for systematics removal.\n"
	  "         -i file   Apply intrapixel correction from 'file'.\n"
	  "         -m        Allow for meridian offset in rms and weighting.\n"
	  "         -mm       Same, only also removes it.\n"
	  "         -n        Do not renormalise median to reference magnitude.\n"
	  "         -s level  Override saturation level to 'level'.\n"
	  "         -V file   Use 'instrument version' table from 'file' (MEarth only).\n"
	  "         -u mag    Set upper(,lower) mag limit for systematics correction.\n\n"
	  "Output:\n"
	  "         -g file   Writes good frames list to 'file'.\n"
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
  int domerid = 0;
  int norenorm = 0;
  int polydeg = -1;

  int len, maxflen, fspc;
  float *medbuf1 = (float *) NULL, *medbuf2, medsat, medlim;
  long nmedsat, nmedlim, nstartot;

  long cflagmed, star;
  long *tmpmed = (long *) NULL;

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
  while((c = getopt(argc, argv, "a:df:g:i:mno:pqs:u:vV:")) != -1)
    switch(c) {
    case 'a':
      aperture = (int) strtol(optarg, &ep, 0);
      if(*ep != '\0' || aperture < 0)
	fatal(1, "invalid aperture: %s", optarg);
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

    meflist[mef].degree = polydeg;
    meflist[mef].aperture = aperture;
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
    if(read_ref(inf, &(meflist[mef]), diffmode, satlev, errstr))
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
      if(verbose && isatty(1))
	printf("\r Reading %*s (%*ld of %*ld)", maxflen, fnlist[f], fspc, f+1, fspc, nf);

      if(read_cat(fnlist[f], f, mef, &(meflist[mef]), &buf,
		  dointra, &(intralist[mef]),
		  doinstvers, instverslist, ninstvers,
		  diffmode, satlev, errstr))
	fatal(1, "read_cat: %s: %s", fnlist[f], errstr);
    }

    if(verbose && isatty(1))
      printf("\n");

    /* Sort out averages */
    meflist[mef].avsigma = sqrtf(meflist[mef].avsigma / nf);
    meflist[mef].avskyfit = sqrtf(meflist[mef].avskyfit / nf);
    meflist[mef].avapcor /= nf;
    meflist[mef].avscint /= nf;

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
    if(sysulim < 0)
      meflist[mef].sysulim = -1.0;  /* calculate it later */
    else 
      meflist[mef].sysulim = meflist[mef].zp - sysulim;  /* user-supplied */

    if(sysllim < 0)
      meflist[mef].sysllim = -1.0;  /* calculate it later */
    else 
      meflist[mef].sysllim = meflist[mef].zp - sysllim;  /* user-supplied */

    /* Call into the main part of the program */
    if(lightcurves(&buf, &(meflist[mef]), norenorm, errstr))
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

      if(write_lc(inf, outf, &buf, &(meflist[mef]), errstr))
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
		     char *errstr) {
  int status = 0, col, ncols;

  char *ttype[] = { "x", "y", "medflux", "rms", "chisq", "nchisq",
		    "class", "bflag", "cflag", "sflag", "pointer",
		    "apflux", "aprms", "apoffsets", "apradius",
		    "hjd", "flux", "fluxerr", "xlc", "ylc", "airmass", "ha",
		    "weight", "peak", "flags",
		    "ra", "dec" };
  char *tform[] = { "1E", "1E", "1E", "1E", "1E", "1J",
		    "1I", "1I", "1J", "1J", "1J",
		    "", "", "", "1E",
		    "", "", "", "", "", "", "",
		    "", "", "",
		    "1D", "1D" };
  char *tunit[] = { "pixels", "pixels", "mag", "mag", "", "",
		    "", "", "", "", "",
		    "mag", "mag", "mag", "pixels",
		    "days", "mag", "mag", "pixels", "pixels", "", "radians",
		    "", "counts", "",
		    "radians", "radians" };
  char *tdisp[] = { "F8.2", "F8.2", "F7.4", "F7.4", "F10.1", "I4",
		    "I2", "I2", "I4", "I8", "I8",
		    "F7.4", "F7.4", "F7.4", "F4.2",
		    "F14.6", "F7.4", "F7.4", "F8.2", "F8.2", "F6.4", "F9.6",
		    "F9.0", "F5.0", "I3",
		    "F9.6", "F9.6" };
  char kbuf[FLEN_KEYWORD];
  char cbuf[FLEN_COMMENT];
  char tabuf[FLEN_VALUE], tsbuf[FLEN_VALUE];
  char tfbuf[FLEN_VALUE], tdbuf[FLEN_VALUE], tbbuf[FLEN_VALUE];

  long pt, star;

  struct lc_point *lcbuf = (struct lc_point *) NULL;

  double *epos = (double *) NULL;

  float *xbuf = (float *) NULL, *ybuf, *medbuf, *rmsbuf, *chibuf, *apbuf;
  double *rabuf = (double *) NULL, *decbuf;
  long *nchibuf = (long *) NULL, *ptrbuf, *cfbuf, *sfbuf;
  short *clsbuf = (short *) NULL, *bfbuf;
  float *apmedbuf = (float *) NULL, *aprmsbuf, *apoffbuf;
  float *fluxbuf = (float *) NULL, *fluxerrbuf, *xlcbuf, *ylcbuf, *airbuf, *habuf, *wtbuf;
  float *peakbuf;
  double *hjdbuf = (double *) NULL;
  unsigned char *flagbuf = (unsigned char *) NULL;

  int ikey, nkeys, kclass;
  char card[FLEN_CARD];
  long r, frow, rblksz, soff;

  long satflag;
  unsigned char flags;

  int ap, ap1, ap2;
  int allast;

  int iseg;

  /* Generate tform specifier for fluxes and errors */
  snprintf(tabuf, sizeof(tabuf), "%dE", NFLUX);
  snprintf(tsbuf, sizeof(tabuf), "%ldE", mefinfo->nseg*NFLUX);
  snprintf(tdbuf, sizeof(tdbuf), "%ldD", mefinfo->nf);
  snprintf(tfbuf, sizeof(tfbuf), "%ldE", mefinfo->nf);
  snprintf(tbbuf, sizeof(tbbuf), "%ldB", mefinfo->nf);

  tform[11] = tabuf;
  tform[12] = tabuf;
  tform[13] = tsbuf;
  tform[15] = tdbuf;
  tform[16] = tfbuf;
  tform[17] = tfbuf;
  tform[18] = tfbuf;
  tform[19] = tfbuf;
  tform[20] = tfbuf;
  tform[21] = tfbuf;
  tform[22] = tfbuf;
  tform[23] = tfbuf;
  tform[24] = tbbuf;

  /* Create table */
  ncols = sizeof(ttype) / sizeof(ttype[0]);

  ffcrtb(fits, BINARY_TBL, mefinfo->nstars, ncols, ttype, tform, tunit, "", &status);
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
  ffpkyf(fits, "UMLIM", mefinfo->zp - mefinfo->sysulim, 4,
	 "Upper mag limit for fit", &status);
  ffpkyf(fits, "LMLIM", mefinfo->zp - mefinfo->sysllim, 4,
	 "Lower mag limit for fit", &status);
  ffpkyj(fits, "POLYDEG", mefinfo->degree,
	 "Polynomial degree in fit", &status);
  ffpkyj(fits, "APSEL", mefinfo->aperture,
	 "Aperture used (0 = automatic)", &status);
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
    ffpkyg(fits, kbuf, mefinfo->frames[pt].mjd, 7, cbuf, &status);
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

    /* Calculate Earth's heliocentric position at this MJD */
    getearth(mefinfo->mjdref + mefinfo->frames[pt].mjd, epos + 3*pt);
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
  xbuf = (float *) malloc(6 * rblksz * sizeof(float));
  rabuf = (double *) malloc(2 * rblksz * sizeof(double));
  nchibuf = (long *) malloc(4 * rblksz * sizeof(long));
  clsbuf = (short *) malloc(2 * rblksz * sizeof(short));
  apmedbuf = (float *) malloc((2+mefinfo->nseg) * rblksz * NFLUX * sizeof(float));
  fluxbuf = (float *) malloc(8 * rblksz * mefinfo->nf * sizeof(float));
  hjdbuf = (double *) malloc(rblksz * mefinfo->nf * sizeof(double));
  flagbuf = (unsigned char *) malloc(rblksz * mefinfo->nf * sizeof(unsigned char));
  if(!xbuf || !rabuf || !nchibuf || !clsbuf || !apmedbuf || !fluxbuf || !hjdbuf || !flagbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  ybuf = xbuf + rblksz;
  medbuf = xbuf + 2 * rblksz;
  rmsbuf = xbuf + 3 * rblksz;
  chibuf = xbuf + 4 * rblksz;
  apbuf = xbuf + 5 * rblksz;

  decbuf = rabuf + rblksz;

  bfbuf = clsbuf + rblksz;

  ptrbuf = nchibuf + rblksz;
  cfbuf = nchibuf + 2 * rblksz;
  sfbuf = nchibuf + 3 * rblksz;

  aprmsbuf = apmedbuf + rblksz * NFLUX;
  apoffbuf = apmedbuf + 2 * rblksz * NFLUX;

  fluxerrbuf = fluxbuf + rblksz * mefinfo->nf;
  xlcbuf = fluxbuf + 2 * rblksz * mefinfo->nf;
  ylcbuf = fluxbuf + 3 * rblksz * mefinfo->nf;
  airbuf = fluxbuf + 4 * rblksz * mefinfo->nf;
  habuf = fluxbuf + 5 * rblksz * mefinfo->nf;
  wtbuf = fluxbuf + 6 * rblksz * mefinfo->nf;
  peakbuf = fluxbuf + 7 * rblksz * mefinfo->nf;

  /* Loop through all stars - read in the lightcurve points one
   * at a time and write out blocks of 'rblksz' objects.
   */
  r = 0;
  frow = 1;

  ap1 = (mefinfo->aperture ? mefinfo->aperture-1 : 0);
  ap2 = (mefinfo->aperture ? mefinfo->aperture : NFLUX);

  for(star = 0; star < mefinfo->nstars; star++) {
    /* Fill in buffers */
    xbuf[r] = mefinfo->stars[star].x;
    ybuf[r] = mefinfo->stars[star].y;
    medbuf[r] = (mefinfo->stars[star].med > 0.0 ?
		 mefinfo->zp - mefinfo->stars[star].med : -999.0);
    rmsbuf[r] = mefinfo->stars[star].rms;
    chibuf[r] = mefinfo->stars[star].chisq;
    nchibuf[r] = mefinfo->stars[star].nchisq;
    clsbuf[r] = mefinfo->stars[star].cls;
    bfbuf[r] = mefinfo->stars[star].bflag;
    cfbuf[r] = mefinfo->stars[star].cflag;
    ptrbuf[r] = mefinfo->stars[star].ptr;

    for(ap = 0; ap < NFLUX; ap++) {
      apmedbuf[r*NFLUX+ap] = -999.0;
      aprmsbuf[r*NFLUX+ap] = -999.0;

      for(iseg = 0; iseg < mefinfo->nseg; iseg++)
	apoffbuf[(r*NFLUX+ap)*mefinfo->nseg+iseg] = -999.0;
    }

    for(ap = ap1; ap < ap2; ap++) {
      apmedbuf[r*NFLUX+ap] = (mefinfo->stars[star].medflux[ap] > 0.0 ?
			      mefinfo->zp - mefinfo->stars[star].medflux[ap] : -999.0);
      aprmsbuf[r*NFLUX+ap] = mefinfo->stars[star].sigflux[ap];

      for(iseg = 0; iseg < mefinfo->nseg; iseg++)
	apoffbuf[(r*NFLUX+ap)*mefinfo->nseg+iseg] = mefinfo->stars[star].segs[iseg].corr[ap];
    }

    apbuf[r] = mefinfo->stars[star].apradius;
    rabuf[r] = mefinfo->stars[star].ra;
    decbuf[r] = mefinfo->stars[star].dec;

    /* Get lightcurve */
    if(buffer_fetch_object(buf, lcbuf, 0, mefinfo->nf, star, 0, errstr))
      goto error;

    /* Fill in buffer */
    soff = r * mefinfo->nf;

    satflag = 0;
    for(pt = 0; pt < mefinfo->nf; pt++) {
      flags = 0;

      if(lcbuf[pt].flux != 0.0) {
	fluxbuf[soff+pt] = mefinfo->zp - lcbuf[pt].flux;

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
	  fluxerrbuf[soff+pt] = lcbuf[pt].fluxerrcom;
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
      wtbuf[soff+pt] = lcbuf[pt].wt;
      peakbuf[soff+pt] = lcbuf[pt].peak;

      flagbuf[soff+pt] = flags;

      /* Calculate HJD (as UTC) */
      hjdbuf[soff+pt] = mefinfo->mjdref + mefinfo->frames[pt].mjd +
	                hjdcorr(epos + 3*pt,
				mefinfo->stars[star].ra,
				mefinfo->stars[star].dec);
    }

    sfbuf[r] = satflag;

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
      ffpclj(fits, 10, frow, 1, r, sfbuf, &status);
      ffpclj(fits, 11, frow, 1, r, ptrbuf, &status);
      ffpcne(fits, 12, frow, 1, r * NFLUX, apmedbuf, -999.0, &status);
      ffpcne(fits, 13, frow, 1, r * NFLUX, aprmsbuf, -999.0, &status);
      ffpcne(fits, 14, frow, 1, mefinfo->nseg * r * NFLUX, apoffbuf, -999.0, &status);
      ffpcle(fits, 15, frow, 1, r, apbuf, &status);
      ffpcnd(fits, 16, frow, 1, r * mefinfo->nf, hjdbuf, -999.0, &status);
      ffpcne(fits, 17, frow, 1, r * mefinfo->nf, fluxbuf, -999.0, &status);
      ffpcne(fits, 18, frow, 1, r * mefinfo->nf, fluxerrbuf, -999.0, &status);
      ffpcne(fits, 19, frow, 1, r * mefinfo->nf, xlcbuf, -999.0, &status);
      ffpcne(fits, 20, frow, 1, r * mefinfo->nf, ylcbuf, -999.0, &status);
      ffpcne(fits, 21, frow, 1, r * mefinfo->nf, airbuf, -999.0, &status);
      ffpcne(fits, 22, frow, 1, r * mefinfo->nf, habuf, -999.0, &status);
      ffpcne(fits, 23, frow, 1, r * mefinfo->nf, wtbuf, -999.0, &status);
      ffpcne(fits, 24, frow, 1, r * mefinfo->nf, peakbuf, -999.0, &status);
      ffpclb(fits, 25, frow, 1, r * mefinfo->nf, flagbuf, &status);
      ffpcld(fits, 26, frow, 1, r, rabuf, &status);
      ffpcld(fits, 27, frow, 1, r, decbuf, &status);
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
    ffpclj(fits, 10, frow, 1, r, sfbuf, &status);
    ffpclj(fits, 11, frow, 1, r, ptrbuf, &status);
    ffpcne(fits, 12, frow, 1, r * NFLUX, apmedbuf, -999.0, &status);
    ffpcne(fits, 13, frow, 1, r * NFLUX, aprmsbuf, -999.0, &status);
    ffpcne(fits, 14, frow, 1, mefinfo->nseg * r * NFLUX, apoffbuf, -999.0, &status);
    ffpcle(fits, 15, frow, 1, r, apbuf, &status);
    ffpcnd(fits, 16, frow, 1, r * mefinfo->nf, hjdbuf, -999.0, &status);
    ffpcne(fits, 17, frow, 1, r * mefinfo->nf, fluxbuf, -999.0, &status);
    ffpcne(fits, 18, frow, 1, r * mefinfo->nf, fluxerrbuf, -999.0, &status);
    ffpcne(fits, 19, frow, 1, r * mefinfo->nf, xlcbuf, -999.0, &status);
    ffpcne(fits, 20, frow, 1, r * mefinfo->nf, ylcbuf, -999.0, &status);
    ffpcne(fits, 21, frow, 1, r * mefinfo->nf, airbuf, -999.0, &status);
    ffpcne(fits, 22, frow, 1, r * mefinfo->nf, habuf, -999.0, &status);
    ffpcne(fits, 23, frow, 1, r * mefinfo->nf, wtbuf, -999.0, &status);
    ffpcne(fits, 24, frow, 1, r * mefinfo->nf, peakbuf, -999.0, &status);
    ffpclb(fits, 25, frow, 1, r * mefinfo->nf, flagbuf, &status);
    ffpcld(fits, 26, frow, 1, r, rabuf, &status);
    ffpcld(fits, 27, frow, 1, r, decbuf, &status);
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
  free((void *) rabuf);
  rabuf = (double *) NULL;
  free((void *) nchibuf);
  nchibuf = (long *) NULL;
  free((void *) clsbuf);
  clsbuf = (short *) NULL;
  free((void *) apmedbuf);
  apmedbuf = (float *) NULL;
  free((void *) fluxbuf);
  fluxbuf = (float *) NULL;
  free((void *) hjdbuf);
  hjdbuf = (double *) NULL;
  free((void *) flagbuf);
  flagbuf = (unsigned char *) NULL;

  return(0);

 error:
  if(lcbuf)
    free((void *) lcbuf);
  if(epos)
    free((void *) epos);
  if(xbuf)
    free((void *) xbuf);
  if(rabuf)
    free((void *) rabuf);
  if(nchibuf)
    free((void *) nchibuf);
  if(clsbuf)
    free((void *) clsbuf);
  if(apmedbuf)
    free((void *) apmedbuf);
  if(fluxbuf)
    free((void *) fluxbuf);
  if(hjdbuf)
    free((void *) hjdbuf);
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
