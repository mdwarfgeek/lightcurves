#include <sys/types.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#include "lightcurves.h"

#include "cvtunit.h"
#include "fitsutil.h"
#include "util.h"

/* Flux keywords and aperture sizes in terms of rcore.  These have no declared
   size so NFLUX can be changed without having to keep modifying this file. */
static char *flux_keys_32[] = { "Core_flux",
				"Core2_flux",
				"Core3_flux",
				"Core4_flux" };
static char *flux_keys_80[] = { "Aper_flux_3",
				"Aper_flux_4",
				"Aper_flux_5",
				"Aper_flux_6" };
float flux_apers[] = { 1.0,
		       M_SQRT2,
		       2.0,
		       2.0*M_SQRT2 };
static char *apcor_keys_32[] = { "APCOR",
				 "APCOR2",
				 "APCOR3",
				 "APCOR4" };
static char *apcor_keys_80[] = { "APCOR3",
				 "APCOR4",
				 "APCOR5",
				 "APCOR6" };

static struct {
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
  { "Ic",          0.05 },
  { "stromgren u", 0.51 },
  { "stromgren v", 0.26 },
  { "stromgren b", 0.15 },
  { "stromgren y", 0.10 },
  { "i+z",         0.10 },  /* MEarth */
  { "I_Burke",     0.05 }
};

/* Reconstruct internal state from a lightcurve file.  This routine
 * should work in the normal case, but we still need the reference
 * fluxes for difference imaging - NOT YET IMPLEMENTED.
 */

int read_lc (fitsfile *fits, struct lc_mef *mefinfo,
	     char *errstr) {
  int status = 0, anynull;

  char *colnames[] = { "x", "y", "class", "pointer", "bflag", "cflag",
		       "medflux", "rms", "offsets", "apnum", "apradius", "compok",
		       "ra", "dec", "pmra", "pmdec", "refmag" };
  unsigned char tpap[] = { 0, 0, 0, 0, 0, 0,
			   1, 1, 1, 0, 0, 0,
			   0, 0, 0, 0, 0 };

  char cn[FLEN_VALUE+1];
  int *gcols = (int *) NULL;
  int tcol, icol, ntmpl, ncols;

  int ap, ap1, ap2, napcol;

  int cats_are_80 = 0;

  float exptime, pedestal = 0, skylev, skynoise, rcore, gain;
  float magzpt;
  long numzpt;
  int haveimgzpt;
  float apcor[NFLUX], percorr;
  char filter[FLEN_VALUE];
  float airmass = 1.0, extinct = 0.0;
  int l1, l2, i, ilim;

  int noexp = 0;

  double mjd;

  float *apbuf = (float *) NULL, *refmagbuf, *pmabuf, *pmdbuf;
  double *xbuf = (double *) NULL, *ybuf, *rabuf, *decbuf;
  short *clsbuf = (short *) NULL, *bfbuf, *apnumbuf;
  long *ptrbuf = (long *) NULL, *cfbuf;
  unsigned char *compokbuf = (unsigned char *) NULL;
  float *apmedbuf = (float *) NULL, *aprmsbuf, *apoffbuf;

  struct lc_star *stars = (struct lc_star *) NULL;
  long nmeas, nrowmast;

  float umlim, lmlim, apcor7;
  long degree;

  long nrows, r, rr, roff, remain, rread, rblksz;

  int iseg;
  long iver, jtmp;
  char kbuf[FLEN_KEYWORD];

  /* Read number of rows */
  ffgnrw(fits, &nrows, &status);
  if(status) {
    fitsio_err(errstr, status, "could not get table dimensions");
    goto error;
  }

  /* Get header information */
  ffgkyj(fits, "NMEAS", &nmeas, (char *) NULL, &status);
  ffgkyj(fits, "NROWMAST", &nrowmast, (char *) NULL, &status);
  ffgkyd(fits, "MJDBASE", &(mefinfo->mjdref), (char *) NULL, &status);
  ffgkye(fits, "SATMAG", &(mefinfo->satmag), (char *) NULL, &status);
  ffgkye(fits, "FLIM", &(mefinfo->refflim), (char *) NULL, &status);
  ffgkye(fits, "ZP", &(mefinfo->zp), (char *) NULL, &status);
  ffgkye(fits, "UMLIM", &umlim, (char *) NULL, &status);
  ffgkyj(fits, "POLYDEG", &degree, (char *) NULL, &status);
  ffgkyj(fits, "APSEL", &(mefinfo->aperture), (char *) NULL, &status);
  ffgkyj(fits, "APMODE", &(mefinfo->apselmode), (char *) NULL, &status);
  ffgkyj(fits, "DOMERID", &(mefinfo->domerid), (char *) NULL, &status);
  ffgkyj(fits, "NSEGME", &(mefinfo->nseg), (char *) NULL, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgkyj: NMEAS");
    goto error;
  }

  /* Allocate and read segment information */
  mefinfo->segs = (struct lc_segment *) malloc(mefinfo->nseg * sizeof(struct lc_segment));
  if(!mefinfo->segs) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  for(iseg = 0; iseg < mefinfo->nseg; iseg++) {
    snprintf(kbuf, sizeof(kbuf), "SEGV%d", iseg+1);
    ffgkyj(fits, kbuf, &iver, (char *) NULL, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgkyj: %s", kbuf);
      goto error;
    }

    if(iver > 0) {
      mefinfo->segs[iseg].instvers = (struct instvers *) malloc(sizeof(struct instvers));
      if(!mefinfo->segs[iseg].instvers) {
	report_syserr(errstr, "malloc");
	goto error;
      }

      mefinfo->segs[iseg].instvers->iver = iver;

      snprintf(kbuf, sizeof(kbuf), "SEGD%d", iseg+1);
      ffgkyj(fits, kbuf, &(mefinfo->segs[iseg].instvers->date), (char *) NULL, &status);
      if(status) {
	fitsio_err(errstr, status, "ffgkyj: %s", kbuf);
	goto error;
      }
    }
    else
      mefinfo->segs[iseg].instvers = (struct instvers *) NULL;

    snprintf(kbuf, sizeof(kbuf), "SEGA%d", iseg+1);
    ffgkyj(fits, kbuf, &jtmp, (char *) NULL, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgkyj: %s", kbuf);
      goto error;
    }

    mefinfo->segs[iseg].iang = jtmp;
  }

  /* New and thus optional */
  ffgkye(fits, "LMLIM", &lmlim, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    lmlim = umlim+USEMAG;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: LMLIM");
    goto error;
  }

  ffgkyl(fits, "THEOSKY", &(mefinfo->theosky), (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    mefinfo->theosky = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkyl: THEOSKY");
    goto error;
  }

  /* Try for field angle - old files will not have it, in which case
   * fallback to the old, incorrect, HA-based method */
  ffgkye(fits, "REFFANG", &(mefinfo->reffang), (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    mefinfo->havefang = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: REFFANG");
    goto error;
  }
  else
    mefinfo->havefang = 1;

  mefinfo->sysulim = mefinfo->zp - umlim;
  mefinfo->sysllim = mefinfo->zp - lmlim;
  mefinfo->degree = degree;

  /* Simple test for 80-column catalogue */
  ffgkye(fits, "APCOR7", &apcor7, (char *) NULL, &status);
  if(status == KEY_NO_EXIST)
    status = 0;
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: APCOR7");
    goto error;
  }
  else
    cats_are_80 = 1;

  /* Read keywords for photometry */
  ffgkye(fits, "EXPTIME", &exptime, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkye(fits, "EXPOSED", &exptime, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkye(fits, "EXP_TIME", &exptime, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
        status = 0;
        exptime = 0;

        if(verbose) {
          printf("Warning: no exposure time found\n");
          mefinfo->warned |= WARNED_EXPTIME;
        }
      }
      else if(status) {
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

  ffgkye(fits, "PEDESTAL", &pedestal, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    pedestal = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: PEDESTAL");
    goto error;
  }

  ffgkye(fits, "SKYLEVEL", &skylev, (char *) NULL, &status);
  if(status == KEY_NO_EXIST && nrows == 0) {
    status = 0;
    skylev = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: SKYLEVEL");
    goto error;
  }

  ffgkye(fits, "SKYNOISE", &skynoise, (char *) NULL, &status);
  if(status == KEY_NO_EXIST && nrows == 0) {
    status = 0;
    skynoise = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: SKYNOISE");
    goto error;
  }

  ffgkye(fits, "RCORE", &rcore, (char *) NULL, &status);
  if(status == KEY_NO_EXIST && nrows == 0) {
    status = 0;
    rcore = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: RCORE");
    goto error;
  }

  ffgkye(fits, "GAIN", &gain, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkye(fits, "HIERARCH ESO DET OUT1 GAIN", &gain, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      ffgkye(fits, "EGAIN", &gain, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
	status = 0;
	gain = 1.0;  /* !!! */

        if(verbose) {
          printf("Warning: using default gain = %.1f\n", gain);
          mefinfo->warned |= WARNED_GAIN;
        }
      }
      else if(status) {
	fitsio_err(errstr, status, "ffgkye: EGAIN");
	goto error;
      }
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: HIERARCH ESO DET OUT1 GAIN");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: GAIN");
    goto error;
  }

  /* How many standards in image? */
  ffgkyj(fits, "IMNZPT", &numzpt, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    numzpt = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkyj: IMNZPT");
    goto error;
  }

  if(numzpt > 0) {
    /* Try to get image ZP solution */
    ffgkye(fits, "IMGZPT", &magzpt, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      haveimgzpt = 0;
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkyj: IMGZPT");
      goto error;
    }
    else if(magzpt == 99.0)  /* flag value */
      haveimgzpt = 0;
    else
      haveimgzpt = 1;
  }
  else
    haveimgzpt = 0;
  
  if(!haveimgzpt) {
    ffgkye(fits, "MAGZPT", &magzpt, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkye(fits, "ZMAG", &magzpt, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
        status = 0;
        magzpt = 25.0;
        
        if(verbose) {
          printf("Warning: using default magzpt = %.1f\n", magzpt);
          mefinfo->warned |= WARNED_MAGZPT;
        }
      }
      else if(status) {
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
  }

  for(icol = 0; icol < NFLUX; icol++) {
    ffgkye(fits, cats_are_80 ? apcor_keys_80[icol] : apcor_keys_32[icol],
	   &(apcor[icol]), (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      apcor[icol] = 1.0;
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: %s", 
		 cats_are_80 ? apcor_keys_80[icol] : apcor_keys_32[icol]);
      goto error;
    }
    else {
      apcor[icol] = powf(10.0, 0.4 * apcor[icol]);
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

  /* Get filter name for plots */
  ffgkys(fits, "FILTER", filter, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkys(fits, "WFFBAND", filter, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkys(fits, "HIERARCH ESO INS FILT1 NAME", filter, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
	status = 0;
	ffgkys(fits, "FILTER2", filter, (char *) NULL, &status);
	if(status == KEY_NO_EXIST) {
	  status = 0;
	  ffgkys(fits, "INSFILTE", filter, (char *) NULL, &status);
          if(status == KEY_NO_EXIST) {
            status = 0;
            filter[0] = '\0';
          }
	  else if(status) {
	    fitsio_err(errstr, status, "ffgkye: INSFILTE");
	    goto error;
	  }
	}
	else if(status) {
	  fitsio_err(errstr, status, "ffgkye: FILTER2");
	  goto error;
	}
      }
      else if(status) {
	fitsio_err(errstr, status, "ffgkye: HIERARCH ESO INS FILT1 NAME");
	goto error;
      }
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: WFFBAND");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: FILTER");
    goto error;
  }

  /* Copy */
  strncpy(mefinfo->filter, filter, sizeof(mefinfo->filter)-1);
  mefinfo->filter[sizeof(mefinfo->filter)-1] = '\0';

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

  ffgkyd(fits, "MJD-OBS", &mjd, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkyd(fits, "MJD", &mjd, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkyd(fits, "JD", &mjd, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
        status = 0;
        mjd = J2K;

        if(verbose) {
          printf("Warning: no MJDs available, timestamps will be garbage\n");
          mefinfo->warned |= WARNED_TIME;
        }
      }
      else if(status) {
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

  /* Correct to midpoint of observation */
  mjd += 0.5 * exptime / 86400.0;

  if(mefinfo->aperture) {
    ap1 = mefinfo->aperture-1;
    ap2 = mefinfo->aperture-1;
    napcol = 1;
  }
  else {
    ap1 = 0;
    ap2 = NFLUX;
    napcol = ap2-ap1 + 1;
  }

  /* Decide how many columns */
  ntmpl = sizeof(colnames) / sizeof(colnames[0]);

  ncols = 0;
  for(tcol = 0; tcol < ntmpl; tcol++) {
    ncols++;

    if(tpap[tcol])
      ncols += ap2-ap1;
  }

  /* Allocate arrays */
  gcols = (int *) malloc(ncols * sizeof(int));
  if(!gcols) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Get column numbers */
  icol = 0;
  for(tcol = 0; tcol < ntmpl; tcol++) {
    ffgcno(fits, CASEINSEN, colnames[tcol], gcols+icol, &status);
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
	ffgcno(fits, CASEINSEN, cn, gcols+icol, &status);
	if(status == COL_NOT_UNIQUE)
	  status = 0;  /* ignore */
	else if(status) {
	  fitsio_err(errstr, status, "ffgcno: %s", cn);
	  goto error;
	}

	icol++;
      }
  }
  
  /* Get block size for row I/O */
  ffgrsz(fits, &rblksz, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgrsz");
    goto error;
  }
  
  /* Allocate column buffers */
  apbuf = (float *) malloc(4 * rblksz * sizeof(float));
  xbuf = (double *) malloc(4 * rblksz * sizeof(double));
  clsbuf = (short *) malloc(3 * rblksz * sizeof(short));
  ptrbuf = (long *) malloc(2 * rblksz * sizeof(long));
  compokbuf = (unsigned char *) malloc(rblksz * sizeof(unsigned char));
  apmedbuf = (float *) malloc((2+mefinfo->nseg) * rblksz * napcol * sizeof(float));
  if(!apbuf || !xbuf || !clsbuf || !ptrbuf || !compokbuf || !apmedbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  refmagbuf = apbuf + rblksz;
  pmabuf = apbuf + 2 * rblksz;
  pmdbuf = apbuf + 3 * rblksz;

  ybuf = xbuf + rblksz;
  rabuf = xbuf + 2 * rblksz;
  decbuf = xbuf + 3 * rblksz;

  bfbuf = clsbuf + rblksz;
  apnumbuf = clsbuf + 2 * rblksz;

  cfbuf = ptrbuf + rblksz;

  aprmsbuf = apmedbuf + napcol * rblksz;
  apoffbuf = apmedbuf + 2 * napcol * rblksz;

  /* Allocate memory for catalogue stars */
  if(nrows > 0) {
    stars = (struct lc_star *) malloc(nrows * sizeof(struct lc_star));
    if(!stars) {
      report_syserr(errstr, "malloc");
      goto error;
    }
  }
  else
    stars = (struct lc_star *) NULL;

  for(rr = 0; rr < nrows; rr++) {
    stars[rr].segs = (struct lc_star_segment *) malloc(mefinfo->nseg * sizeof(struct lc_star_segment));
    if(!stars[rr].segs) {
      report_syserr(errstr, "malloc");
      goto error;
    }
  }

  /* Read catalogue */
  roff = 0L;
  remain = nrows;

  while(remain > 0) {
    rread = (remain > rblksz ? rblksz : remain);
    
    icol = 0;
    ffgcvd(fits, gcols[icol++], roff + 1, 1, rread, -999.0, xbuf, &anynull, &status);
    ffgcvd(fits, gcols[icol++], roff + 1, 1, rread, -999.0, ybuf, &anynull, &status);
    ffgcvi(fits, gcols[icol++], roff + 1, 1, rread, 0, clsbuf, &anynull, &status);
    ffgcvj(fits, gcols[icol++], roff + 1, 1, rread, 0, ptrbuf, &anynull, &status);
    ffgcvi(fits, gcols[icol++], roff + 1, 1, rread, 0, bfbuf, &anynull, &status);
    ffgcvj(fits, gcols[icol++], roff + 1, 1, rread, 0, cfbuf, &anynull, &status);

    for(ap = 0; ap < napcol; ap++)
      ffgcve(fits, gcols[icol++], roff + 1, 1, rread, -999.0, apmedbuf+ap*rblksz, &anynull,
	     &status);

    for(ap = 0; ap < napcol; ap++)
      ffgcve(fits, gcols[icol++], roff + 1, 1, rread, -999.0, aprmsbuf+ap*rblksz, &anynull,
	     &status);

    for(ap = 0; ap < napcol; ap++)
      ffgcve(fits, gcols[icol++], roff + 1, 1, mefinfo->nseg * rread, -999.0, apoffbuf+ap*mefinfo->nseg*rblksz, &anynull,
	     &status);

    ffgcvi(fits, gcols[icol++], roff + 1, 1, rread, 0, apnumbuf, &anynull, &status);
    ffgcve(fits, gcols[icol++], roff + 1, 1, rread, -999.0, apbuf, &anynull, &status);
    ffgcvb(fits, gcols[icol++], roff + 1, 1, rread, 0, compokbuf, &anynull, &status);
    ffgcvd(fits, gcols[icol++], roff + 1, 1, rread, -999.0, rabuf, &anynull,
	   &status);
    ffgcvd(fits, gcols[icol++], roff + 1, 1, rread, -999.0, decbuf, &anynull,
	   &status);
    ffgcve(fits, gcols[icol++], roff + 1, 1, rread, -999.0, pmabuf, &anynull, &status);
    ffgcve(fits, gcols[icol++], roff + 1, 1, rread, -999.0, pmdbuf, &anynull, &status);
    ffgcve(fits, gcols[icol++], roff + 1, 1, rread, -999.0, refmagbuf, &anynull, &status);
    if(status) {
      fitsio_err(errstr, status, "ffgcv");
      goto error;
    }
    
    for(r = 0; r < rread; r++) {
      rr = roff + r;

      stars[rr].ptr = ptrbuf[r];
      stars[rr].x = xbuf[r];
      stars[rr].y = ybuf[r];
      stars[rr].ra = rabuf[r];
      stars[rr].dec = decbuf[r];
      stars[rr].cls = clsbuf[r];
      stars[rr].bflag = bfbuf[r];
      stars[rr].cflag = cfbuf[r];

      if(pmabuf[r] == -999.0 && pmdbuf[r] == -999.0) {
        stars[rr].pmra = 0;
        stars[rr].pmdec = 0;
        stars[rr].havepm = 0;
      }
      else {
        stars[rr].pmra = pmabuf[r];
        stars[rr].pmdec = pmdbuf[r];
        stars[rr].havepm = 1;
      }

      /* UTC epoch is used rather than the correct TDB, not likely to matter. */
      source_star(&(stars[rr].src), 
                  stars[rr].ra, stars[rr].dec,
                  stars[rr].pmra, stars[rr].pmdec,
                  0.0, 0.0,
                  2000.0 + (mjd-J2K)/JYR);

      if(mefinfo->aperture) {
	stars[rr].medflux[mefinfo->aperture-1] = (apmedbuf[r] > 0.0 ? mefinfo->zp - apmedbuf[r] : -999.0);
	stars[rr].sigflux[mefinfo->aperture-1] = aprmsbuf[r];

	for(iseg = 0; iseg < mefinfo->nseg; iseg++)
          stars[rr].segs[iseg].corr[mefinfo->aperture-1] = apoffbuf[r*mefinfo->nseg+iseg];
      }
      else {
	for(ap = ap1; ap < ap2; ap++) {
	  stars[rr].medflux[ap] = (apmedbuf[(ap-ap1+1)*rblksz+r] > 0.0 ? mefinfo->zp - apmedbuf[(ap-ap1+1)*rblksz+r] : -999.0);
	  stars[rr].sigflux[ap] = aprmsbuf[(ap-ap1+1)*rblksz+r];
	  for(iseg = 0; iseg < mefinfo->nseg; iseg++)
	    stars[rr].segs[iseg].corr[ap] = apoffbuf[((ap-ap1+1)*rblksz+r)*mefinfo->nseg+iseg];
	}
      }

      stars[rr].iap = apnumbuf[r]-1;
      stars[rr].apradius = apbuf[r];
      stars[rr].compok = compokbuf[r];
      stars[rr].refmag = (refmagbuf[r] > 0.0 ? mefinfo->zp - refmagbuf[r] : -999.0);
    }
    
    roff += rread;
    remain -= rread;
  }

  mefinfo->stars = stars;
  mefinfo->nstars = nrows;
  mefinfo->nrows = nrowmast;  /* so we can check for mismatches */

  mefinfo->refexp = exptime;
  mefinfo->refsigma = mefinfo->theosky ? sqrtf(skylev / gain) : skynoise;
  mefinfo->refextinct = noexp ? 0.0 : extinct;
  mefinfo->refairmass = airmass;
  mefinfo->refmagzpt = magzpt;
  mefinfo->refgain = gain;
  mefinfo->refrcore = rcore;

  memcpy(&(mefinfo->apcor), apcor, sizeof(apcor));
  mefinfo->percorr = percorr;

  /* Free workspace */
  free((void *) gcols);
  gcols = (int *) NULL;
  free((void *) apbuf);
  apbuf = (float *) NULL;
  free((void *) xbuf);
  xbuf = (double *) NULL;
  free((void *) clsbuf);
  clsbuf = (short *) NULL;
  free((void *) ptrbuf);
  ptrbuf = (long *) NULL;
  free((void *) compokbuf);
  compokbuf = (unsigned char *) NULL;
  free((void *) apmedbuf);
  apmedbuf = (float *) NULL;

  return(0);

 error:
  if(gcols)
    free((void *) gcols);
  if(apbuf)
    free((void *) apbuf);
  if(xbuf)
    free((void *) xbuf);
  if(clsbuf)
    free((void *) clsbuf);
  if(ptrbuf)
    free((void *) ptrbuf);
  if(compokbuf)
    free((void *) compokbuf);
  if(apmedbuf)
    free((void *) apmedbuf);
  if(stars)
    free((void *) stars);

  return(1);
}

int read_ref (fitsfile *fits, struct lc_mef *mefinfo,
	      int diffmode, float satlev,
	      float sysllim, float sysulim,
	      int outcls, int wantoutcls,
	      char *errstr) {
  int status = 0;

  char *colnames[6] = { "X_coordinate", "Y_coordinate", "Peak_height",
			"Classification", "Areal_7_profile", "Skylev" };
  char *optcolnames[2] = { "PMRA", "PMDec" };
  int gcols[8+NFLUX], col, collim, optcollim;

  struct lc_star *stars = (struct lc_star *) NULL;

  double *xbuf = (double *) NULL, *ybuf;
  float *allfluxbuf = (float *) NULL, *pkhtbuf, *clsbuf, *a7buf, *locskybuf, *fluxbuf;
  float *pmabuf, *pmdbuf;
  float *sattmp = (float *) NULL;
  long nsattmp;

  struct wcs_info wcs;
  double fang;

  float skylev, pedestal = 0, skynoise, exptime, rcore, gain, percorr;
  float magzpt;
  long numzpt;
  int haveimgzpt;
  float apcor[NFLUX];

  long nrows, nstars, rblksz, roff, rin, rout, rrin, remain, rread;
  float satflux;

  char filter[FLEN_VALUE];
  float airmass = 1.0, extinct = 0.0;
  int l1, l2, i, ilim;

  double mjd;

  int noexp = 0;

  char inst[FLEN_VALUE], tel[FLEN_VALUE];
  float scatcoeff = 0.0;
  double v[3], dvdt[3], xi, xn, rr;

  int cats_are_80 = 0;

  float tmp;

  /* Read number of rows */
  ffgnrw(fits, &nrows, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgnrw");
    goto error;
  }
  
  /* Simple test for the new 80-column format */
  ffgcno(fits, CASEINSEN, flux_keys_80[0], &cats_are_80, &status);
  if(status == COL_NOT_FOUND || status == COL_NOT_UNIQUE)
    status = 0;
  else if(status) {
    fitsio_err(errstr, status, "ffgcno: %s", flux_keys_80[0]);
    goto error;
  }
  else {
    cats_are_80 = 1;

    colnames[5] = "Sky_level";  /* WHY were these changed? */
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
    ffgcno(fits, CASEINSEN, cats_are_80 ? flux_keys_80[col] : flux_keys_32[col],
	   &(gcols[collim+col]), &status);
    if(status == COL_NOT_UNIQUE)
      status = 0;  /* ignore */
    else if(status) {
      fitsio_err(errstr, status, "ffgcno: %s",
		 cats_are_80 ? flux_keys_80[col] : flux_keys_32[col]);
      goto error;
    }
  }

  /* Read WCS info */
  if(read_wcs(fits, &wcs, 0, errstr))
    goto error;

  fang = atan2(wcs.b, wcs.a);

  if(satlev < 0) {
    /* Get saturation level - tries SATLEV first, on MEarth this
     * is a better estimate, based on the non-linearity curve.
     */
    ffgkye(fits, "SATLEV", &satlev, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkye(fits, "SATURATE", &satlev, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
	status = 0;
	satlev = 65535;
	
	if(verbose > 1) {
	  printf("Warning: using default satlev = %.1f\n", satlev);
          mefinfo->warned |= WARNED_SATLEV;
        }
      }
      else if(status) {
	fitsio_err(errstr, status, "ffgkye: SATURATE");
	goto error;
      }
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: SATLEV");
      goto error;
    }
  }

  /* Read keywords for photometry */
  ffgkye(fits, "EXPTIME", &exptime, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkye(fits, "EXPOSED", &exptime, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkye(fits, "EXP_TIME", &exptime, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
        status = 0;
        exptime = 0;

        if(verbose) {
          printf("Warning: no exposure time found\n");
          mefinfo->warned |= WARNED_EXPTIME;
        }
      }
      else if(status) {
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

  ffgkye(fits, "PEDESTAL", &pedestal, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    pedestal = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: PEDESTAL");
    goto error;
  }

  ffgkye(fits, "SKYLEVEL", &skylev, (char *) NULL, &status);
  if(status == KEY_NO_EXIST && nrows == 0) {
    status = 0;
    skylev = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: SKYLEVEL");
    goto error;
  }

  ffgkye(fits, "SKYNOISE", &skynoise, (char *) NULL, &status);
  if(status == KEY_NO_EXIST && nrows == 0) {
    status = 0;
    skynoise = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: SKYNOISE");
    goto error;
  }

  ffgkye(fits, "RCORE", &rcore, (char *) NULL, &status);
  if(status == KEY_NO_EXIST && nrows == 0) {
    status = 0;
    rcore = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: RCORE");
    goto error;
  }

  ffgkye(fits, "GAIN", &gain, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkye(fits, "HIERARCH ESO DET OUT1 GAIN", &gain, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkye(fits, "EGAIN", &gain, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
	status = 0;
	gain = 1.0;  /* !!! */

        if(verbose) {
          printf("Warning: using default gain = %.1f\n", gain);
          mefinfo->warned |= WARNED_GAIN;
        }
      }
      else if(status) {
	fitsio_err(errstr, status, "ffgkye: EGAIN");
	goto error;
      }
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: HIERARCH ESO DET OUT1 GAIN");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: GAIN");
    goto error;
  }

  /* How many standards in image? */
  ffgkyj(fits, "IMNZPT", &numzpt, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    numzpt = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkyj: IMNZPT");
    goto error;
  }

  if(numzpt > 0) {
    /* Try to get image ZP solution */
    ffgkye(fits, "IMGZPT", &magzpt, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      haveimgzpt = 0;
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkyj: IMGZPT");
      goto error;
    }
    else if(magzpt == 99.0)  /* flag value */
      haveimgzpt = 0;
    else
      haveimgzpt = 1;
  }
  else
    haveimgzpt = 0;
  
  if(!haveimgzpt) {
    ffgkye(fits, "MAGZPT", &magzpt, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkye(fits, "ZMAG", &magzpt, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
        status = 0;
        magzpt = 25.0;
        
        if(verbose) {
          printf("Warning: using default magzpt = %.1f\n", magzpt);
          mefinfo->warned |= WARNED_MAGZPT;
        }
      }
      else if(status) {
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
  }

  for(col = 0; col < NFLUX; col++) {
    ffgkye(fits, cats_are_80 ? apcor_keys_80[col] : apcor_keys_32[col],
	   &(apcor[col]), (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      apcor[col] = 1.0;
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: %s", 
		 cats_are_80 ? apcor_keys_80[col] : apcor_keys_32[col]);
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
  ffgkys(fits, "FILTER", filter, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkys(fits, "WFFBAND", filter, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkys(fits, "HIERARCH ESO INS FILT1 NAME", filter, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
	status = 0;
	ffgkys(fits, "FILTER2", filter, (char *) NULL, &status);
	if(status == KEY_NO_EXIST) {
	  status = 0;
	  ffgkys(fits, "INSFILTE", filter, (char *) NULL, &status);
          if(status == KEY_NO_EXIST) {
            status = 0;
            filter[0] = '\0';
          }
	  else if(status) {
	    fitsio_err(errstr, status, "ffgkye: INSFILTE");
	    goto error;
	  }
	}
	else if(status) {
	  fitsio_err(errstr, status, "ffgkye: FILTER2");
	  goto error;
	}
      }
      else if(status) {
	fitsio_err(errstr, status, "ffgkye: HIERARCH ESO INS FILT1 NAME");
	goto error;
      }
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: WFFBAND");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: FILTER");
    goto error;
  }

  /* Copy */
  strncpy(mefinfo->filter, filter, sizeof(mefinfo->filter)-1);
  mefinfo->filter[sizeof(mefinfo->filter)-1] = '\0';

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
  else {
    mefinfo->zp = magzpt - (airmass - 1.0)*extinct;
    if(exptime)
      mefinfo->zp += 2.5 * log10f(exptime);
  }

  /* Calculate syslim for this frame */
  if(sysulim < 0)
    mefinfo->sysulim = -1.0;  /* calculate it later */
  else 
    mefinfo->sysulim = mefinfo->zp - sysulim;  /* user-supplied */
  
  if(sysllim < 0)
    mefinfo->sysllim = -1.0;  /* calculate it later */
  else 
    mefinfo->sysllim = mefinfo->zp - sysllim;  /* user-supplied */

  ffgkyd(fits, "MJD-OBS", &mjd, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkyd(fits, "MJD", &mjd, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkyd(fits, "JD", &mjd, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
        status = 0;
        mjd = J2K;

        if(verbose) {
          printf("Warning: no MJDs available, timestamps will be garbage\n");
          mefinfo->warned |= WARNED_TIME;
        }
      }
      else if(status) {
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

  /* Correct to midpoint of observation */
  mjd += 0.5 * exptime / 86400.0;

  /* Telescope/instrument-specific kludges */
  ffgkys(fits, "INSTRUME", inst, (char *) NULL, &status);
  ffgkys(fits, "TELESCOP", tel, (char *) NULL, &status);
  if(status == KEY_NO_EXIST)
    status = 0;
  else if(status) {
    fitsio_err(errstr, status, "ffgkys: INSTRUME/TELESCOP");
    goto error;
  }
  else {
    /* ESO WFI needs a scattered light correction too */
    if(!strcasecmp(inst, "WFI") && !strcasecmp(tel, "MPI-2.2"))
      scatcoeff = -1.5;
    
    /* WFCAM has incorrect gain in fits headers */
    if(!strcasecmp(inst, "WFCAM") && !strcasecmp(tel, "UKIRT"))
      gain /= 1.2;
  }

  /* Get block size for row I/O */
  ffgrsz(fits, &rblksz, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgrsz");
    goto error;
  }
  
  /* Allocate column buffers */
  xbuf = (double *) malloc(2 * rblksz * sizeof(double));
  allfluxbuf = (float *) malloc((6+NFLUX) * rblksz * sizeof(float));
  if(!xbuf || !allfluxbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }
  
  ybuf = xbuf + rblksz;

  pkhtbuf = allfluxbuf + NFLUX * rblksz;
  clsbuf = pkhtbuf + rblksz;
  a7buf = pkhtbuf + 2 * rblksz;
  locskybuf = pkhtbuf + 3 * rblksz;
  pmabuf = pkhtbuf + 4*rblksz;
  pmdbuf = pkhtbuf + 5*rblksz;
  
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
    
    ffgcvd(fits, gcols[0], roff + 1L, 1L, rread, 0.0, xbuf, (int *) NULL, &status);
    ffgcvd(fits, gcols[1], roff + 1L, 1L, rread, 0.0, ybuf, (int *) NULL, &status);
    ffgcve(fits, gcols[2], roff + 1L, 1L, rread, 0.0, pkhtbuf, (int *) NULL, &status);
    ffgcve(fits, gcols[3], roff + 1L, 1L, rread, 0.0, clsbuf, (int *) NULL, &status);
    ffgcve(fits, gcols[4], roff + 1L, 1L, rread, 0.0, a7buf, (int *) NULL, &status);
    ffgcve(fits, gcols[5], roff + 1L, 1L, rread, 0.0, locskybuf, (int *) NULL, &status);

    if(gcols[6])
      ffgcve(fits, gcols[6], roff + 1L, 1L, rread, 0.0, pmabuf, (int *) NULL, &status);
    if(gcols[7])
      ffgcve(fits, gcols[7], roff + 1L, 1L, rread, 0.0, pmdbuf, (int *) NULL, &status);

    for(col = 0; col < NFLUX; col++) {
      ffgcve(fits, gcols[8+col], roff + 1L, 1L, rread, 0.0, allfluxbuf + col*rblksz,
	     (int *) NULL, &status);
    }

    if(status) {
      fitsio_err(errstr, status, "ffgcve");
      goto error;
    }
    
    for(rin = 0; rin < rread; rin++) {
      rrin = roff + rin;

      stars[rout].ptr = rrin + 1;
      stars[rout].x = xbuf[rin];
      stars[rout].y = ybuf[rin];

      wcs_xy2vec(&wcs, xbuf[rin], ybuf[rin], v);

      rr = sqrt(v[0]*v[0]+v[1]*v[1]);

      if(gcols[6] && gcols[7] &&
         pmabuf[rin] != 0.0 && pmdbuf[rin] != 0.0) {
        stars[rout].pmra = pmabuf[rin];
        stars[rout].pmdec = pmdbuf[rin];
        stars[rout].havepm = 1;

        if(rr > 0) {
          dvdt[0] = -(pmabuf[rin]*v[1] + pmdbuf[rin]*v[0]*v[2])/(JYR*rr);
          dvdt[1] =  (pmabuf[rin]*v[0] - pmdbuf[rin]*v[1]*v[2])/(JYR*rr);
        }
        else {
          dvdt[0] = 0;
          dvdt[1] = 0;
        }

        dvdt[2] = pmdbuf[rin]*rr/JYR;
      }
      else {
        stars[rout].pmra = 0.0;
        stars[rout].pmdec = 0.0;
        stars[rout].havepm = 0;

        dvdt[0] = 0;
        dvdt[1] = 0;
        dvdt[2] = 0;
      }

      /* UTC epoch is used rather than the correct TDB, not likely to matter. */
      source_star_vec(&(stars[rout].src), v, dvdt, 0, 2000.0 + (mjd-J2K)/JYR);

      stars[rout].ra = atan2(v[1], v[0]);
      if(stars[rout].ra < 0.0)
	stars[rout].ra += TWOPI;

      stars[rout].dec = atan2(v[2], rr);

      stars[rout].cls = lrintf(clsbuf[rin]);
      stars[rout].bflag = (a7buf[rin] < 0.0 ? 1 : 0);
      stars[rout].cflag = 0;

      stars[rout].used = 0;

      stars[rout].apradius = 1.0;  /* default = rcore */

      stars[rout].ref.sky = locskybuf[rin]-pedestal;
      stars[rout].ref.peak = pkhtbuf[rin]+locskybuf[rin]-pedestal;
      
      if(stars[rout].ref.peak > 0.95*satlev)
	stars[rout].ref.satur = 1;
      else
	stars[rout].ref.satur = 0;
      
      wcs_xy2tp(&wcs, xbuf[rin], ybuf[rin], &xi, &xn);

      for(col = 0; col < NFLUX; col++) {
	fluxbuf = allfluxbuf + col*rblksz;
	
	/* Apply scattered light correction */
	if(scatcoeff != 0.0) {
	  fluxbuf[rin] *= powf(10.0, 0.4 * scatcoeff * (xi*xi + xn*xn) * RAD_TO_DEG * RAD_TO_DEG);
	}

	stars[rout].ref.aper[col].flux = fluxbuf[rin] * apcor[col] * percorr;
	stars[rout].ref.aper[col].fluxvar = fabsf(fluxbuf[rin]) * apcor[col] / gain;
	/* sky contribution ? only affects normalisation I think but not sure */

	/* Store reference magnitude for this star */
	if(col == REFAP) {
          if(stars[rout].ref.aper[col].flux > 0.0)
            stars[rout].refmag = 2.5 * log10f(stars[rout].ref.aper[col].flux);
          else
            stars[rout].refmag = -999.0;
        }
      }

      if(pkhtbuf[rin]+locskybuf[rin]-pedestal > 0.95*satlev) {
	sattmp[nsattmp] = stars[rout].ref.aper[REFAP].flux;
	nsattmp++;
      }

      rout++;
    }
    
    roff += rread;
    remain -= rread;
  }

  /* Free workspace */
  free((void *) xbuf);
  xbuf = (double *) NULL;
  free((void *) allfluxbuf);
  allfluxbuf = (float *) NULL;

  /* Determine saturation level robustly - lower quartile */
  if(nsattmp > 0)
    satflux = fquickselect(sattmp, nsattmp/4, nsattmp);

  /* Decide which are usable as comparison stars */
  systematic_select(stars, nrows, mefinfo->sysllim, mefinfo->sysulim);

  /* Consolidate the list if requested */
  if(wantoutcls) {
    rout = 0;

    for(rin = 0; rin < nrows; rin++)
      if(stars[rin].compok ||
	 stars[rin].cls == outcls) {
	/* Keep */
	memmove(&(stars[rout]), &(stars[rin]), sizeof(struct lc_star));
	rout++;
      }

    nstars = rout;

    /* Release any unused star list space.  No free on zero is a kludge
       to avoid creating a null pointer. */
    if(nstars > 0 && nstars < nrows) {
      stars = (struct lc_star *) realloc(stars, nstars * sizeof(struct lc_star));
      if(!stars) {
	report_syserr(errstr, "realloc");
	goto error;
      }
    }
  }
  else
    nstars = nrows;

  mefinfo->stars = stars;
  mefinfo->nstars = nstars;
  mefinfo->nrows = nrows;
  if(satflux > 0.0)
    mefinfo->satmag = mefinfo->zp - 2.5 * log10f(satflux);
  else
    mefinfo->satmag = -999.0;
  mefinfo->reffang = fang;
  mefinfo->havefang = 1;
  mefinfo->refexp = exptime;
  mefinfo->refextinct = noexp ? 0.0 : extinct;
  mefinfo->refairmass = airmass;
  mefinfo->refmagzpt = magzpt;
  mefinfo->refsigma = mefinfo->theosky ? sqrtf(skylev / gain) : skynoise;

  tmp = 5.0 * sqrtf(M_PI * rcore * rcore) * skynoise * apcor[0];
  if(tmp > 0.0)
    mefinfo->refflim = mefinfo->zp - 2.5 * log10f(tmp);
  else
    mefinfo->refflim = -999.0;

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
  if(allfluxbuf)
    free((void *) allfluxbuf);
  if(sattmp)
    free((void *) sattmp);

  return(1);
}

int read_cat (char *catfile, int iframe, int mef, struct lc_mef *mefinfo,
	      struct buffer_info *buf,
	      struct dtai_table *dtab,
	      struct iers_table *itab,
	      struct jpleph_table *jtab,
	      struct jpleph_table *ttab,
	      int dointra, struct intra *icorr,
	      int doinstvers, struct instvers *instverslist, int ninstvers,
	      int diffmode, float satlev,
	      char *errstr) {
  fitsfile *fits;
  int status = 0;

  char *colnames[7] = { "Isophotal_flux", "X_coordinate", "Y_coordinate",
                        "Peak_height", "Areal_1_profile", "Skylev", "Skyrms" };
  char *optcolnames[1] = { "Bad_pixels" };
  int gcols[8+NFLUX], col, collim, optcollim;

  struct lc_point *points = (struct lc_point *) NULL;

  double *xbuf = (double *) NULL, *ybuf;
  float *allfluxbuf = (float *) NULL;
  float *isobuf, *pkhtbuf, *areal1buf, *locskybuf, *skyrmsbuf, *badpixbuf;
  float *fluxbuf;

  struct wcs_info wcs;
  double fang;

  float seeing, ellipt, pedestal, skylev, skynoise, exptime, rcore, gain, percorr;
  float skyvar, area, tmp, expfac;
  double mjd;
  int lmjd, iy, im, id;

  float apcor[NFLUX];

  long nrows, rblksz, roff, remain, routoff, rin, rrin, rout, rrout, rread;
  float flux, fluxvar, locsky, skyfiterr, peak;
  unsigned char measured;

  char inst[FLEN_VALUE], tel[FLEN_VALUE];
  float scatcoeff = 0.0;
  double xi, xn;

  char latstr[FLEN_VALUE] = { '\0' }, lonstr[FLEN_VALUE] = { '\0' };
  char heightstr[FLEN_VALUE] = { '\0' };
  double lat, lon, height = 0;
  unsigned char doairm;

  struct observer obs;
  int utcdn;
  double utcdf, dtt;

  double s[3], pr, seq[3];

  float diam, sc, err, var, avskyfiterr, avscint;

  float *zpoffbuf[NFLUX];
  long nzpoff[NFLUX];
  float avzpoff, avzprms;

  float *skyfiterrbuf = (float *) NULL;
  long navskyfiterr, navscint;

  int cats_are_80 = 0;
  char *ep;

  long split_iexp = 0, split_nexp = -1, rtstat = -1;
  float tamb = -999, humid = -999, press = -999, skytemp = -999;
  int iha = 0;

  char filter[FLEN_VALUE];
  float zpcorr, airmass = 1.0, extinct = 0.0;
  float magzpt;
  long numzpt;
  int haveimgzpt;
  int l1, l2, i, ilim, noexp = 0;

  long schpri;
  float schcad;
  char schtype[FLEN_VALUE];

  struct instvers *instvers = (struct instvers *) NULL;

  long cadencenum = -999;

  /* Init */
  for(col = 0; col < NFLUX; col++) {
    zpoffbuf[col] = (float *) NULL;
    nzpoff[col] = 0;
  }

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

  /* Check for mismatches */
  if(nrows != mefinfo->nrows) {
    report_err(errstr,
	       "number of stars does not match: %d != %d",
	       nrows, mefinfo->nrows);
    goto error;
  }

  /* Simple test for the new 80-column format */
  ffgcno(fits, CASEINSEN, flux_keys_80[0], &cats_are_80, &status);
  if(status == COL_NOT_FOUND || status == COL_NOT_UNIQUE)
    status = 0;
  else if(status) {
    fitsio_err(errstr, status, "ffgcno: %s", flux_keys_80[0]);
    goto error;
  }
  else {
    cats_are_80 = 1;

    colnames[5] = "Sky_level";  /* these were changed in 80-col */
    colnames[6] = "Sky_rms";
    optcolnames[0] = "Error_bit_flag";
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
    ffgcno(fits, CASEINSEN, cats_are_80 ? flux_keys_80[col] : flux_keys_32[col],
	   &(gcols[collim+col]), &status);
    if(status == COL_NOT_UNIQUE)
      status = 0;  /* ignore */
    else if(status) {
      fitsio_err(errstr, status, "ffgcno: %s", 
		 cats_are_80 ? flux_keys_80[col] : flux_keys_32[col]);
      goto error;
    }
  }

  /* Read WCS info */
  if(read_wcs(fits, &wcs, 0, errstr))
    goto error;

  fang = atan2(wcs.b, wcs.a);

  if(satlev < 0) {
    /* Get saturation level - tries SATLEV first.  On MEarth this
     * is a better estimate, based on the non-linearity curve.
     */
    ffgkye(fits, "SATLEV", &satlev, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkye(fits, "SATURATE", &satlev, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
	status = 0;
	satlev = 65535;

	if(verbose > 1 && !diffmode &&
           !(mefinfo->warned & WARNED_SATLEV)) {
	  printf("Warning: using default satlev = %.1f\n", satlev);
          mefinfo->warned |= WARNED_SATLEV;
        }
      }
      else if(status) {
	fitsio_err(errstr, status, "ffgkye: SATURATE");
	goto error;
      }
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: SATLEV");
      goto error;
    }
  }

  /* Read keywords for photometry */
  ffgkye(fits, "EXPTIME", &exptime, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkye(fits, "EXPOSED", &exptime, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkye(fits, "EXP_TIME", &exptime, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
        status = 0;
        exptime = 0;

        if(verbose && !(mefinfo->warned & WARNED_EXPTIME)) {
          printf("Warning: no exposure time found\n");
          mefinfo->warned |= WARNED_EXPTIME;
        }
      }
      else if(status) {
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

  if(exptime && mefinfo->refexp)
    expfac = mefinfo->refexp / exptime;
  else
    expfac = 1.0;

  ffgkye(fits, "SEEING", &seeing, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    seeing = -1.0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: SEEING");
    goto error;
  }

  ffgkye(fits, "ELLIPTIC", &ellipt, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ellipt = -999.0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: ELLIPTIC");
    goto error;
  }

  ffgkye(fits, "PEDESTAL", &pedestal, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    pedestal = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: PEDESTAL");
    goto error;
  }

  ffgkye(fits, "SKYLEVEL", &skylev, (char *) NULL, &status);
  if(status == KEY_NO_EXIST && nrows == 0) {
    status = 0;
    skylev = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: SKYLEVEL");
    goto error;
  }

  ffgkye(fits, "SKYNOISE", &skynoise, (char *) NULL, &status);
  if(status == KEY_NO_EXIST && nrows == 0) {
    status = 0;
    skynoise = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: SKYNOISE");
    goto error;
  }

  ffgkye(fits, "RCORE", &rcore, (char *) NULL, &status);
  if(status == KEY_NO_EXIST && nrows == 0) {
    status = 0;
    rcore = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: RCORE");
    goto error;
  }

  if(verbose > 0 && rcore != mefinfo->refrcore &&
     !(mefinfo->warned & WARNED_RCORE)) {
    printf("Warning: rcore does not match reference: %.3f != %.3f\n",
	   rcore, mefinfo->refrcore);
    mefinfo->warned |= WARNED_RCORE;
  }

  ffgkye(fits, "GAIN", &gain, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkye(fits, "HIERARCH ESO DET OUT1 GAIN", &gain, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkye(fits, "EGAIN", &gain, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
	status = 0;
	gain = 1.0;  /* !!! */

        if(verbose && !(mefinfo->warned & WARNED_GAIN)) {
          printf("Warning: using default gain = %.1f\n", gain);
          mefinfo->warned |= WARNED_GAIN;
        }
      }
      else if(status) {
	fitsio_err(errstr, status, "ffgkye: EGAIN");
	goto error;
      }
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: HIERARCH ESO DET OUT1 GAIN");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: GAIN");
    goto error;
  }

  /* How many standards in image? */
  ffgkyj(fits, "IMNZPT", &numzpt, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    numzpt = 0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkyj: IMNZPT");
    goto error;
  }

  if(numzpt > 0) {
    /* Try to get image ZP solution */
    ffgkye(fits, "IMGZPT", &magzpt, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      haveimgzpt = 0;
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkyj: IMGZPT");
      goto error;
    }
    else if(magzpt == 99.0)  /* flag value */
      haveimgzpt = 0;
    else
      haveimgzpt = 1;
  }
  else
    haveimgzpt = 0;
  
  if(!haveimgzpt) {
    ffgkye(fits, "MAGZPT", &magzpt, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkye(fits, "ZMAG", &magzpt, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
        status = 0;
        magzpt = 25.0;
        
        if(verbose && !(mefinfo->warned & WARNED_MAGZPT)) {
          printf("Warning: using default magzpt = %.1f\n", magzpt);
          mefinfo->warned |= WARNED_MAGZPT;
        }
      }
      else if(status) {
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
  }

  for(col = 0; col < NFLUX; col++) {
    ffgkye(fits, cats_are_80 ? apcor_keys_80[col] : apcor_keys_32[col],
	   &(apcor[col]), (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      apcor[col] = mefinfo->apcor[col];  /* as a backup */
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: %s", 
		 cats_are_80 ? apcor_keys_80[col] : apcor_keys_32[col]);
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
  ffgkys(fits, "FILTER", filter, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkys(fits, "WFFBAND", filter, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkys(fits, "HIERARCH ESO INS FILT1 NAME", filter, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
	status = 0;
	ffgkys(fits, "FILTER2", filter, (char *) NULL, &status);
	if(status == KEY_NO_EXIST) {
	  status = 0;
	  ffgkys(fits, "INSFILTE", filter, (char *) NULL, &status);
          if(status == KEY_NO_EXIST) {
            status = 0;
            filter[0] = '\0';
          }
	  else if(status) {
	    fitsio_err(errstr, status, "ffgkye: INSFILTE");
	    goto error;
	  }
	}
	else if(status) {
	  fitsio_err(errstr, status, "ffgkye: FILTER2");
	  goto error;
	}
      }
      else if(status) {
	fitsio_err(errstr, status, "ffgkye: HIERARCH ESO INS FILT1 NAME");
	goto error;
      }
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkye: WFFBAND");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: FILTER");
    goto error;
  }

  /* Copy */
  strncpy(mefinfo->filter, filter, sizeof(mefinfo->filter)-1);
  mefinfo->filter[sizeof(mefinfo->filter)-1] = '\0';

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

  /* Correction required to account for differences in exposure and extinction */
  if(noexp)
    zpcorr = 0.0;
  else {
    zpcorr = (mefinfo->refairmass - 1.0)*mefinfo->refextinct - (airmass - 1.0)*extinct;
    if(exptime && mefinfo->refexp)
      zpcorr += 2.5 * log10f(exptime/mefinfo->refexp);
  }

  ffgkyd(fits, "MJD-OBS", &mjd, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkyd(fits, "MJD", &mjd, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkyd(fits, "JD", &mjd, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
        status = 0;
        mjd = J2K;

        if(verbose && !(mefinfo->warned & WARNED_TIME)) {
          printf("Warning: no MJDs available, timestamps will be garbage\n");
          mefinfo->warned |= WARNED_TIME;
        }
      }
      else if(status) {
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

  /* Correct to midpoint of observation */
  mjd += 0.5 * exptime / 86400.0;

  /* Get telescope location for airmass calculation */
  doairm = 1;
  
  if(!jtab)  /* can't do without JPL at the moment */
    doairm = 0;

  ffgkys(fits, "LATITUDE", latstr, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkys(fits, "OBSLAT", latstr, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkys(fits, "LAT-OBS", latstr, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
	status = 0;
	ffgkys(fits, "SITELAT", latstr, (char *) NULL, &status);
	if(status == KEY_NO_EXIST) {
          status = 0;
          ffgkys(fits, "HIERARCH ESO TEL GEOLAT", latstr, (char *) NULL,
                 &status);
          if(status == KEY_NO_EXIST) {
            status = 0;
            doairm = 0;
          }
          else if(status) {
            fitsio_err(errstr, status, "ffgkys: HIERARCH ESO TEL GEOLAT");
            goto error;
          }
	}
	else if(status) {
	  fitsio_err(errstr, status, "ffgkys: SITELAT");
	  goto error;
	}
      }
      else if(status) {
	fitsio_err(errstr, status, "ffgkys: LAT-OBS");
	goto error;
      }
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkys: OBSLAT");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkys: LATITUDE");
    goto error;
  }

  if(latstr[0]) {
    lat = strtod(latstr, &ep) * DEG_TO_RAD;
    if(*ep != '\0') {
      if(base60_to_10(latstr, &ep, ":", UNIT_DEG, &lat, UNIT_RAD))
	report_err(errstr, "could not understand latitude: %s", latstr);
    }
  }
  
  ffgkys(fits, "LONGITUD", lonstr, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkys(fits, "OBSLONG", lonstr, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkys(fits, "LONG-OBS", lonstr, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
	status = 0;
	ffgkys(fits, "SITELONG", lonstr, (char *) NULL, &status);
	if(status == KEY_NO_EXIST) {
          status = 0;
          ffgkys(fits, "HIERARCH ESO TEL GEOLON", lonstr, (char *) NULL,
                 &status);
          if(status == KEY_NO_EXIST) {
            status = 0;
            doairm = 0;
          }
          else if(status) {
            fitsio_err(errstr, status, "ffgkys: HIERARCH ESO TEL GEOLON");
            goto error;
          }
	}
	else if(status) {
	  fitsio_err(errstr, status, "ffgkys: SITELONG");
	  goto error;
	}
      }
      else if(status) {
	fitsio_err(errstr, status, "ffgkys: LONG-OBS");
	goto error;
      }
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkys: OBSLONG");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkys: LONGITUD");
    goto error;
  }

  if(lonstr[0]) {
    lon = strtod(lonstr, &ep) * DEG_TO_RAD;
    if(*ep != '\0') {
      if(base60_to_10(lonstr, &ep, ":", UNIT_DEG, &lon, UNIT_RAD))
	report_err(errstr, "could not understand longitude: %s", latstr);
    }
  }
  
  ffgkys(fits, "HEIGHT", heightstr, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkys(fits, "OBSALT", heightstr, (char *) NULL, &status);
    if(status == KEY_NO_EXIST) {
      status = 0;
      ffgkys(fits, "ALT-OBS", heightstr, (char *) NULL, &status);
      if(status == KEY_NO_EXIST) {
	status = 0;
	ffgkys(fits, "SITEALT", heightstr, (char *) NULL, &status);
	if(status == KEY_NO_EXIST) {
          status = 0;
          ffgkys(fits, "HIERARCH ESO TEL GEOELEV", heightstr, (char *) NULL,
                 &status);
          if(status == KEY_NO_EXIST) {
            status = 0;
            height = 0;
          }
          else if(status) {
            fitsio_err(errstr, status, "ffgkys: HIERARCH ESO TEL GEOELEV");
            goto error;
          }
	}
	else if(status) {
	  fitsio_err(errstr, status, "ffgkys: SITEALT");
	  goto error;
	}
      }
      else if(status) {
	fitsio_err(errstr, status, "ffgkys: ALT-OBS");
	goto error;
      }
    }
    else if(status) {
      fitsio_err(errstr, status, "ffgkys: OBSALT");
      goto error;
    }
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkys: HEIGHT");
    goto error;
  }

  if(heightstr[0]) {
    height = strtod(heightstr, &ep);
    if(*ep != '\0' && *ep != 'm') {
      report_err(errstr, "could not understand height: %s", heightstr);
      goto error;
    }
  }

  /* Check for telescope-specific things */
  ffgkys(fits, "INSTRUME", inst, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    ffgkys(fits, "DETECTOR", inst, (char *) NULL, &status);
    if(status == KEY_NO_EXIST)
      status = 0;
    else if(status) {
      fitsio_err(errstr, status, "ffgkys: INSTRUME/DETECTOR");
      goto error;
    }
  }
  ffgkys(fits, "TELESCOP", tel, (char *) NULL, &status);
  if(status == KEY_NO_EXIST)
    status = 0;
  else if(status) {
    fitsio_err(errstr, status, "ffgkys: TELESCOP");
    goto error;
  }
  else {
    if(!strcasecmp(inst, "WFC") && !strcasecmp(tel, "INT"))
      diam = 2540;
    else if(!strcasecmp(inst, "WFI") && !strcasecmp(tel, "MPI-2.2")) {
      scatcoeff = -1.5;  /* ESO WFI needs a scattered light correction too */
      diam = 2200;
    }
    else if(!strcasecmp(inst, "WFCAM") && !strcasecmp(tel, "UKIRT")) {
      gain /= 1.2;  /* WFCAM has incorrect gain in fits headers */
      diam = 3803;
    }
    else if(!strcasecmp(inst, "CCDMosaThin1"))
      diam = 3934;  /* KPNO Mosaic */
    else if(!strcasecmp(inst, "Mosaic2"))
      diam = 3934;  /* CTIO Mosaic II */
    else if(!strcasecmp(inst, "Apogee Alta"))
      diam = 400;  /* default for MEarth - it has the actual value in APTDIA */
  }

  /* Try to read telescope aperture diameter */
  ffgkye(fits, "APTDIA", &diam, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    diam = -1.0;  /* flag as dunno */
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: APTDIA");
    goto error;
  }

  /* MEarth-specific: exposure grouping */
  ffgkyj(fits, "IEXP", &split_iexp, (char *) NULL, &status);
  ffgkyj(fits, "NEXP", &split_nexp, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    split_iexp = 0;
    split_nexp = -1;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkyj: IEXP/NEXP");
    goto error;
  }

  /* MEarth-specific: weather parameters */
  ffgkye(fits, "TEMPERAT", &tamb, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    tamb = -999;
  }

  ffgkye(fits, "HUMIDITY", &humid, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    humid = -999;
  }

  ffgkye(fits, "PRESSURE", &press, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    press = -999;
  }

  ffgkye(fits, "SKYTEMP", &skytemp, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    skytemp = -999;
  }

  if(status) {
    fitsio_err(errstr, status, "ffgkye: TEMPERAT/HUMIDITY/PRESSURE/SKYTEMP");
    goto error;
  }

  /* MEarth-specific: realtime status */
  ffgkyj(fits, "RTSTAT", &rtstat, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    rtstat = -1;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkyj: RTSTAT");
    goto error;
  }

  /* MEarth-specific: scheduling information */
  ffgkyj(fits, "SCHPRI", &schpri, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    schpri = -1;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkyj: SCHPRI");
    goto error;
  }

  ffgkye(fits, "SCHCAD", &schcad, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    schcad = -999.0;
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkye: SCHCAD");
    goto error;
  }

  ffgkys(fits, "SCHTYPE", schtype, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    schtype[0] = '\0';
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkys: SCHTYPE");
    goto error;
  }

  /* MEarth-specific: figure out "instrument version" */
  if(doinstvers) {
    /* First, we need the "night of" date.  This is essentially the
     * same algorithm as the observing system uses (we don't use the
     * filename because it's not guaranteed to have been maintained)
     * but it's not a perfect replica because we use MJD here and
     * the observing system uses zone time.
     */
    lmjd = floor(mjd + lon*RAD_TO_HR/24 - 0.5);
    mjd2date(lmjd, &iy, &im, &id);
    id += iy*10000 + im*100;

    /* Now look for it in the table */
    for(i = 0; i < ninstvers; i++)
      if(id < instverslist[i].date)
	break;

    if(i > 0)
      instvers = instverslist + i-1;
  }

  /* Kepler-specific: cadence number */
  ffgkyj(fits, "CADENCE", &cadencenum, (char *) NULL, &status);
  if(status == KEY_NO_EXIST) {
    status = 0;
    cadencenum = -999;  /* valid numbers should be positive */
  }
  else if(status) {
    fitsio_err(errstr, status, "ffgkyj: CADENCE");
    goto error;
  }

  if(doairm) {
    /* Create observer structure */
    observer_init(&obs, lon, lat, height);

    refract_const(tamb != -999 ? 273.16+tamb : 283.0,
		  humid != -999 ? humid/100 : 0.5,
		  press != -999 ? press : 1013.25*exp(-height/(29.3*273.0)),
		  0.80,  /* wavelength, microns */
		  height,
		  obs.refco);
  }
  else
    observer_geoc(&obs);  /* geocentric corrections only */

  if(dtab) {
    /* Look up delta(TT) */
    utcdn = floor(mjd);
    utcdf = mjd - utcdn;
    
    dtt = dtai(dtab, utcdn, utcdf) + DTT;
  }
  else
    dtt = DTT;  /* assumes TAI=UTC */
  
  if(jtab) {
    /* Compute time-dependent quantities */
    observer_update(&obs, jtab, ttab, itab, mjd, dtt, OBSERVER_UPDATE_ALL);
  }

  /* Get block size for row I/O */
  ffgrsz(fits, &rblksz, &status);
  if(status) {
    fitsio_err(errstr, status, "ffgrsz");
    goto error;
  }
  
  /* Allocate column buffers */
  xbuf = (double *) malloc(2 * rblksz * sizeof(double));
  allfluxbuf = (float *) malloc((6+NFLUX) * rblksz * sizeof(float));
  if(!xbuf || !allfluxbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }
  
  ybuf = xbuf + rblksz;

  isobuf = allfluxbuf + NFLUX * rblksz;
  pkhtbuf = isobuf + rblksz;
  areal1buf = isobuf + 2 * rblksz;
  locskybuf = isobuf + 3 * rblksz;
  skyrmsbuf = isobuf + 4 * rblksz;
  badpixbuf = isobuf + 5 * rblksz;

  /* Allocate memory for lightcurve points */
  points = (struct lc_point *) malloc(rblksz * sizeof(struct lc_point));
  if(!points) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Allocate memory for medians */
  for(col = 0; col < NFLUX; col++) {
    zpoffbuf[col] = (float *) malloc(nrows * sizeof(float));
    if(!zpoffbuf[col]) {
      report_syserr(errstr, "malloc");
      goto error;
    }
  }

  skyfiterrbuf = (float *) malloc(nrows * sizeof(float));
  if(!skyfiterrbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Read catalogue */
  roff = 0L;
  remain = nrows;
  routoff = 0L;

  navskyfiterr = 0;

  avscint = 0;
  navscint = 0;

  while(remain > 0) {
    rread = (remain > rblksz ? rblksz : remain);
    
    ffgcve(fits, gcols[0], roff + 1L, 1L, rread, 0.0, isobuf, (int *) NULL, &status);
    ffgcvd(fits, gcols[1], roff + 1L, 1L, rread, 0.0, xbuf, (int *) NULL, &status);
    ffgcvd(fits, gcols[2], roff + 1L, 1L, rread, 0.0, ybuf, (int *) NULL, &status);
    ffgcve(fits, gcols[3], roff + 1L, 1L, rread, 0.0, pkhtbuf, (int *) NULL, &status);
    ffgcve(fits, gcols[4], roff + 1L, 1L, rread, 0.0, areal1buf, (int *) NULL, &status);
    ffgcve(fits, gcols[5], roff + 1L, 1L, rread, 0.0, locskybuf, (int *) NULL, &status);
    ffgcve(fits, gcols[6], roff + 1L, 1L, rread, 0.0, skyrmsbuf, (int *) NULL, &status);

    if(gcols[7])
      ffgcve(fits, gcols[7], roff + 1L, 1L, rread, 0.0, badpixbuf, (int *) NULL, &status);

    for(col = 0; col < NFLUX; col++) {
      ffgcve(fits, gcols[8+col], roff + 1L, 1L, rread, 0.0, allfluxbuf + col*rblksz,
	     (int *) NULL, &status);
    }

    rout = 0;

    for(rin = 0; rin < rread; rin++) {
      rrin = roff + rin;
      rrout = routoff + rout;

      /* Do we want this row?  The list is guaranteed to be sorted,
	 so a simple comparison is adequate. */
      if(rrout >= mefinfo->nstars ||
	 mefinfo->stars[rrout].ptr != rrin+1)
	continue;

      /* Figure out if the source was measured.  Unfortunately, zero
         is used as the flag value and is a valid result, so we test
         a bunch of quantities to make it more robust.  These are all
         computed in update, which checks first if the rcore radius
         aperture is entirely on the frame.  The photometry (phopt)
         does no checks so sources with the aperture partly in the
         frame can have measurements in Core_flux and the others.
         The decision was made to discard such cases, a partial
         measurement is probably not useful and we don't have the
         other parameters such as x, y, peak height to go with it
         anyway.  This is a change from old versions of the program
         which just checked for a flux and assumed the object was
         on the frame if there was one.
         NaN values are also trapped.  This can happen if there were
         NaNs in the input image.  In such cases, assume all the
         measurements are bad and proceed as if not measured. */
      if(isfinite(isobuf[rin]) &&
         isfinite(pkhtbuf[rin]) &&
         isfinite(areal1buf[rin]) &&
         (isobuf[rin] || pkhtbuf[rin] || areal1buf[rin]))  /* one non-zero */
        measured = 1;
      else
        measured = 0;

      /* Store x, y positions */
      points[rout].x = xbuf[rin];
      points[rout].y = ybuf[rin];

      if(jtab) {
        /* Apply space motion */
        source_place(&obs, &(mefinfo->stars[rrout].src),
                     jtab, TR_MOTION, s, NULL, &pr);
        
        /* Modified BJD(TDB) */
        points[rout].bjd = obs.tt + (obs.dtdb + bary_delay(&obs, s, pr)) / DAY;
      }
      else
        points[rout].bjd = -999.0;  /* can't compute */

      /* Compute airmass */
      if(doairm) {
	/* Observed place - parallax left out, we don't have it anyway */
	observer_ast2obs(&obs, s, NULL, pr, TR_TO_OBS_AZ);

	points[rout].airmass = v_airmass(s);

	mt_x_v(obs.phm, s, seq);
	points[rout].ha = -atan2(seq[1], seq[0]);
	
	if(points[rout].ha > 0)
	  iha = 1;  /* "average" HA in some sense - will NOT work for wide-field */
      }
      else {
	points[rout].airmass = -999.0;  /* flag unusability */
	points[rout].ha = -999.0;
      }

      if(diffmode) {
	locsky = mefinfo->stars[rrout].ref.sky + locskybuf[rin] - pedestal;
	peak = mefinfo->stars[rrout].ref.peak + pkhtbuf[rin]+locskybuf[rin] - pedestal; /* ?? */

	points[rout].satur = mefinfo->stars[rrout].ref.satur;
      }
      else {
	locsky = locskybuf[rin] - pedestal;
	peak = pkhtbuf[rin]+locskybuf[rin] - pedestal;
	
	if(peak > 0.95*satlev || mefinfo->stars[rrout].ref.satur)
	  points[rout].satur = 1;
	else
	  points[rout].satur = 0;
      }

      points[rout].sky = locsky;
      points[rout].peak = peak;

      /* Accumulate sky fitting error */
      skyfiterr = 0.0;  /* no good way to get this yet */

      skyfiterrbuf[navskyfiterr] = skyfiterr;
      navskyfiterr++;

      /* Accumulate counts of frames with bad pixels */
      if(gcols[7] && badpixbuf[rin] > 0.0) {
	points[rout].conf = 1;
	mefinfo->stars[rrout].cflag++;
      }
      else
	points[rout].conf = 0;

      wcs_xy2tp(&wcs, xbuf[rin], ybuf[rin], &xi, &xn);

      for(col = 0; col < NFLUX; col++) {
	fluxbuf = allfluxbuf + col*rblksz;

	/* Apply scattered light correction */
	if(scatcoeff != 0.0) {
	  fluxbuf[rin] *= powf(10.0, 0.4 * scatcoeff * (xi*xi + xn*xn) * RAD_TO_DEG * RAD_TO_DEG);
	}

	area = M_PI * rcore * rcore * flux_apers[col] * flux_apers[col];

        if(mefinfo->theosky)
          skyvar = skylev * area / gain +
                   skyfiterr * skyfiterr * area * area;
        else
          skyvar = skynoise * skynoise * area +
                   skyfiterr * skyfiterr * area * area;

	if(diffmode) {
	  flux = fluxbuf[rin] * mefinfo->apcor[col] * mefinfo->percorr;
	  fluxvar = fabsf(fluxbuf[rin]) * mefinfo->apcor[col] / gain + skyvar;

	  if(mefinfo->stars[rrout].ref.aper[col].flux == 0.0)
	    flux = 0.0;
	  else
	    flux += mefinfo->stars[rrout].ref.aper[col].flux;

	  fluxvar += mefinfo->stars[rrout].ref.aper[col].fluxvar;
	}
	else {
	  flux = fluxbuf[rin] * apcor[col] * percorr;
	  fluxvar = fabsf(fluxbuf[rin]) * apcor[col] / gain + skyvar;
	}

	if(measured && flux > 0.0) {
	  points[rout].aper[col].flux = 2.5 * log10f(flux);

	  if(!diffmode)
	    points[rout].aper[col].flux -= zpcorr;

	  /* Apply intrapixel correction if requested */
	  if(dointra)
	    points[rout].aper[col].flux += calc_intra(xbuf[rin], ybuf[rin], icorr);

          /* If suitable, accumulate for median */
          if(mefinfo->stars[rrout].compok &&
             !points[rout].satur) {
            *(zpoffbuf[col] + nzpoff[col])
              = points[rout].aper[col].flux - mefinfo->stars[rrout].refmag;
            nzpoff[col]++;
          }

          /* Relative variance */
          var = fluxvar / (flux*flux);

	  /* Scintillation */
	  if(diam > 0.0 && points[rout].airmass > 0.0 && exptime) {
	    sc = 0.09 * powf(diam / 10.0, -2.0 / 3.0) *
 	                powf(points[rout].airmass, 3.0 / 2.0) *
 	                expf(-height / 8000.0) /
	                sqrtf(2*exptime);

	    avscint += 2.5 * M_LOG10E * sc;  /* flux to mag */
	    navscint++;

	    var += sc*sc;
	  }

          /* Convert to magnitudes.  Bright side error bar, see
             Naylor et al. 2002, MNRAS, 335, 291 section 10. */
          err = 2.5*log10f(1.0 + sqrtf(var));

	  /* Convert to "virtual variance" and stash result */
	  points[rout].aper[col].fluxvar = err*err;
	  points[rout].aper[col].fluxvarcom = points[rout].aper[col].fluxvar;
	}
	else {
	  points[rout].aper[col].flux = 0.0;
	  points[rout].aper[col].fluxvar = 0.0;
	  points[rout].aper[col].fluxvarcom = 0.0;
	}

	points[rout].aper[col].wt = 0;
      }

      rout++;
    }

    if(rout > 0) {
      /* Write out those */
      if(buffer_put_frame(buf, points, routoff, rout, iframe, errstr))
	goto error;
    }

    if(status) {
      fitsio_err(errstr, status, "ffgcve");
      goto error;
    }
    
    roff += rread;
    remain -= rread;
    routoff += rout;
  }

  /* Free workspace */
  free((void *) xbuf);
  xbuf = (double *) NULL;
  free((void *) allfluxbuf);
  allfluxbuf = (float *) NULL;

  free((void *) points);
  points = (struct lc_point *) NULL;

  /* Check we got everything */
  if(routoff != mefinfo->nstars) {
    report_err(errstr,
	       "number of stars does not match: %d != %d",
	       routoff, mefinfo->nstars);
    goto error;
  }

  /* Accumulate averages of noise contributions */
  if(mefinfo->theosky)
    tmp = skylev * expfac / gain;
  else
    tmp = skynoise * skynoise * expfac;
  
  if(diffmode)
    mefinfo->avsigma += tmp + mefinfo->refsigma * mefinfo->refsigma;
  else
    mefinfo->avsigma += tmp;

  if(navskyfiterr > 0) {
    fmedsig(skyfiterrbuf, navskyfiterr, &avskyfiterr, (float *) NULL);
    mefinfo->avskyfit += avskyfiterr * avskyfiterr * expfac;
  }

  mefinfo->avapcor += apcor[0];

  if(navscint > 0)
    mefinfo->avscint += avscint / navscint;

  /* Store this frame MJD and seeing */
  mefinfo->frames[iframe].mjd = mjd;
  mefinfo->frames[iframe].exptime = exptime;
  mefinfo->frames[iframe].seeing = seeing;
  mefinfo->frames[iframe].ellipt = ellipt;
  mefinfo->frames[iframe].skylev = skylev-pedestal;
  mefinfo->frames[iframe].skynoise = skynoise;
  mefinfo->frames[iframe].fang = fang;

  if(mefinfo->havefang)
    mefinfo->frames[iframe].iang = lrintf(ranorm(fang - mefinfo->reffang)/M_PI) % 2;
  else
    mefinfo->frames[iframe].iang = iha;  /* old, incorrect, method */

  mefinfo->frames[iframe].split_iexp = split_iexp;
  mefinfo->frames[iframe].split_nexp = split_nexp;

  mefinfo->frames[iframe].tamb = tamb;
  mefinfo->frames[iframe].humid = humid;
  mefinfo->frames[iframe].press = press;
  mefinfo->frames[iframe].skytemp = skytemp;

  mefinfo->frames[iframe].rtstat = rtstat;
  mefinfo->frames[iframe].schpri = schpri;
  mefinfo->frames[iframe].schcad = schcad;
  memcpy(mefinfo->frames[iframe].schtype, schtype, sizeof(mefinfo->frames[iframe].schtype));
  mefinfo->frames[iframe].cadencenum = cadencenum;

  mefinfo->frames[iframe].zpdiff = magzpt - mefinfo->refmagzpt;

  /* Initialise these (extinc is cumulative) */
  mefinfo->frames[iframe].offset = 0;
  mefinfo->frames[iframe].rms = -999.0;
  mefinfo->frames[iframe].extinc = 0;
  mefinfo->frames[iframe].sigm = 0;

  mefinfo->frames[iframe].instvers = instvers;

  /* Compute first estimate of ZP offset and init systematic fit
     with an educated guess. */
  for(col = 0; col < NFLUX; col++) {
    if(nzpoff[col] > 0) {
      fmedsig(zpoffbuf[col], nzpoff[col], &avzpoff, &avzprms);

      mefinfo->frames[iframe].sys[col].xbar = 0;
      mefinfo->frames[iframe].sys[col].ybar = 0;
      memset(mefinfo->frames[iframe].sys[col].coeff,
             0, sizeof(mefinfo->frames[iframe].sys[col].coeff));
      mefinfo->frames[iframe].sys[col].coeff[0] = avzpoff;

      memset(mefinfo->frames[iframe].sys[col].cov,
             0, sizeof(mefinfo->frames[iframe].sys[col].cov));
      mefinfo->frames[iframe].sys[col].degree = 0;

      mefinfo->frames[iframe].sys[col].medoff = 0;
      mefinfo->frames[iframe].sys[col].sigoff = avzprms;
      mefinfo->frames[iframe].sys[col].sigm = 0;
      mefinfo->frames[iframe].sys[col].npt = nzpoff[col];
    }

    free((void *) zpoffbuf[col]);
    zpoffbuf[col] = (float *) NULL;
  }

  free((void *) skyfiterrbuf);
  skyfiterrbuf = (float *) NULL;

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
  if(allfluxbuf)
    free((void *) allfluxbuf);
  if(points)
    free((void *) points);

  for(col = 0; col < NFLUX; col++)
    if(zpoffbuf[col])
      free((void *) zpoffbuf[col]);

  if(skyfiterrbuf)
    free((void *) skyfiterrbuf);

  return(1);
}
