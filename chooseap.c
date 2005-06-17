#include <stdio.h>
#include <stdlib.h>

#include "lightcurves.h"

#include "floatmath.h"
#include "util.h"

int chooseap (struct buffer_info *buf, struct lc_mef *mefinfo,
	      struct lc_point *ptbuf, float *medbuf, char *errstr) {
  long pt, star;
  struct lc_point *refbuf = (struct lc_point *) NULL;
  float *corbuf = (float *) NULL, medcor;
  long ncor;

  float avapcor[NFLUX];
  long nav[NFLUX];

  int aper, useaper;
  float rms, rmsmin;

  /* Allocate buffers */
  refbuf = (struct lc_point *) malloc(mefinfo->nstars * sizeof(struct lc_point));
  corbuf = (float *) malloc(mefinfo->nstars * sizeof(float));
  if(!refbuf || !corbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Initialise means */
  for(aper = 0; aper < NFLUX; aper++) {
    avapcor[aper] = 0.0;
    nav[aper] = 0;
  }

  if(verbose)
    printf(" Computing aperture corrections\n");

  /* Compute empirical aperture corrections */
  for(pt = 0; pt < mefinfo->nf; pt++) {
    /* Read in reference aperture */
    if(buffer_fetch_frame(buf, refbuf, 0, mefinfo->nstars, pt, 0, errstr))
      goto error;

    /* Loop through comparison apertures */
    for(aper = 1; aper < NFLUX; aper++) {
      if(buffer_fetch_frame(buf, ptbuf, 0, mefinfo->nstars, pt, aper, errstr))
	goto error;

      /* Compute median aperture correction */
      ncor = 0;
      for(star = 0; star < mefinfo->nstars; star++)
	if(refbuf[star].flux > 0.0 && ptbuf[star].flux > 0.0 &&
	   refbuf[star].fluxerr > 0.0 && ptbuf[star].fluxerr > 0.0 &&
	   !refbuf[star].satur && !ptbuf[star].satur &&
	   mefinfo->stars[star].sigflux[aper] > 0 &&
	   !mefinfo->stars[star].bflag &&
	   mefinfo->stars[star].cls == -1) {
	  corbuf[ncor] = ptbuf[star].flux - refbuf[star].flux;
	  ncor++;
	}

      if(ncor > 1) {
	medsig(corbuf, ncor, &medcor, (float *) NULL);

	/* Accumulate mean */
	avapcor[aper] += medcor;
	nav[aper]++;
      }
    }
  }

  free((void *) refbuf);
  refbuf = (struct lc_point *) NULL;
  free((void *) corbuf);
  corbuf = (float *) NULL;

  /* Compute mean aperture corrections */
  for(aper = 1; aper < NFLUX; aper++) {
    avapcor[aper] /= nav[aper];

    if(verbose)
      printf("  Aperture correction %d->1 = %.4f using %ld frames\n",
	     aper+1, avapcor[aper], nav[aper]);
  }

  /* Loop through all stars */
  for(star = 0; star < mefinfo->nstars; star++) {
    /* Choose the aperture with the lowest RMS scatter */
    useaper = -1;
    rmsmin = 0.0;
    for(aper = 0; aper < NFLUX; aper++) {
      rms = mefinfo->stars[star].sigflux[aper];

      if(useaper < 0 || rms < rmsmin) {
	useaper = aper;
	rmsmin = rms;
      }
    }

    /* Sanity check */
    if(useaper < 0) {
      report_err(errstr, "could not find an aperture for star %ld", star+1);
      goto error;
    }

    /* Reshuffle */
    if(useaper > 0) {
      /* Read in measurements for this star in 'useaper' */
      if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, useaper, errstr))
	goto error;

      /* Apply aperture correction */
      for(pt = 0; pt < mefinfo->nf; pt++)
	if(ptbuf[pt].flux > 0.0 && ptbuf[pt].fluxerr > 0.0)
	  ptbuf[pt].flux -= avapcor[useaper];

      /* Update aperture size */
      mefinfo->stars[star].apradius *= flux_apers[useaper];

      /* Write out into aperture 0 */
      if(buffer_put_object(buf, ptbuf, 0, mefinfo->nf, star, 0, errstr))
	goto error;
    }
  }

  return(0);

 error:
  if(refbuf)
    free((void *) refbuf);
  if(corbuf)
    free((void *) corbuf);

  return(1);
}

