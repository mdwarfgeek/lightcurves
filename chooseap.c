#include <stdio.h>
#include <stdlib.h>

#include "lightcurves.h"

#include "floatmath.h"
#include "util.h"

/* Parameters: k-sigma for clipping and max number of iterations */
#define SIGCLIP   5
#define NITERMAX  5

int chooseap (struct buffer_info *buf, struct lc_mef *mefinfo,
	      struct lc_point *ptbuf, float *medbuf, char *errstr) {
  int aper, iter;
  long star, opt, m, nmed;

  float *tmpbuf = (float *) NULL;

  float medflux, sigflux, cliplow, cliphigh;
  float loge, area, skycont, readcont, photons, theor, off;
  float medoff, sigoff;

  float medofflist[NFLUX], fluxbound[NFLUX];
  float flux, area1, area2;

  int useaper;

  /* Precalculate parameters */
  loge = log10f(expf(1.0));
  readcont = mefinfo->refreadnois * mefinfo->refreadnois;

  /* Allocate buffers */
  tmpbuf = (float *) malloc(mefinfo->nf * sizeof(float));
  if(!tmpbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  for(aper = 0; aper < NFLUX; aper++) {
    /* Calculate quadrature contribution to RMS from aperture loss */
    area = M_PI * mefinfo->refrcore * mefinfo->refrcore *
                  flux_apers[aper] * flux_apers[aper];

    skycont = area * (mefinfo->refgain * mefinfo->refgain *
		      mefinfo->avsigma * mefinfo->avsigma);

    cliphigh = mefinfo->satclip[aper];
    cliplow = cliphigh - 2.0;

    opt = 0;
    for(star = 0; star < mefinfo->nstars; star++) {
      medflux = mefinfo->stars[star].medflux[aper];
      sigflux = mefinfo->stars[star].sigflux[aper];

      if(medflux > cliplow && medflux < cliphigh &&
	 mefinfo->stars[star].cls == -1) {
	photons = powf(10.0, 0.4 * medflux) * mefinfo->refgain;
	theor = 2.5 * loge * sqrtf(photons + skycont + readcont) / photons;

	off = sigflux*sigflux - theor*theor;

	medbuf[opt] = off;
	opt++;
      }
    }

    nmed = opt;
    memcpy(tmpbuf, medbuf, nmed * sizeof(float));
    medsig(tmpbuf, nmed, &medoff, &sigoff);

    for(iter = 0; iter < NITERMAX; iter++) {
      opt = 0;
      for(m = 0; m < nmed; m++)
	if(fabsf(medbuf[m] - medoff) < SIGCLIP * sigoff) {
	  medbuf[opt] = medbuf[m];
	  opt++;
	}

      nmed = opt;
      memcpy(tmpbuf, medbuf, nmed * sizeof(float));
      medsig(tmpbuf, nmed, &medoff, &sigoff);
    }

    if(verbose)
      printf("Aperture %d offset %.4f +/- %.4f\n", aper+1, sqrtf(medoff), sqrtf(sigoff));

    medofflist[aper] = medoff;
  }

  /* Calculate cross-over points.  Assumes the apertures are ordered by
   * increasing size.
   */
  fluxbound[0] = 0.0;

  for(aper = 1; aper < NFLUX; aper++) {
    if(medofflist[aper-1] > medofflist[aper]) {
      area1 = M_PI * mefinfo->refrcore * mefinfo->refrcore *
	             flux_apers[aper-1] * flux_apers[aper-1];
      area2 = M_PI * mefinfo->refrcore * mefinfo->refrcore *
                     flux_apers[aper] * flux_apers[aper];

      flux = 2.5 * loge * mefinfo->avsigma * sqrtf((area2 - area1) /
						   (medofflist[aper-1] -
						    medofflist[aper]));
      fluxbound[aper] = 2.5 * log10f(flux);
    }
    else {
      /* Don't bother, then */
      fluxbound[aper] = 0.0;
    }

    if(verbose)
      printf("%d -> %d: %.2f\n", aper-1, aper, mefinfo->zp - fluxbound[aper]);
  }

  free((void *) tmpbuf);
  tmpbuf = (float *) NULL;

  /* Choose the relevant aperture for each object */
  for(star = 0; star < mefinfo->nstars; star++) {
    medflux = mefinfo->stars[star].medflux[0];

    useaper = 0;  /* default to the smallest */
    if(medflux > 0.0) {
      for(aper = 0; aper < NFLUX; aper++) {
	if(medflux > fluxbound[aper])
	  useaper = aper;
	else
	  break;
      }
    }

    /* Reshuffle */
    if(useaper > 0) {
      /* Read in measurements for this star in 'useaper' */
      if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, useaper, errstr))
	goto error;

      /* Write out into aperture 0 */
      if(buffer_put_object(buf, ptbuf, 0, mefinfo->nf, star, 0, errstr))
	goto error;
    }
  }

  return(0);

 error:
  if(tmpbuf)
    free((void *) tmpbuf);

  return(1);
}

