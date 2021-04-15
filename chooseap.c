#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lightcurves.h"

#include "util.h"

int chooseap (struct buffer_info *buf, struct lc_mef *mefinfo,
	      struct lc_point *ptbuf, float *medbuf, int doapcor, char *errstr) {
  long pt, star;
  float *corbuf = (float *) NULL, medcor;
  long ncor;

  float avapcor[NFLUX];
  long nav[NFLUX];

  int aper, useaper;
  float rms, rmsmin;

  /* Initialise means */
  for(aper = 0; aper < NFLUX; aper++) {
    avapcor[aper] = 0.0;
    nav[aper] = 0;
  }

  if(mefinfo->aperture == 0 && doapcor) {
    /* Allocate buffers */
    corbuf = (float *) malloc(mefinfo->nstars * sizeof(float));
    if(!corbuf) {
      report_syserr(errstr, "malloc");
      goto error;
    }
    
    if(verbose)
      printf(" Computing aperture corrections\n");
    
    /* Compute empirical aperture corrections */
    for(pt = 0; pt < mefinfo->nf; pt++) {
      /* Read */
      if(buffer_fetch_frame(buf, ptbuf, 0, mefinfo->nstars, pt, errstr))
	goto error;
      
      /* Loop through comparison apertures */
      for(aper = 0; aper < NFLUX; aper++) {
	if(aper == REFAP)
	  continue;

	/* Compute median aperture correction */
	ncor = 0;
	for(star = 0; star < mefinfo->nstars; star++)
	  if(ptbuf[star].aper[REFAP].flux > 0.0 &&
	     ptbuf[star].aper[aper].flux > 0.0 &&
	     ptbuf[star].aper[REFAP].fluxvar > 0.0 &&
	     ptbuf[star].aper[aper].fluxvar > 0.0 &&
	     !ptbuf[star].satur &&
	     mefinfo->stars[star].sigflux[aper] >= 0 &&
	     mefinfo->stars[star].medflux[aper] >= mefinfo->sysllim &&
	     mefinfo->stars[star].medflux[aper] <= mefinfo->sysulim &&
	     mefinfo->stars[star].cls == -1) {
	    corbuf[ncor] = ptbuf[star].aper[aper].flux - ptbuf[star].aper[REFAP].flux;
	    ncor++;
	  }
	
	if(ncor > 1) {
	  fmedsig(corbuf, ncor, &medcor, (float *) NULL);
	  
	  /* Accumulate mean */
	  if(fabsf(medcor) < 1) {  /* it can't ever be that big => junk */
	    avapcor[aper] += medcor;
	    nav[aper]++;
	  }
	}
      }
    }
    
    free((void *) corbuf);
    corbuf = (float *) NULL;
    
    /* Compute mean aperture corrections */
    for(aper = 0; aper < NFLUX; aper++) {
      if(aper == REFAP)
	continue;

      if(nav[aper] >= 1)
	avapcor[aper] /= nav[aper];
      else
	avapcor[aper] = 0.0;
      
      if(verbose)
	printf("  Aperture correction %d->%d = %.4f using %ld frames\n",
	       aper+1, REFAP+1, avapcor[aper], nav[aper]);
    }
  }    

  /* Loop through all stars */
  for(star = 0; star < mefinfo->nstars; star++) {
    if(mefinfo->aperture == 0) {
      /* Choose the aperture with the lowest RMS scatter */
      useaper = -1;
      rmsmin = 0.0;
      for(aper = 0; aper < NFLUX; aper++) {
	rms = mefinfo->stars[star].sigflux[aper];
	
        /* Uses less than or equal to comparison so we keep increasing the
           aperture size if the rms are all equal.  Larger aperture should
           usually yield lower systematics in such cases.  Also causes the
           1 data point case (zero rms) to use the largest aperture, which
           is desirable for MEarth RT analysis. */
	if(useaper < 0 || rms <= rmsmin) {
	  useaper = aper;
	  rmsmin = rms;
	}
      }
      
      /* Sanity check */
      if(useaper < 0) {
	report_err(errstr, "could not find an aperture for star %ld", star+1);
	goto error;
      }
    }
    else
      useaper = mefinfo->aperture-1;

    mefinfo->stars[star].iap = useaper;
    mefinfo->stars[star].apradius *= flux_apers[useaper];

    /* Apply aperture corrections */
    if(doapcor) {
      /* Read in measurements for this star */
      if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, errstr))
	goto error;
      
      for(aper = 0; aper < NFLUX; aper++) {
	if(aper == REFAP)
	  continue;
	
	/* Apply aperture correction */
	for(pt = 0; pt < mefinfo->nf; pt++)
	  if(ptbuf[pt].aper[aper].flux > 0.0 && ptbuf[pt].aper[aper].fluxvar > 0.0)
	    ptbuf[pt].aper[aper].flux -= avapcor[aper];
	
	/* Correct median flux */
	mefinfo->stars[star].medflux[aper] -= avapcor[aper];
      }
      
      /* Write out */
      if(buffer_put_object(buf, ptbuf, 0, mefinfo->nf, star, errstr))
	goto error;
    }
  }

  return(0);

 error:
  if(corbuf)
    free((void *) corbuf);

  return(1);
}

