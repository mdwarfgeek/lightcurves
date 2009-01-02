#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>

#include "lightcurves.h"

#include "cvtunit.h"
#include "floatmath.h"
#include "util.h"

#define NITER 3

int lightcurves (struct buffer_info *buf, struct lc_mef *mefinfo,
		 int noapsel, int norenorm, char *errstr) {
  struct lc_point *ptbuf = (struct lc_point *) NULL;
  float *medbuf = (float *) NULL;
  long nmedbuf;

  int iter, degree;

  struct systematic_fit *sysbuf = (struct systematic_fit *) NULL;

  long star, meas, nfluxuse, pt, opt;

  float medflux, sigflux, rmsflux, frameoff, framerms;
  float tmp, chisq;
  long nchisq;

  float medoff, sigoff;

  /* Allocate temporary workspace for calculating medians */
  nmedbuf = MAX(mefinfo->nf, mefinfo->nstars);

  ptbuf = (struct lc_point *) malloc(nmedbuf * sizeof(struct lc_point));
  medbuf = (float *) malloc(nmedbuf * sizeof(float));
  if(!ptbuf || !medbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Allocate buffer for systematics fits */
  sysbuf = (struct systematic_fit *) malloc(mefinfo->nf * sizeof(struct systematic_fit));
  if(!sysbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Loop through the various flux measures */
  nfluxuse = (noapsel ? 1 : NFLUX);

  for(meas = 0; meas < nfluxuse; meas++) {
    /* Apply polynomial correction if requested */
    if(mefinfo->degree >= 0) {
      if(verbose)
	printf("\r Processing aperture %ld of %d",
	       meas+1, NFLUX);
      
      /* Iterate between computing object median and rms, and
       * successive frame corrections.  The polynomial fit
       * (if any) is only done on the last iteration, using
       * a simple constant for the others, so we don't
       * disappear up our own exhaust pipe.
       *
       * The purpose of the earlier iterations is to refine
       * the rms estimate used for weighting in the least
       * squares analysis, so that it reflects the actual
       * error in that object and not the frame-to-frame
       * offsets (which are the dominant effect otherwise
       * on data taken in non-photometric conditions).
       */
      for(iter = 0; iter < NITER; iter++) {
	degree = (iter == NITER-1 ? mefinfo->degree : 0);

	/* Compute per-object median flux */
	for(star = 0; star < mefinfo->nstars; star++) {
	  /* Read in measurements for this star */
	  if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, meas, errstr))
	    goto error;
	  
	  /* Calculate median flux */
	  opt = 0;
	  for(pt = 0; pt < mefinfo->nf; pt++) {
	    if(ptbuf[pt].flux != 0.0) {
	      medbuf[opt] = ptbuf[pt].flux;
	      opt++;
	    }
	  }
	  
	  medsig(medbuf, opt, &medflux, &sigflux);
	  mefinfo->stars[star].medflux[meas] = medflux;
	  mefinfo->stars[star].sigflux[meas] = sigflux;
	}
	
	/* Compute 2-D correction for each frame */
	for(pt = 0; pt < mefinfo->nf; pt++) {
	  /* Read in measurements for this frame */
	  if(buffer_fetch_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	    goto error;
	  
	  /* Perform polynomial fit correction */
	  if(systematic_fit(ptbuf, mefinfo, pt, meas, medbuf, degree, sysbuf+pt,
			    &frameoff, &framerms, errstr))
	    goto error;
	  
	  /* Store frame RMS for normal aperture (meas = 0) */
	  if(meas == 0) {
	    mefinfo->frames[pt].offset = frameoff;
	    mefinfo->frames[pt].rms = framerms;
	    mefinfo->frames[pt].extinc += sysbuf[pt].coeff[0];
	  }
	  
	  /* Perform polynomial fit correction */
	  if(systematic_apply(ptbuf, mefinfo, pt, meas, medbuf, sysbuf, errstr))
	    goto error;
	  
	  /* Write out corrected fluxes */
	  if(buffer_put_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	    goto error;
	}
      }
    }

    if(!noapsel || !norenorm) {
      /* Compute final per-object median flux */
      for(star = 0; star < mefinfo->nstars; star++) {
	/* Read in measurements for this star */
	if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, meas, errstr))
	  goto error;
	
	/* Calculate median flux */
	opt = 0;
	for(pt = 0; pt < mefinfo->nf; pt++) {
	  if(ptbuf[pt].flux != 0.0) {
	    medbuf[opt] = ptbuf[pt].flux;
	    opt++;
	  }
	}
	
	medsig(medbuf, opt, &medflux, &sigflux);
	mefinfo->stars[star].medflux[meas] = medflux;
	mefinfo->stars[star].sigflux[meas] = sigflux;
      }
    }

    if(!norenorm) {
      /* Compute median offset from reference system */
      opt = 0;
      for(star = 0; star < mefinfo->nstars; star++)
	if(mefinfo->stars[star].medflux[meas] != 0.0 &&
	   mefinfo->stars[star].sigflux[meas] != 0.0) {
	  medbuf[opt] = mefinfo->stars[star].medflux[meas] - mefinfo->stars[star].refmag;
	  opt++;
	}
      
      medsig(medbuf, opt, &medoff, &sigoff);
      
      /* Apply offset */
      for(pt = 0; pt < mefinfo->nf; pt++) {
	/* Read in measurements for this frame */
	if(buffer_fetch_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	  goto error;
	
	for(star = 0; star < mefinfo->nstars; star++)
	  ptbuf[star].flux -= medoff;
	
	/* Write out corrected fluxes */
	if(buffer_put_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	  goto error;
      
	if(meas == 0)
	  mefinfo->frames[pt].extinc += medoff;
      }
    }
  }

  if(verbose)
    printf("\n");

  if(!noapsel) {
    /* Combine apertures */
    if(chooseap(buf, mefinfo, ptbuf, medbuf, errstr))
      goto error;
  }

  /* Recompute per-object median flux */
  if(verbose)
    printf(" Computing median fluxes for combined apertures\n");

  for(star = 0; star < mefinfo->nstars; star++) {
    /* Read in measurements for this star */
    if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, 0, errstr))
      goto error;
    
    /* Calculate median flux */
    opt = 0;
    for(pt = 0; pt < mefinfo->nf; pt++) {
      if(ptbuf[pt].flux != 0.0 && ptbuf[pt].fluxerr != 0.0) {
	medbuf[opt] = ptbuf[pt].flux;
	opt++;
      }
    }
    
    if(opt > 0) {
      medsig(medbuf, opt, &medflux, &rmsflux);
      mefinfo->stars[star].med = medflux;

      if(opt > 1) {
	mefinfo->stars[star].rms = rmsflux;
      }
    }
    else {
      mefinfo->stars[star].med = 0.0;
      mefinfo->stars[star].rms = 0.0;
    }
  }

  /* Compute chi-squared */
  if(verbose)
    printf(" Computing chisq statistic\n");

  for(star = 0; star < mefinfo->nstars; star++) {
    /* Read in measurements for this star */
    if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, 0, errstr))
      goto error;

    /* Calculate chi-squared */
    chisq = 0.0;
    nchisq = 0;

    for(pt = 0; pt < mefinfo->nf; pt++) {
      if(ptbuf[pt].flux != 0.0 && ptbuf[pt].fluxerr != 0.0) {
	tmp = ptbuf[pt].flux - mefinfo->stars[star].med;

	chisq += tmp*tmp / (ptbuf[pt].fluxerr * ptbuf[pt].fluxerr);
	nchisq++;
      }
    }

    mefinfo->stars[star].chisq = chisq;
    mefinfo->stars[star].nchisq = nchisq;
  }

  /* Free workspace */
  free((void *) ptbuf);
  ptbuf = (struct lc_point *) NULL;
  free((void *) medbuf);
  medbuf = (float *) NULL;
  free((void *) sysbuf);
  sysbuf = (struct systematic_fit *) NULL;

  return(0);

 error:
  if(!ptbuf)
    free((void *) ptbuf);
  if(!medbuf)
    free((void *) medbuf);
  if(!sysbuf)
    free((void *) sysbuf);

  return(1);
}
