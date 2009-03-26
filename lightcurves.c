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
		 int norenorm, char *errstr) {
  struct lc_point *ptbuf = (struct lc_point *) NULL;
  float *medbuf1 = (float *) NULL, *medbuf2;
  long nmedbuf;

  float med1, med2, corr;

  int iter, degree;

  struct systematic_fit *sysbuf = (struct systematic_fit *) NULL;

  long star, meas, nfluxuse, pt, opt1, opt2;

  float medflux, sigflux, rmsflux, frameoff, framerms;
  float tmp, chisq;
  long nchisq;
  long framenpt;

  float medoff, sigoff;
  float framesigm;

  /* Allocate temporary workspace for calculating medians */
  nmedbuf = MAX(mefinfo->nf, mefinfo->nstars);

  ptbuf = (struct lc_point *) malloc(nmedbuf * sizeof(struct lc_point));
  medbuf1 = (float *) malloc(2 * nmedbuf * sizeof(float));
  if(!ptbuf || !medbuf1) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  medbuf2 = medbuf1 + nmedbuf;

  /* Allocate buffer for systematics fits */
  sysbuf = (struct systematic_fit *) malloc(mefinfo->nf * sizeof(struct systematic_fit));
  if(!sysbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  /* Loop through the various flux measures */
  nfluxuse = (mefinfo->doapsel ? NFLUX : 1);

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
	  opt1 = 0;
	  opt2 = 0;
	  for(pt = 0; pt < mefinfo->nf; pt++) {
	    if(ptbuf[pt].flux != 0.0) {
	      if(mefinfo->domerid && ptbuf[pt].ha > 0.0) {
		medbuf2[opt2] = ptbuf[pt].flux;
		opt2++;
	      }
	      else {
		medbuf1[opt1] = ptbuf[pt].flux;
		opt1++;
	      }
	    }
	  }

	  if(mefinfo->domerid && opt1 > 0 && opt2 > 0) {
	    medsig(medbuf1, opt1, &med1, (float *) NULL);
	    medsig(medbuf2, opt2, &med2, (float *) NULL);

	    corr = 0.5*(med2-med1);

	    opt1 = 0;
	    for(pt = 0; pt < mefinfo->nf; pt++) {
	      if(ptbuf[pt].flux != 0.0) {
		if(ptbuf[pt].ha > 0.0)
		  medbuf1[opt1] = ptbuf[pt].flux - corr;
		else
		  medbuf1[opt1] = ptbuf[pt].flux + corr;

		opt1++;
	      }
	    }
	  }
	  else {
	    opt1 = 0;
	    for(pt = 0; pt < mefinfo->nf; pt++)
	      if(ptbuf[pt].flux != 0.0) {
		medbuf1[opt1] = ptbuf[pt].flux;
		opt1++;
	      }
	    }

	  medsig(medbuf1, opt1, &medflux, &sigflux);
	  mefinfo->stars[star].medflux[meas] = medflux;
	  mefinfo->stars[star].sigflux[meas] = sigflux;
	}
	
	/* Compute 2-D correction for each frame */
	for(pt = 0; pt < mefinfo->nf; pt++) {
	  /* Read in measurements for this frame */
	  if(buffer_fetch_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	    goto error;
	  
	  /* Perform polynomial fit correction */
	  if(systematic_fit(ptbuf, mefinfo, pt, meas, medbuf1, degree, sysbuf+pt,
			    &frameoff, &framerms, &framesigm, &framenpt, errstr))
	    goto error;
	  
	  /* Trap for no points in fit - discard in this case */
	  if(framenpt > 0) {
	    /* Store frame RMS for normal aperture (meas = 0) */
	    if(meas == 0) {
	      mefinfo->frames[pt].offset = frameoff;
	      mefinfo->frames[pt].rms = framerms;
	      mefinfo->frames[pt].extinc += sysbuf[pt].coeff[0];
	      mefinfo->frames[pt].sigm = framesigm;
	    }

	    /* Perform polynomial fit correction */
	    if(systematic_apply(ptbuf, mefinfo, pt, meas, medbuf1, sysbuf,
				framesigm, errstr))
	      goto error;
	    
	    /* Write out corrected fluxes */
	    if(buffer_put_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	      goto error;
	  }
	}
      }
    }

    if(mefinfo->doapsel || !norenorm) {
      /* Compute final per-object median flux */
      for(star = 0; star < mefinfo->nstars; star++) {
	/* Read in measurements for this star */
	if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, meas, errstr))
	  goto error;
	
	/* Calculate median flux */
	opt1 = 0;
	for(pt = 0; pt < mefinfo->nf; pt++) {
	  if(ptbuf[pt].flux != 0.0) {
	    medbuf1[opt1] = ptbuf[pt].flux;
	    opt1++;
	  }
	}
	
	medsig(medbuf1, opt1, &medflux, &sigflux);
	mefinfo->stars[star].medflux[meas] = medflux;
	mefinfo->stars[star].sigflux[meas] = sigflux;
      }
    }

    if(!norenorm) {
      /* Compute median offset from reference system */
      opt1 = 0;
      for(star = 0; star < mefinfo->nstars; star++)
	if(mefinfo->stars[star].medflux[meas] != 0.0 &&
	   mefinfo->stars[star].sigflux[meas] != 0.0) {
	  medbuf1[opt1] = mefinfo->stars[star].medflux[meas] - mefinfo->stars[star].refmag;
	  opt1++;
	}
      
      medsig(medbuf1, opt1, &medoff, &sigoff);
      
      /* Apply offset */
      for(pt = 0; pt < mefinfo->nf; pt++) {
	/* Read in measurements for this frame */
	if(buffer_fetch_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	  goto error;
	
	for(star = 0; star < mefinfo->nstars; star++)
	  if(ptbuf[star].flux != 0.0)
	    ptbuf[star].flux -= medoff;
	
	/* Write out corrected fluxes */
	if(buffer_put_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	  goto error;
      
	if(meas == 0)
	  mefinfo->frames[pt].extinc += medoff;
      }

      /* Compute final per-object median flux */
      for(star = 0; star < mefinfo->nstars; star++) {
	/* Read in measurements for this star */
	if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, meas, errstr))
	  goto error;
	
	/* Calculate median flux */
	opt1 = 0;
	for(pt = 0; pt < mefinfo->nf; pt++) {
	  if(ptbuf[pt].flux != 0.0) {
	    medbuf1[opt1] = ptbuf[pt].flux;
	    opt1++;
	  }
	}
	
	medsig(medbuf1, opt1, &medflux, &sigflux);
	mefinfo->stars[star].medflux[meas] = medflux;
	mefinfo->stars[star].sigflux[meas] = sigflux;
      }
    }
  }

  if(verbose)
    printf("\n");

  if(mefinfo->doapsel) {
    /* Combine apertures */
    if(chooseap(buf, mefinfo, ptbuf, medbuf1, errstr))
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
    opt1 = 0;
    for(pt = 0; pt < mefinfo->nf; pt++) {
      if(ptbuf[pt].flux != 0.0 && ptbuf[pt].fluxerrcom != 0.0) {
	medbuf1[opt1] = ptbuf[pt].flux;
	opt1++;
      }
    }
    
    if(opt1 > 0) {
      medsig(medbuf1, opt1, &medflux, &rmsflux);
      mefinfo->stars[star].med = medflux;

      if(opt1 > 1) {
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
      if(ptbuf[pt].flux != 0.0 && ptbuf[pt].fluxerrcom != 0.0) {
	tmp = ptbuf[pt].flux - mefinfo->stars[star].med;

	chisq += tmp*tmp / (ptbuf[pt].fluxerrcom * ptbuf[pt].fluxerrcom);
	nchisq++;
      }
    }

    mefinfo->stars[star].chisq = chisq;
    mefinfo->stars[star].nchisq = nchisq;
  }

  /* Free workspace */
  free((void *) ptbuf);
  ptbuf = (struct lc_point *) NULL;
  free((void *) medbuf1);
  medbuf1 = (float *) NULL;
  free((void *) sysbuf);
  sysbuf = (struct systematic_fit *) NULL;

  return(0);

 error:
  if(!ptbuf)
    free((void *) ptbuf);
  if(!medbuf1)
    free((void *) medbuf1);
  if(!sysbuf)
    free((void *) sysbuf);

  return(1);
}

int lightcurves_append (struct buffer_info *buf, struct lc_mef *mefinfo,
			char *errstr) {
  struct lc_point *ptbuf = (struct lc_point *) NULL;
  float *medbuf = (float *) NULL;
  long nmedbuf;

  struct systematic_fit *sysbuf = (struct systematic_fit *) NULL;

  long star, meas, nfluxuse, pt, opt;
  float medflux, rmsflux, chisq, tmp;
  long nchisq;
  long framenpt;

  float frameoff, framerms;

  int aper, useaper;
  float diff, diffmin;

  float framesigm;

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
  nfluxuse = (mefinfo->doapsel ? NFLUX : 1);

  for(meas = 0; meas < nfluxuse; meas++) {
    /* Apply polynomial correction if requested */
    if(mefinfo->degree >= 0) {
      if(verbose)
	printf("\r Processing aperture %ld of %d",
	       meas+1, NFLUX);
      
      /* Compute 2-D correction for each frame */
      for(pt = 0; pt < mefinfo->nf; pt++) {
	/* Read in measurements for this frame */
	if(buffer_fetch_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	  goto error;
	
	/* Perform polynomial fit correction */
	if(systematic_fit(ptbuf, mefinfo, pt, meas, medbuf, mefinfo->degree, sysbuf+pt,
			  &frameoff, &framerms, &framesigm, &framenpt, errstr))
	  goto error;
	
	/* Trap for no points in fit - discard in this case */
	if(framenpt > 0) {
	  /* Store frame RMS for normal aperture (meas = 0) */
	  if(meas == 0) {
	    mefinfo->frames[pt].offset = frameoff;
	    mefinfo->frames[pt].rms = framerms;
	    mefinfo->frames[pt].extinc += sysbuf[pt].coeff[0];
	    mefinfo->frames[pt].sigm = framesigm;
	  }
	  
	  /* Perform polynomial fit correction */
	  if(systematic_apply(ptbuf, mefinfo, pt, meas, medbuf, sysbuf,
			      framesigm, errstr))
	    goto error;
	  
	  /* Write out corrected fluxes */
	  if(buffer_put_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	    goto error;
	}
      }
    }
  }

  if(verbose)
    printf("\n");

  /* Apply old aperture selections */
  if(mefinfo->doapsel) {
    /* Loop through all stars */
    for(star = 0; star < mefinfo->nstars; star++) {
      /* Which one was it? */
      useaper = -1;
      diffmin = 0;
      for(aper = 0; aper < NFLUX; aper++) {
	diff = fabsf(mefinfo->stars[star].apradius - flux_apers[aper]);
	
	if(useaper < 0 || diff < diffmin) {
	  useaper = aper;
	  diffmin = diff;
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
	
	/* Write out into aperture 0 */
	if(buffer_put_object(buf, ptbuf, 0, mefinfo->nf, star, 0, errstr))
	  goto error;
      }
    }
  }

  /* Recompute per-object median flux (NEW POINTS ONLY) */
  if(verbose)
    printf(" Computing median fluxes for combined apertures\n");

  for(star = 0; star < mefinfo->nstars; star++) {
    /* Read in measurements for this star */
    if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, 0, errstr))
      goto error;
    
    /* Calculate median flux */
    opt = 0;
    for(pt = 0; pt < mefinfo->nf; pt++) {
      if(ptbuf[pt].flux != 0.0 && ptbuf[pt].fluxerrcom != 0.0) {
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

  /* Compute chi-squared (NEW POINTS ONLY) */
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
      if(ptbuf[pt].flux != 0.0 && ptbuf[pt].fluxerrcom != 0.0) {
	tmp = ptbuf[pt].flux - mefinfo->stars[star].med;

	chisq += tmp*tmp / (ptbuf[pt].fluxerrcom * ptbuf[pt].fluxerrcom);
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
