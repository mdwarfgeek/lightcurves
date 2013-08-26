#include <sys/types.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <unistd.h>

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

  int iter, degree;

  struct systematic_fit *sysbuf = (struct systematic_fit *) NULL;

  long star, meas, meas1, meas2, pt, opt1, opt2;
  int used;

  float medflux, sigflux, rmsflux;
  float tmp, chisq;
  long nchisq;

  float medoff, sigoff;

  int iseg, found, haveref;
  float medref, corr;

  float *medcorr = (float *) NULL;
  
  /* Decide how many segments and allocate them */
  mefinfo->segs = (struct lc_segment *) NULL;
  mefinfo->nseg = 0;

  for(pt = 0; pt < mefinfo->nf; pt++) {
    /* Have we seen this segment before? */
    found = -1;
    for(iseg = 0; iseg < mefinfo->nseg; iseg++)
      /* This test is a bit complicated.  The first part handles version
       * numbers, including a NULL = NULL comparison and then the numeric
       * comparison if not null.
       * The second bit handles meridian sides if that is enabled.
       */
      if(((!mefinfo->frames[pt].instvers && !mefinfo->segs[iseg].instvers) ||
	  (mefinfo->frames[pt].instvers && mefinfo->segs[iseg].instvers &&
	   mefinfo->frames[pt].instvers->iver == mefinfo->segs[iseg].instvers->iver)) &&
	 (!mefinfo->domerid || (mefinfo->frames[pt].iang == mefinfo->segs[iseg].iang)))
	found = iseg;

    if(found < 0) {
      /* Add a new one */
      mefinfo->nseg++;
      mefinfo->segs = (struct lc_segment *) realloc(mefinfo->segs,
						    mefinfo->nseg * sizeof(struct lc_segment));
      if(!mefinfo->segs) {
	report_syserr(errstr, "realloc");
	goto error;
      }

      mefinfo->segs[mefinfo->nseg-1].instvers = mefinfo->frames[pt].instvers;
      mefinfo->segs[mefinfo->nseg-1].iang = mefinfo->frames[pt].iang;

      mefinfo->frames[pt].iseg = mefinfo->nseg-1;
    }
    else
      mefinfo->frames[pt].iseg = found;
  }

  for(star = 0; star < mefinfo->nstars; star++) {
    mefinfo->stars[star].segs = (struct lc_star_segment *) malloc(mefinfo->nseg * sizeof(struct lc_star_segment));
    if(!mefinfo->stars[star].segs) {
      report_syserr(errstr, "malloc");
      goto error;
    }

    /* Initialize correction to zero */
    for(iseg = 0; iseg < mefinfo->nseg; iseg++)
      for(meas = 0; meas < NFLUX; meas++)
	mefinfo->stars[star].segs[iseg].corr[meas] = 0;
  }

  /* Allocate temporary workspace for corrections */
  medcorr = (float *) malloc(mefinfo->nseg * sizeof(float));
  if(!medcorr) {
    report_syserr(errstr, "malloc");
    goto error;
  }

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
  meas1 = (mefinfo->aperture ? mefinfo->aperture-1 : 0);
  meas2 = (mefinfo->aperture ? mefinfo->aperture : NFLUX);

  for(meas = meas1; meas < meas2; meas++) {
    /* Apply polynomial correction if requested */
    if(mefinfo->degree >= 0) {
      if(verbose && isatty(1))
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

	/* Compute per-object, per-segment median flux */
	for(star = 0; star < mefinfo->nstars; star++) {
	  /* Read in measurements for this star */
	  if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, meas, errstr))
	    goto error;
	  
	  /* Compute segment medians and global median */
	  opt1 = 0;
	  medref = 0;
	  haveref = 0;

	  for(iseg = 0; iseg < mefinfo->nseg; iseg++) {
	    /* Segment median */
	    opt2 = 0;

	    for(pt = 0; pt < mefinfo->nf; pt++)
	      if(ptbuf[pt].flux != 0.0 && mefinfo->frames[pt].iseg == iseg) {
		medbuf2[opt2] = ptbuf[pt].flux;
		opt2++;
	      }

	    if(opt2 > 0) {
	      medsig(medbuf2, opt2, &medflux, &sigflux);
	      if(iseg == 0) {
		medref = medflux;
		haveref = 1;
	      }

	      if(haveref)
		corr = medflux - medref;
	      else
		corr = 0;

	      if(mefinfo->domerid > 1)
		mefinfo->stars[star].segs[iseg].corr[meas] += corr;
	      else
		mefinfo->stars[star].segs[iseg].corr[meas] = corr;
	    }
	    else
	      corr = 0;

	    /* Global median, after correction */
	    for(pt = 0; pt < mefinfo->nf; pt++)
	      if(ptbuf[pt].flux != 0.0 && mefinfo->frames[pt].iseg == iseg) {
		if(mefinfo->domerid > 1) {
		  ptbuf[pt].flux -= corr;

		  medbuf1[opt1] = ptbuf[pt].flux;
		}
		else
		  medbuf1[opt1] = ptbuf[pt].flux - corr;

		opt1++;
	      }
	  }

	  medsig(medbuf1, opt1, &medflux, &sigflux);
          mefinfo->stars[star].medflux[meas] = medflux;
          mefinfo->stars[star].sigflux[meas] = sigflux;

	  if(mefinfo->domerid > 1)
            /* Write out measurements for this star */
            if(buffer_put_object(buf, ptbuf, 0, mefinfo->nf, star, meas, errstr))
              goto error;
	}

	/* Compute and take off median correction in each segment - this
	 * forces any "common-mode" components present in every object to
	 * go into the frame corrections instead.
	 */
	for(iseg = 0; iseg < mefinfo->nseg; iseg++) {
	  /* Median correction */
	  opt1 = 0;
	  
	  for(star = 0; star < mefinfo->nstars; star++)
	    if(mefinfo->stars[star].medflux[meas] != 0.0 &&
	       mefinfo->stars[star].sigflux[meas] != 0.0)
	      /* Use only high s/n once we have the information */
	      if(iter == 0 ||
		 (mefinfo->stars[star].medflux[meas] >= mefinfo->sysllim &&
		  mefinfo->stars[star].medflux[meas] <= mefinfo->sysulim)) {
		medbuf1[opt1] = mefinfo->stars[star].segs[iseg].corr[meas];
		opt1++;
	      }
	  
	  if(opt1 > 0) {
	    medsig(medbuf1, opt1, &(medcorr[iseg]), (float *) NULL);

	    if(verbose > 2)
	      printf("medcorr[%d]=%.4f\n", iseg, medcorr[iseg]);

	    /* Take it off */
	    for(star = 0; star < mefinfo->nstars; star++)
	      mefinfo->stars[star].segs[iseg].corr[meas] -= medcorr[iseg];
	  }
	  else
	    medcorr[iseg] = 0;
	}
	
	/* Compute 2-D correction for each frame */
	for(pt = 0; pt < mefinfo->nf; pt++) {
	  /* Read in measurements for this frame */
	  if(buffer_fetch_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	    goto error;

	  /* Take off median correction from actual points if necessary */
	  if(mefinfo->domerid > 1)
	    for(star = 0; star < mefinfo->nstars; star++)
	      if(ptbuf[star].flux != 0.0)
		ptbuf[star].flux += medcorr[mefinfo->frames[pt].iseg];
	  
	  /* Perform polynomial fit correction */
	  if(systematic_fit(ptbuf, mefinfo, pt, meas, medbuf1, degree,
			    sysbuf+pt, errstr))
	    goto error;
	  
	  /* Trap for no points in fit - discard in this case */
	  if(sysbuf[pt].npt > 0) {
	    /* Store frame RMS for normal aperture */
	    if(meas == meas1) {
	      mefinfo->frames[pt].offset = sysbuf[pt].medoff;
	      mefinfo->frames[pt].rms = sysbuf[pt].sigoff;
	      mefinfo->frames[pt].extinc += sysbuf[pt].coeff[0];
	      mefinfo->frames[pt].sigm = sqrtf(sysbuf[pt].cov[0][0]);
	    }

	    /* Perform polynomial fit correction */
	    if(systematic_apply(ptbuf, mefinfo, pt, meas, medbuf1, sysbuf, errstr))
	      goto error;
	  }	    

	  if(sysbuf[pt].npt > 0 || mefinfo->domerid > 1) {
	    /* Write out corrected fluxes */
	    if(buffer_put_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	      goto error;
	  }
	}
      }
    }

    /* Compute final per-object median flux */
    for(star = 0; star < mefinfo->nstars; star++) {
      /* Read in measurements for this star */
      if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, meas, errstr))
	goto error;
	
      /* Calculate median flux and set "used" flag for comparison stars */
      opt1 = 0;
      used = 0;
      for(pt = 0; pt < mefinfo->nf; pt++) {
	if(ptbuf[pt].flux != 0.0) {
	  medbuf1[opt1] = ptbuf[pt].flux;
	  opt1++;
	}
	if(ptbuf[pt].wt > 0)
	  used = 1;
      }
      
      medsig(medbuf1, opt1, &medflux, &sigflux);
      mefinfo->stars[star].medflux[meas] = medflux;
      mefinfo->stars[star].sigflux[meas] = sigflux;
      mefinfo->stars[star].used += used;
    }

    if(!norenorm) {
      /* Compute median offset from reference system */
      opt1 = 0;
      for(star = 0; star < mefinfo->nstars; star++)
	if(mefinfo->stars[star].medflux[meas] != 0.0 &&
	   mefinfo->stars[star].sigflux[meas] != 0.0 &&
	   mefinfo->stars[star].medflux[meas] >= mefinfo->sysllim &&
	   mefinfo->stars[star].medflux[meas] <= mefinfo->sysulim) {
	  medbuf1[opt1] = mefinfo->stars[star].medflux[meas] - mefinfo->stars[star].refmag;
	  opt1++;
	}
      
      medsig(medbuf1, opt1, &medoff, &sigoff);
      
      if(verbose)
	printf("  Calibration offset: %.3f %.3f\n", medoff, sigoff);

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
      
	if(meas == meas1)
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

  if(verbose && isatty(1))
    printf("\n");

  /* Aperture correct and choose default aperture */
  if(chooseap(buf, mefinfo, ptbuf, medbuf1, errstr))
    goto error;

  /* Compute median and chi squared in chosen aperture */
  if(verbose)
    printf(" Computing median fluxes and chisq for chosen apertures\n");

  for(star = 0; star < mefinfo->nstars; star++) {
    /* Read in measurements for this star */
    if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, mefinfo->stars[star].iap, errstr))
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

    /* Calculate median x,y positions */
    for(iseg = 0; iseg < mefinfo->nseg; iseg++) {
      opt1 = 0;

      for(pt = 0; pt < mefinfo->nf; pt++) {
	if(ptbuf[pt].flux != 0.0 && mefinfo->frames[pt].iseg == iseg) {
	  medbuf1[opt1] = ptbuf[pt].x;
	  medbuf2[opt1] = ptbuf[pt].y;
	  opt1++;
	}
      }
      
      if(opt1 > 0) {
	medsig(medbuf1, opt1,
	       &(mefinfo->stars[star].segs[iseg].medx),
	       &(mefinfo->stars[star].segs[iseg].sigx));
	medsig(medbuf2, opt1,
	       &(mefinfo->stars[star].segs[iseg].medy),
	       &(mefinfo->stars[star].segs[iseg].sigy));
      }
    }
  }

  if(verbose)
    printf(" Recomputing astrometric transformations\n");

  /* Compute positioning error for each frame */
  for(pt = 0; pt < mefinfo->nf; pt++) {
    /* Read in measurements for this frame */
    if(buffer_fetch_frame(buf, ptbuf, 0, mefinfo->nstars, pt, 0, errstr))
      goto error;

    /* Compute transformation to ref. system */
    if(xytoxy(ptbuf, mefinfo, mefinfo->frames[pt].tr, errstr))
      goto error;

    /* OLD: */
    opt1 = 0;

    for(star = 0; star < mefinfo->nstars; star++)
      if(ptbuf[star].flux != 0.0) {
	medbuf1[opt1] = ptbuf[star].x - mefinfo->stars[star].segs[mefinfo->frames[pt].iseg].medx;
	medbuf2[opt1] = ptbuf[star].y - mefinfo->stars[star].segs[mefinfo->frames[pt].iseg].medy;
	opt1++;
      }
  
    if(opt1 > 0) {
      medsig(medbuf1, opt1, &(mefinfo->frames[pt].xoff), &(mefinfo->frames[pt].xsig));
      medsig(medbuf2, opt1, &(mefinfo->frames[pt].yoff), &(mefinfo->frames[pt].ysig));
    }
  }

  /* Compute median positioning errors (should be zero) and rms */
  for(iseg = 0; iseg < mefinfo->nseg; iseg++) {
    opt1 = 0;
    
    for(pt = 0; pt < mefinfo->nf; pt++) {
      if(mefinfo->frames[pt].iseg == iseg) {
	medbuf1[opt1] = mefinfo->frames[pt].xoff;
	medbuf2[opt1] = mefinfo->frames[pt].yoff;
	opt1++;
      }
    }
    
    medsig(medbuf1, opt1, &(mefinfo->segs[iseg].medxoff), &(mefinfo->segs[iseg].sigxoff));
    medsig(medbuf2, opt1, &(mefinfo->segs[iseg].medyoff), &(mefinfo->segs[iseg].sigyoff));
  }

  /* Free workspace */
  free((void *) ptbuf);
  ptbuf = (struct lc_point *) NULL;
  free((void *) medbuf1);
  medbuf1 = (float *) NULL;
  free((void *) sysbuf);
  sysbuf = (struct systematic_fit *) NULL;
  free((void *) medcorr);
  medcorr = (float *) NULL;

  return(0);

 error:
  if(ptbuf)
    free((void *) ptbuf);
  if(medbuf1)
    free((void *) medbuf1);
  if(sysbuf)
    free((void *) sysbuf);
  if(medcorr)
    free((void *) medcorr);

  return(1);
}

int lightcurves_append (struct buffer_info *buf, struct lc_mef *mefinfo,
			char *errstr) {
  struct lc_point *ptbuf = (struct lc_point *) NULL;
  float *medbuf = (float *) NULL;
  long nmedbuf;

  struct systematic_fit *sysbuf = (struct systematic_fit *) NULL;

  long star, meas, meas1, meas2, pt, opt;
  float medflux, rmsflux, chisq, tmp;
  long nchisq;

  int aper, useaper;
  float diff, diffmin;

  int iseg, found, imr;

  /* Decide segments for all new points */
  for(pt = 0; pt < mefinfo->nf; pt++) {
    /* Have we seen this segment before? */
    found = -1;
    imr = mefinfo->nseg-1;
    for(iseg = 0; iseg < mefinfo->nseg; iseg++) {
      /* This test is a bit complicated.  The first part handles version
       * numbers, including a NULL = NULL comparison and then the numeric
       * comparison if not null.
       * The second bit handles meridian sides if that is enabled.
       */
      if(((!mefinfo->frames[pt].instvers && !mefinfo->segs[iseg].instvers) ||
	  (mefinfo->frames[pt].instvers && mefinfo->segs[iseg].instvers &&
	   mefinfo->frames[pt].instvers->iver == mefinfo->segs[iseg].instvers->iver)) &&
	 (!mefinfo->domerid || (mefinfo->frames[pt].iang == mefinfo->segs[iseg].iang)))
	found = iseg;

      /* Accumulate most recent on correct meridian side */
      if(!mefinfo->domerid || (mefinfo->frames[pt].iang == mefinfo->segs[iseg].iang))
	imr = iseg;
    }

    if(found < 0)
      found = imr;  /* assume most recent */

    mefinfo->frames[pt].iseg = found;
  }

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
  meas1 = (mefinfo->aperture ? mefinfo->aperture-1 : 0);
  meas2 = (mefinfo->aperture ? mefinfo->aperture : NFLUX);

  for(meas = meas1; meas < meas2; meas++) {
    /* Apply polynomial correction if requested */
    if(mefinfo->degree >= 0) {
      if(verbose && isatty(1))
	printf("\r Processing aperture %ld of %d",
	       meas+1, NFLUX);
      
      /* Compute 2-D correction for each frame */
      for(pt = 0; pt < mefinfo->nf; pt++) {
	/* Read in measurements for this frame */
	if(buffer_fetch_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	  goto error;
	
	/* Apply any meridian flip offsets */
	if(mefinfo->domerid > 1) 
	  for(star = 0; star < mefinfo->nstars; star++) 
	    if(ptbuf[star].flux != 0.0)
	      ptbuf[star].flux -= mefinfo->stars[star].segs[mefinfo->frames[pt].iseg].corr[meas];

	/* Perform polynomial fit correction */
	if(systematic_fit(ptbuf, mefinfo, pt, meas, medbuf, mefinfo->degree,
			  sysbuf+pt, errstr))
	  goto error;
	
	/* Trap for no points in fit - discard in this case */
	if(sysbuf[pt].npt > 0) {
	  /* Store frame RMS for normal aperture (meas = 0) */
	  if(meas == meas1) {
	    mefinfo->frames[pt].offset = sysbuf[pt].medoff;
	    mefinfo->frames[pt].rms = sysbuf[pt].sigoff;
	    mefinfo->frames[pt].extinc += sysbuf[pt].coeff[0];
	    mefinfo->frames[pt].sigm = sqrtf(sysbuf[pt].cov[0][0]);
	  }
	  
	  /* Perform polynomial fit correction */
	  if(systematic_apply(ptbuf, mefinfo, pt, meas, medbuf, sysbuf, errstr))
	    goto error;
	}	  
	
	if(sysbuf[pt].npt > 0 || mefinfo->domerid > 1) {
	  /* Write out corrected fluxes */
	  if(buffer_put_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	    goto error;
	}
      }
    }
  }

  if(verbose && isatty(1))
    printf("\n");

  if(mefinfo->aperture == 0) {
    /* Apply old aperture selections */
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

      mefinfo->stars[star].iap = useaper;
    }
  }
  else {
    for(star = 0; star < mefinfo->nstars; star++)
      mefinfo->stars[star].iap = mefinfo->aperture-1;
  }

  /* Recompute per-object median flux (NEW POINTS ONLY) */
  if(verbose)
    printf(" Computing median fluxes\n");

  for(star = 0; star < mefinfo->nstars; star++) {
    /* Read in measurements for this star */
    if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, mefinfo->stars[star].iap, errstr))
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
    if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, mefinfo->stars[star].iap, errstr))
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

  if(verbose)
    printf(" Recomputing astrometric transformations\n");

  /* Compute positioning error for each frame */
  for(pt = 0; pt < mefinfo->nf; pt++) {
    /* Read in measurements for this frame */
    if(buffer_fetch_frame(buf, ptbuf, 0, mefinfo->nstars, pt, 0, errstr))
      goto error;

    /* Compute transformation to ref. system */
    if(xytoxy(ptbuf, mefinfo, mefinfo->frames[pt].tr, errstr))
      goto error;
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
  if(ptbuf)
    free((void *) ptbuf);
  if(medbuf)
    free((void *) medbuf);
  if(sysbuf)
    free((void *) sysbuf);

  return(1);
}
