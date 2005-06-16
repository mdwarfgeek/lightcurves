#include <stdio.h>
#include <stdlib.h>

#include <cpgplot.h>

#include "lightcurves.h"

#include "floatmath.h"
#include "util.h"

/* Maximum number of components to fit and remove */
#define CITERMAX  5

/* Maximum number of iterations per component */
#define NITERMAX  100

/* Convergence criterion on coefficients */
#define CONVERGED 1.0e-4

int pcasys (struct buffer_info *buf, struct lc_point *ptbuf, struct lc_mef *mefinfo,
	    long meas, float *medbuf, char *errstr) {
  int citer, iter;
  long star, pt, opt;
  float fval, wt, sdc, scc, cval, aval, cdiff, adiff, meancorr;
  float ctmp, atmp;

  float *medflux = (float *) NULL, *sigflux;
  float *ccoeff = (float *) NULL, *acoeff;

  float xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0, xrange, yrange;

#ifdef DEBUG
  float tmpx[2], tmpy[2];
  char title[1024];
#endif

  /* Calculate pixel ranges */
  for(star = 0; star < mefinfo->nstars; star++) {
    if(star == 0 || mefinfo->stars[star].x < xmin)
      xmin = mefinfo->stars[star].x;
    if(star == 0 || mefinfo->stars[star].x > xmax)
      xmax = mefinfo->stars[star].x;
    if(star == 0 || mefinfo->stars[star].y < ymin)
      ymin = mefinfo->stars[star].y;
    if(star == 0 || mefinfo->stars[star].y > ymax)
      ymax = mefinfo->stars[star].y;
  }

  xrange = xmax - xmin;
  yrange = ymax - ymin;

  /* Allocate buffers */
  medflux = (float *) malloc(2 * mefinfo->nstars * sizeof(float));
  ccoeff = (float *) malloc((mefinfo->nstars + mefinfo->nf) * sizeof(float));
  if(!medflux || !ccoeff) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  sigflux = medflux + mefinfo->nstars;
  acoeff = ccoeff + mefinfo->nstars;

  /* Calculate per-object median flux */
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
    
    medsig(medbuf, opt, &(medflux[star]), &(sigflux[star]));
  }

  /* Loop through iterations of successive systematic components */
  for(citer = 0; citer < CITERMAX; citer++) {
    /* Initialise coefficients */
    for(star = 0; star < mefinfo->nstars; star++)
      ccoeff[star] = 0.0;
    for(pt = 0; pt < mefinfo->nf; pt++)
      acoeff[pt] = 1.0;

    /* Calculate per-object median flux */
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
      
      if(opt > 1)
	medsig(medbuf, opt, &(medflux[star]), &(sigflux[star]));
      else
	sigflux[star] = 0.0;
    }

    /* Loop through iterations */
    for(iter = 0; iter < NITERMAX; iter++) {
      /* Compute c for each star */
      cdiff = 0.0;

      for(star = 0; star < mefinfo->nstars; star++) {
	/* Read in measurements for this star */
	if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, meas, errstr))
	  goto error;

	/* Calculate weight */
	if(sigflux[star] > 0.0) {
	  wt = 1.0 / (sigflux[star]*sigflux[star]);
	  
	  /* Accumulate sums */
	  sdc = 0.0;
	  scc = 0.0;
	  for(pt = 0; pt < mefinfo->nf; pt++) {
	    if(ptbuf[pt].flux != 0.0) {
	      fval = ptbuf[pt].flux - medflux[star];
	      
	      sdc += fval * acoeff[pt] * wt;
	      scc += acoeff[pt] * acoeff[pt] * wt;
	    }
	  }
	  
	  if(scc != 0.0)
	    cval = sdc / scc;
	  else
	    cval = 0.0;
	}
	else
	  cval = 0.0;

	ctmp = cval - ccoeff[star];
	cdiff += ctmp*ctmp;
	ccoeff[star] = cval;
      }

      cdiff = sqrtf(cdiff / mefinfo->nstars);

      /* Now compute a for each frame */
      adiff = 0.0;
      meancorr = 0.0;

      for(pt = 0; pt < mefinfo->nf; pt++) {
	/* Read in measurements for this frame */
	if(buffer_fetch_frame(buf, ptbuf, 0, mefinfo->nstars, pt, meas, errstr))
	  goto error;

	/* Accumulate sums */
	sdc = 0.0;
	scc = 0.0;
	for(star = 0; star < mefinfo->nstars; star++) {
	  if(ptbuf[star].flux != 0.0 && sigflux[star] > 0.0) {
	    fval = ptbuf[star].flux - medflux[star];
	    wt = 1.0 / (sigflux[star]*sigflux[star]);
	    
	    sdc += fval * ccoeff[star] * wt;
	    scc += ccoeff[star] * ccoeff[star] * wt;
	  }
	}

	if(scc != 0.0)
	  aval = sdc / scc;
	else
	  aval = 0.0;

	atmp = aval - acoeff[pt];
	adiff += atmp*atmp;
	acoeff[pt] = aval;

	/* Compute mean correction over everything */
	for(star = 0; star < mefinfo->nstars; star++)
	  meancorr += aval * ccoeff[star];
      }

      adiff = sqrtf(adiff / mefinfo->nf);
      meancorr /= (mefinfo->nf * mefinfo->nstars);

      printf("Component %d iteration %d adiff %.4f cdiff %.4f meancorr %.4f\n",
	     citer+1, iter+1, adiff, cdiff, meancorr);

      if(iter > 0 && adiff < CONVERGED && cdiff < CONVERGED)
	break;
    }

    /* Give up if not any good */
    if(adiff >= CONVERGED || cdiff >= CONVERGED)
      break;

#ifdef DEBUG
    /* Plot before */
    snprintf(title, sizeof(title),
	     "Component %d input objects", citer+1);

    cpgenv(xmin - 0.05 * xrange, xmax + 0.05 * xrange,
	   ymin - 0.05 * yrange, ymax + 0.05 * yrange,
	   0, 0);
    cpglab("x", "y", title);

    for(star = 0; star < mefinfo->nstars; star++) {
      if(sigflux[star] > 0.0 &&
	 !mefinfo->stars[star].bflag &&
	 mefinfo->stars[star].cls == -1) {
	tmpx[0] = mefinfo->stars[star].x;
	tmpy[0] = mefinfo->stars[star].y;
	tmpx[1] = mefinfo->stars[star].x;
	tmpy[1] = mefinfo->stars[star].y + yrange * sigflux[star];
	cpgpt(1, tmpx, tmpy, 22);
	cpgline(2, tmpx, tmpy);
      }
    }

    /* Plot after */
    snprintf(title, sizeof(title),
	     "Component %d corrected", citer+1);

    cpgenv(xmin - 0.05 * xrange, xmax + 0.05 * xrange,
	   ymin - 0.05 * yrange, ymax + 0.05 * yrange,
	   0, 0);
    cpglab("x", "y", title);
#endif

    /* Remove them and re-calculate median flux */
    for(star = 0; star < mefinfo->nstars; star++) {
      /* Read in measurements for this star */
      if(buffer_fetch_object(buf, ptbuf, 0, mefinfo->nf, star, meas, errstr))
	goto error;
      
      /* Remove and calculate median flux */
      opt = 0;
      for(pt = 0; pt < mefinfo->nf; pt++) {
	if(ptbuf[pt].flux != 0.0) {
	  ptbuf[pt].flux -= acoeff[pt]*ccoeff[star];

	  medbuf[opt] = ptbuf[pt].flux;
	  opt++;
	}
      }
      
      medsig(medbuf, opt, &(medflux[star]), &(sigflux[star]));

#ifdef DEBUG
      /* Plot */
      if(sigflux[star] > 0.0 &&
	 !mefinfo->stars[star].bflag &&
	 mefinfo->stars[star].cls == -1) {
	tmpx[0] = mefinfo->stars[star].x;
	tmpy[0] = mefinfo->stars[star].y;
	tmpx[1] = mefinfo->stars[star].x;
	tmpy[1] = mefinfo->stars[star].y + yrange * sigflux[star];
	cpgpt(1, tmpx, tmpy, 22);
	cpgline(2, tmpx, tmpy);
      }
#endif

      /* Write out */
      if(buffer_put_object(buf, ptbuf, 0, mefinfo->nf, star, meas, errstr))
	goto error;
    }
  }

  free((void *) medflux);
  medflux = (float *) NULL;
  free((void *) ccoeff);
  ccoeff = (float *) NULL;

  return(0);

 error:
  if(medflux)
    free((void *) medflux);
  if(ccoeff)
    free((void *) ccoeff);

  return(1);
}
