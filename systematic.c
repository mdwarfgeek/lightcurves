#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cpgplot.h>

#include "lightcurves.h"

#include "floatmath.h"
#include "util.h"

/* Parameters: k-sigma for clipping and number of iterations of polynomial fit */
#define SIGCLIP   3
#define NITERMAX  10

/* Clip saturation at 2.5 mags below nominal level from frame */
#define SATCLIP  2.5

/* Use all objects within 4 mags below saturation for fit */
#define USEMAG  4

/* Bin size for plotting */
#define BINSIZE  256  /* pixels */

/* Stop iterations when median and sigma change by less than... */
#define TTHRESH  0.0001  /* 0.1 mmag */

static void polyaccum(float dx, float dy, float wt, float val,
		      double a[50][50], double b[50], int degree) {
  int ox1, oy1, mx1, i1;
  int ox2, oy2, mx2, i2;
  float xp1, yp1, xp2, yp2, tmp;

  i1 = 0;
  for(oy1 = 0, yp1 = 1.0; oy1 <= degree; oy1++, yp1 *= dy) {
    mx1 = degree - oy1;
    for(ox1 = 0, xp1 = 1.0; ox1 <= mx1; ox1++, xp1 *= dx) {
      tmp = xp1 * yp1;

      i2 = 0;
      for(oy2 = 0, yp2 = 1.0; oy2 <= degree; oy2++, yp2 *= dy) {
	mx2 = degree - oy2;
	for(ox2 = 0, xp2 = 1.0; ox2 <= mx2; ox2++, xp2 *= dx) {
	  a[i1][i2] += tmp * xp2 * yp2 * wt;

	  i2++;
	}
      }

      b[i1] += val * tmp * wt;

      i1++;
    }
  }
}

static float polyeval(float dx, float dy, double coeff[50], int degree) {
  float result = 0.0, xp, yp;
  int ox, oy, mx, i;

  i = 0;
  for(oy = 0, yp = 1.0; oy <= degree; oy++, yp *= dy) {
    mx = degree - oy;
    for(ox = 0, xp = 1.0; ox <= mx; ox++, xp *= dx) {
      result += coeff[i] * xp * yp;
      i++;
    }
  }
  
  return(result);
}

int systematic_fit (struct lc_point *data, struct lc_mef *mefinfo, long frame, long meas,
		    float *medbuf, struct systematic_fit *f,
		    float *med_r, float *rms_r, char *errstr) {
  double a[50][50], b[50], coeff[50];
  int ncoeff, ncoeffmax, iter;
  long star, opt;
  float fmin, fmax, siglq;

  int mfirst;

  float medoff, sigoff, newmed, newsig, cxbar, cybar, xbar, ybar, swt;
  float dx, dy, wt, pdx, pdy, val, corr, pval, pcorr;
  int j, k;

  int ilast;

  float xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0, xrange, yrange;

#ifdef DEBUG
  float tmpx[2], tmpy[2], xcord, ycord;
  char title[1024];
#endif

  /* Calculate number of coefficients */
  ncoeffmax = sizeof(coeff) / sizeof(coeff[0]);
  ncoeff = ((mefinfo->degree + 1) * (mefinfo->degree + 2)) / 2;

  if(ncoeff > ncoeffmax) {
    report_err(errstr, "polynomial degree %d too large, maximum is 8", mefinfo->degree);
    goto error;
  }

  /* Calculate the brightest unsaturated magnitude and pixel coordinate
   * ranges for plotting.
   */
  mfirst = 1;
  fmax = 0.0;

  for(star = 0; star < mefinfo->nstars; star++) {
    if(data[star].flux > 0.0 &&          /* Has a flux measurement */
       data[star].fluxerr > 0.0 &&       /* And a reliable error */
       !data[star].satur &&              /* Not saturated */
       !mefinfo->stars[star].bflag &&    /* Not blended */
       !mefinfo->stars[star].cflag &&    /* No bad pixels */
       mefinfo->stars[star].cls == -1) { /* Is classified as stellar */
      if(mfirst || data[star].flux > fmax) {
	fmax = data[star].flux;
	mfirst = 0;
      }
    }

    if(star == 0 || mefinfo->stars[star].x < xmin)
      xmin = mefinfo->stars[star].x;
    if(star == 0 || mefinfo->stars[star].x > xmax)
      xmax = mefinfo->stars[star].x;
    if(star == 0 || mefinfo->stars[star].y < ymin)
      ymin = mefinfo->stars[star].y;
    if(star == 0 || mefinfo->stars[star].y > ymax)
      ymax = mefinfo->stars[star].y;
  }

  /* Allow SATCLIP more mags */
  fmax -= SATCLIP;

  if(mefinfo->syslim >= 0.0)
    /* User override */
    fmax = mefinfo->syslim;

  mefinfo->satclip[meas] = fmax;

  /* Allow all objects between fmax-USEMAG and fmax to counter
   * large numbers of objects with rubbish photometry at the
   * faint end.
   */
  fmin = fmax - USEMAG;

#ifdef DEBUG
  printf("Flux range %.1f %.1f\n", fmin, fmax);
#endif

  /* Calculate lower quartile of sigflux */
  opt = 0;
  for(star = 0; star < mefinfo->nstars; star++) {
    if(data[star].flux > 0.0 &&          /* Has a flux measurement */
       data[star].fluxerr > 0.0 &&       /* And a reliable error */
       mefinfo->stars[star].sigflux[meas] > 0 &&
       data[star].flux > fmin &&
       data[star].flux < fmax &&         /* Not saturated */
       !mefinfo->stars[star].bflag &&    /* Not blended */
       !mefinfo->stars[star].cflag &&    /* No bad pixels */
       mefinfo->stars[star].cls == -1) { /* Is classified as stellar */
      medbuf[opt] = mefinfo->stars[star].sigflux[meas];
      opt++;
    }
  }

  sortfloat(medbuf, opt);

  /* Don't bother if there aren't enough */
  if(opt > 10)
    siglq = medbuf[opt/4];
  else
    siglq = medbuf[opt-1];  /* use everything */

#ifdef DEBUG
  printf("siglq = %f\n", siglq);
#endif

  /* Calculate initial median and sigma offset */
  opt = 0;
  for(star = 0; star < mefinfo->nstars; star++) {
    if(data[star].flux > 0.0 &&          /* Has a flux measurement */
       data[star].fluxerr > 0.0 &&       /* And a reliable error */
       mefinfo->stars[star].sigflux[meas] > 0 &&
       //       mefinfo->stars[star].sigflux[meas] < siglq &&
       data[star].flux > fmin &&
       data[star].flux < fmax &&         /* Not saturated */
       !mefinfo->stars[star].bflag &&    /* Not blended */
       !mefinfo->stars[star].cflag &&    /* No bad pixels */
       mefinfo->stars[star].cls == -1) { /* Is classified as stellar */
      val = data[star].flux - mefinfo->stars[star].medflux[meas];

      medbuf[opt] = val;
      opt++;
    }
  }

  medsig(medbuf, opt, &medoff, &sigoff);

  /* Calculate ranges */
  xrange = xmax - xmin;
  yrange = ymax - ymin;

  /* Initialise master set of coefficients */
  for(k = 0; k < ncoeff; k++)
    coeff[k] = 0.0;

  cxbar = 0.0;
  cybar = 0.0;

  /* Iteratively solve for best-fitting polynomial */
  ilast = 0;

  for(iter = 0; iter < NITERMAX && !ilast; iter++) {
#ifdef DEBUG
    printf("Iteration %d: using median %.4f sigma %.4f nobj %ld\n",
    	   iter+1, medoff, sigoff, opt);
#endif

    /* Initialise matrices */
    for(k = 0; k < ncoeff; k++) {
      for(j = 0; j < ncoeff; j++)
        a[k][j] = 0.0;
      
      b[k] = 0.0;
    }

#ifdef DEBUG
    snprintf(title, sizeof(title),
	     "Frame %ld Iteration %d input objects", frame+1, iter+1);

    cpgenv(xmin - 0.05 * xrange, xmax + 0.05 * xrange,
	   ymin - 0.05 * yrange, ymax + 0.05 * yrange,
	   0, 0);
    cpglab("x", "y", title);
#endif

    /* Calculate means */
    xbar = 0.0;
    ybar = 0.0;
    swt = 0.0;

    for(star = 0; star < mefinfo->nstars; star++) {
      if(data[star].flux > 0.0 &&          /* Has a flux measurement */
	 data[star].fluxerr > 0.0 &&       /* And a reliable error */
	 mefinfo->stars[star].sigflux[meas] > 0 &&
	 //	 mefinfo->stars[star].sigflux[meas] < siglq &&
	 data[star].flux > fmin &&
	 data[star].flux < fmax &&         /* Not saturated */
	 !mefinfo->stars[star].bflag &&    /* Not blended */
	 !mefinfo->stars[star].cflag &&    /* No bad pixels */
	 mefinfo->stars[star].cls == -1) { /* Is classified as stellar */
	pdx = mefinfo->stars[star].x - cxbar;
	pdy = mefinfo->stars[star].y - cybar;
	wt = 1.0 / (mefinfo->stars[star].sigflux[meas] *
		    mefinfo->stars[star].sigflux[meas]);

	pcorr = polyeval(pdx, pdy, coeff, mefinfo->degree);

	val = data[star].flux - mefinfo->stars[star].medflux[meas];
	pval = val - pcorr;

	/* Use the ones within SIGCLIP of the median */
	if(fabsf(pval - medoff) < SIGCLIP * sigoff) {
	  xbar += mefinfo->stars[star].x * wt;
	  ybar += mefinfo->stars[star].y * wt;
	  swt += wt;

#ifdef DEBUG
	  tmpx[0] = mefinfo->stars[star].x;
	  tmpy[0] = mefinfo->stars[star].y;
	  tmpx[1] = mefinfo->stars[star].x;
	  tmpy[1] = mefinfo->stars[star].y + yrange * pval;
	  cpgpt(1, tmpx, tmpy, 22);
	  cpgline(2, tmpx, tmpy);
#endif
	}
      }
    }

    xbar /= swt;
    ybar /= swt;

    /* Accumulate sums for polynomial fit */
    for(star = 0; star < mefinfo->nstars; star++) {
      if(data[star].flux > 0.0 &&          /* Has a flux measurement */
	 data[star].fluxerr > 0.0 &&       /* And a reliable error */
	 mefinfo->stars[star].sigflux[meas] > 0 &&
	 //	 mefinfo->stars[star].sigflux[meas] < siglq &&
	 data[star].flux > fmin &&
	 data[star].flux < fmax &&         /* Not saturated */
	 !mefinfo->stars[star].bflag &&    /* Not blended */
	 !mefinfo->stars[star].cflag &&    /* No bad pixels */
	 mefinfo->stars[star].cls == -1) { /* Is classified as stellar */
	pdx = mefinfo->stars[star].x - cxbar;
	pdy = mefinfo->stars[star].y - cybar;
	dx = mefinfo->stars[star].x - xbar;
	dy = mefinfo->stars[star].y - ybar;
	wt = 1.0 / (mefinfo->stars[star].sigflux[meas] *
		    mefinfo->stars[star].sigflux[meas]);
            
	pcorr = polyeval(pdx, pdy, coeff, mefinfo->degree);

	val = data[star].flux - mefinfo->stars[star].medflux[meas];
	pval = val - pcorr;

	/* Use the ones within SIGCLIP of the median */
	if(fabsf(pval - medoff) < SIGCLIP * sigoff)
	  polyaccum(dx, dy, wt, val, a, b, mefinfo->degree);
      }
    }

    /* Solve for coefficients */
    dsolve(a, b, ncoeff);

#ifdef DEBUG
    /* Plot polynomial surface */
    snprintf(title, sizeof(title),
	     "Frame %ld Iteration %d polynomial", frame+1, iter+1);

    cpgenv(xmin - 0.05 * xrange, xmax + 0.05 * xrange,
	   ymin - 0.05 * yrange, ymax + 0.05 * yrange,
	   0, 0);
    cpglab("x", "y", title);

    for(xcord = xmin; xcord <= xmax; xcord += BINSIZE)
      for(ycord = ymin; ycord <= ymax; ycord += BINSIZE) {
	corr = polyeval(xcord - xbar, ycord - ybar, b, mefinfo->degree);

	tmpx[0] = xcord;
	tmpy[0] = ycord;
	tmpx[1] = xcord;
	tmpy[1] = ycord + yrange * corr;
	cpgpt(1, tmpx, tmpy, 22);
	cpgline(2, tmpx, tmpy);
      }

    snprintf(title, sizeof(title),
	     "Frame %ld Iteration %d corrected", frame+1, iter+1);

    cpgenv(xmin - 0.05 * xrange, xmax + 0.05 * xrange,
	   ymin - 0.05 * yrange, ymax + 0.05 * yrange,
	   0, 0);
    cpglab("x", "y", title);
#endif

    /* Calculate new median and sigma */
    opt = 0;
    for(star = 0; star < mefinfo->nstars; star++) {
      if(data[star].flux > 0.0 &&          /* Has a flux measurement */
	 data[star].fluxerr > 0.0 &&       /* And a reliable error */
	 mefinfo->stars[star].sigflux[meas] > 0 &&
	 //	 mefinfo->stars[star].sigflux[meas] < siglq &&
	 data[star].flux > fmin &&
	 data[star].flux < fmax &&         /* Not saturated */
	 !mefinfo->stars[star].bflag &&    /* Not blended */
	 !mefinfo->stars[star].cflag &&    /* No bad pixels */
	 mefinfo->stars[star].cls == -1) { /* Is classified as stellar */
	pdx = mefinfo->stars[star].x - cxbar;
	pdy = mefinfo->stars[star].y - cybar;
	dx = mefinfo->stars[star].x - xbar;
	dy = mefinfo->stars[star].y - ybar;
	wt = 1.0 / (mefinfo->stars[star].sigflux[meas] *
		    mefinfo->stars[star].sigflux[meas]);
            
	pcorr = polyeval(pdx, pdy, coeff, mefinfo->degree);
	corr = polyeval(dx, dy, b, mefinfo->degree);

	val = data[star].flux - mefinfo->stars[star].medflux[meas];
	pval = val - pcorr;

	/* Use the ones within SIGCLIP of the (old) median */
	if(fabsf(pval - medoff) < SIGCLIP * sigoff) {
	  medbuf[opt] = val - corr;
	  opt++;

#ifdef DEBUG
	  tmpx[0] = mefinfo->stars[star].x;
	  tmpy[0] = mefinfo->stars[star].y;
	  tmpx[1] = mefinfo->stars[star].x;
	  tmpy[1] = mefinfo->stars[star].y + yrange * (val - corr);
	  cpgpt(1, tmpx, tmpy, 22);
	  cpgline(2, tmpx, tmpy);
#endif
	}
      }
    }

    medsig(medbuf, opt, &newmed, &newsig);

    /* Last one? */
    if(fabsf(newmed - medoff) < TTHRESH &&
       fabsf(newsig - sigoff) < TTHRESH)
      ilast = 1;

    /* Store result */
    medoff = newmed;
    sigoff = newsig;
    memcpy(coeff, b, sizeof(coeff));
    cxbar = xbar;
    cybar = ybar;
  }

#ifdef DEBUG
  printf("Final: median %.4f sigma %.4f using %ld objects\n", medoff, sigoff, opt);
#endif

  /* Copy out result */
  f->xbar = cxbar;
  f->ybar = cybar;
  memcpy(f->coeff, coeff, sizeof(f->coeff));

  /* Print coeffs */
  //printf("Coefficients %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n", coeff[0]*1000, coeff[1]*1000, coeff[3]*1000, coeff[2]*1000000, coeff[5]*1000000, coeff[4]*1000000);

  *med_r = medoff;
  *rms_r = sigoff;

  return(0);

 error:
  return(1);
}

int systematic_apply (struct lc_point *data, struct lc_mef *mefinfo, long frame, long meas,
		      float *medbuf, struct systematic_fit *sysbuf, char *errstr) {
  long star, f;
  float dx, dy, corr;

  /* Apply fit */
  for(star = 0; star < mefinfo->nstars; star++) {
    if(data[star].flux > 0.0) {
#if 0
      for(f = 0; f < mefinfo->nf; f++) {
	dx = mefinfo->stars[star].x - sysbuf[f].xbar;
	dy = mefinfo->stars[star].y - sysbuf[f].ybar;

	corr = polyeval(dx, dy, sysbuf[f].coeff, mefinfo->degree);

	medbuf[f] = corr;
      }

      if(hanning(medbuf, mefinfo->nf, errstr))
      	goto error;

      corr = medbuf[frame];
#else
      dx = mefinfo->stars[star].x - sysbuf[frame].xbar;
      dy = mefinfo->stars[star].y - sysbuf[frame].ybar;

      corr = polyeval(dx, dy, sysbuf[frame].coeff, mefinfo->degree);
#endif

      //data[star].flux = 13 + corr;
      data[star].flux -= corr;
    }
  }

  return(0);

 error:
  return(1);
}

