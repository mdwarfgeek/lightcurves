#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef DEBUG
#include <cpgplot.h>
#endif

#include "lightcurves.h"

#include "util.h"

/* Parameters: k-sigma for clipping and number of iterations of polynomial fit */
#define SIGCLIP   3
#define NITERMAX  10

/* Clip saturation at 2.5 mags below nominal level from frame */
#define SATCLIP  2.5

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

static float polyvar(float dx, float dy, double cov[50][50], int degree) {
  int ox1, oy1, mx1, i1;
  int ox2, oy2, mx2, i2;
  float xp1, yp1, xp2, yp2, tmp, result;

  result = 0;

  i1 = 0;
  for(oy1 = 0, yp1 = 1.0; oy1 <= degree; oy1++, yp1 *= dy) {
    mx1 = degree - oy1;
    for(ox1 = 0, xp1 = 1.0; ox1 <= mx1; ox1++, xp1 *= dx) {
      tmp = xp1 * yp1;

      i2 = 0;
      for(oy2 = 0, yp2 = 1.0; oy2 <= degree; oy2++, yp2 *= dy) {
	mx2 = degree - oy2;
	for(ox2 = 0, xp2 = 1.0; ox2 <= mx2; ox2++, xp2 *= dx) {
	  result += cov[i1][i2] * tmp * xp2 * yp2;

	  i2++;
	}
      }

      i1++;
    }
  }
  
  return(result);
}

/* This routine flags potential reference stars.  Called from read_lc
   or read_ref in readfits.c when the star list is read in.  This
   allows unused stars to be stripped out at the earliest possible
   stage in read_ref for performance reasons, when the "-c" command
   line option is given. */

void systematic_select (struct lc_star *stars, long nstars, float fmin, float fmax) {
  long star;

  if(fmax >= 0 && fmin < 0)
    fmin = fmax - USEMAG;

  for(star = 0; star < nstars; star++) {
    if(stars[star].ref.aper[REFAP].flux > 0 &&
       stars[star].ref.aper[REFAP].fluxvar > 0 &&
       stars[star].cls == -1 &&
       (fmin < 0 || stars[star].refmag > fmin) &&
       (fmax < 0 || stars[star].refmag < fmax))
      stars[star].compok = 1;
    else
      stars[star].compok = 0;
  }
}

int systematic_fit (struct lc_point *data, struct lc_mef *mefinfo, long frame, long meas,
		    float *medbuf, int degree, struct systematic_fit *f,
		    char *errstr) {
  double a[50][50], ainv[50][50], b[50], coeff[50], cov[50][50];
  int ncoeff, ncoeffmax, iter;
  long star, opt;
  float fmin, fmax;
  float rmsclip;

  int mfirst;

  float medoff, sigoff, newmed, newsig, cxbar, cybar, xbar, ybar, swt;
  float dx, dy, wt, pdx, pdy, val, corr, pval, pcorr;
  int j, k;

  int ilast;

  float xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;
  float lastvar;  /* kludge to capture sigma of last star used for npt=1 case */

  float chisq = 0.0, newchisq, varscale;

#ifdef DEBUG
  float xrange, yrange;
  float tmpx[2], tmpy[2], xcord, ycord;
  char title[1024];
#endif

  /* Calculate number of coefficients */
  ncoeffmax = sizeof(coeff) / sizeof(coeff[0]);
  ncoeff = ((degree + 1) * (degree + 2)) / 2;

  if(ncoeff > ncoeffmax) {
    report_err(errstr, "polynomial degree %d too large, maximum is 8", degree);
    goto error;
  }

  /* Calculate the brightest unsaturated magnitude and pixel coordinate
   * ranges for plotting.
   */
  mfirst = 1;
  fmax = 0.0;

  for(star = 0; star < mefinfo->nstars; star++) {
    if(data[star].aper[meas].flux > 0.0 &&          /* Has a flux measurement */
       data[star].aper[meas].fluxvar > 0.0 &&       /* And a reliable error */
       !data[star].satur &&                         /* Not saturated */
       mefinfo->stars[star].compok) {               /* OK for comp */
      if(mfirst || mefinfo->stars[star].refmag > fmax) {
	fmax = mefinfo->stars[star].refmag;
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

  if(mefinfo->sysulim >= 0.0)
    /* User override */
    fmax = mefinfo->sysulim;
  else
    /* Store it */
    mefinfo->sysulim = fmax;

  mefinfo->satclip[meas] = fmax;

  /* Allow all objects between fmax-USEMAG and fmax to counter
   * large numbers of objects with rubbish photometry at the
   * faint end.
   */
  fmin = fmax - USEMAG;

  if(mefinfo->sysllim >= 0.0)
    /* User override */
    fmin = mefinfo->sysllim;
  else
    /* Store it */
    mefinfo->sysllim = fmin;
  
#ifdef DEBUG
  printf("Flux range %.1f %.1f\n", fmin, fmax);
#endif

  /* Calculate median and sigma of the rms of the comp stars */
  opt = 0;
  for(star = 0; star < mefinfo->nstars; star++) {
    if(data[star].aper[meas].flux > 0.0 &&          /* Has a flux measurement */
       data[star].aper[meas].fluxvar > 0.0 &&       /* And a reliable error */
       mefinfo->stars[star].sigflux[meas] > 0 &&
       mefinfo->stars[star].refmag > fmin &&
       mefinfo->stars[star].refmag < fmax &&        /* Right mag range */
       mefinfo->stars[star].compok) {               /* OK for comp */
      medbuf[opt] = mefinfo->stars[star].sigflux[meas];
      opt++;
    }
  }

  /* 90th percentile.  Uses ceil so it behaves the same as the original
     did (given what we use this for below), without needing to average
     two elements. */
  if(opt > 5)
    rmsclip = fquickselect(medbuf, ceil(0.9*opt), opt);
  else
    rmsclip = 999.0;  /* I think this should be safe :) */

  /* Calculate initial median and sigma offset */
  opt = 0;
  for(star = 0; star < mefinfo->nstars; star++) {
    if(data[star].aper[meas].flux > 0.0 &&          /* Has a flux measurement */
       data[star].aper[meas].fluxvar > 0.0 &&       /* And a reliable error */
       mefinfo->stars[star].sigflux[meas] > 0 &&
       mefinfo->stars[star].sigflux[meas] < rmsclip &&
       mefinfo->stars[star].refmag > fmin &&
       mefinfo->stars[star].refmag < fmax &&        /* Right mag range */
       mefinfo->stars[star].compok) {               /* OK for comp */
      val = data[star].aper[meas].flux - mefinfo->stars[star].medflux[meas];

      medbuf[opt] = val;
      opt++;
    }
  }

  fmedsig(medbuf, opt, &medoff, &sigoff);

#ifdef DEBUG
  /* Calculate ranges */
  xrange = xmax - xmin;
  yrange = ymax - ymin;
#endif

  /* Initialise master set of coefficients and covariance */
  for(k = 0; k < ncoeff; k++) {
    coeff[k] = 0.0;

    for(j = 0; j < ncoeff; j++)
      cov[k][j] = 0.0;
  }

  cxbar = 0.0;
  cybar = 0.0;

  /* Iteratively solve for best-fitting polynomial */
  ilast = 0;

  lastvar = 0;

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
      if(data[star].aper[meas].flux > 0.0 &&          /* Has a flux measurement */
	 data[star].aper[meas].fluxvar > 0.0 &&       /* And a reliable error */
	 mefinfo->stars[star].sigflux[meas] > 0 &&
	 mefinfo->stars[star].sigflux[meas] < rmsclip &&
	 mefinfo->stars[star].refmag > fmin &&
	 mefinfo->stars[star].refmag < fmax &&        /* Right mag range */
	 mefinfo->stars[star].compok) {               /* OK for comp */
	pdx = mefinfo->stars[star].x - cxbar;
	pdy = mefinfo->stars[star].y - cybar;
	wt = 1.0 / (mefinfo->stars[star].sigflux[meas] *
		    mefinfo->stars[star].sigflux[meas]);

	pcorr = polyeval(pdx, pdy, coeff, degree);

	val = data[star].aper[meas].flux - mefinfo->stars[star].medflux[meas];
	pval = val - pcorr;

	/* Use the ones within SIGCLIP of the median */
	if(sigoff == 0.0 || fabsf(pval - medoff) < SIGCLIP * sigoff) {
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
      if(data[star].aper[meas].flux > 0.0 &&          /* Has a flux measurement */
	 data[star].aper[meas].fluxvar > 0.0 &&       /* And a reliable error */
	 mefinfo->stars[star].sigflux[meas] > 0 &&
	 mefinfo->stars[star].sigflux[meas] < rmsclip &&
	 mefinfo->stars[star].refmag > fmin &&
	 mefinfo->stars[star].refmag < fmax &&        /* Right mag range */
	 mefinfo->stars[star].compok) {               /* OK for comp */
	pdx = mefinfo->stars[star].x - cxbar;
	pdy = mefinfo->stars[star].y - cybar;
	dx = mefinfo->stars[star].x - xbar;
	dy = mefinfo->stars[star].y - ybar;
	wt = 1.0 / (mefinfo->stars[star].sigflux[meas] *
		    mefinfo->stars[star].sigflux[meas]);
            
	pcorr = polyeval(pdx, pdy, coeff, degree);

	val = data[star].aper[meas].flux - mefinfo->stars[star].medflux[meas];
	pval = val - pcorr;

	/* Use the ones within SIGCLIP of the median */
	if(sigoff == 0.0 || fabsf(pval - medoff) < SIGCLIP * sigoff) {
	  polyaccum(dx, dy, wt, val, a, b, degree);
	  data[star].aper[meas].wt = wt;
	}
	else
	  data[star].aper[meas].wt = 0;
      }
      else
	data[star].aper[meas].wt = 0;
    }

    /* Make a copy first */
    memcpy(ainv, a, sizeof(ainv));

    /* Solve for coefficients */
    dsolve(a, b, ncoeff);

    /* Invert to covariance matrix */
    dmatinv(ainv, ncoeff);

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
	corr = polyeval(xcord - xbar, ycord - ybar, b, degree);

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
    newchisq = 0;
    opt = 0;
    for(star = 0; star < mefinfo->nstars; star++) {
      if(data[star].aper[meas].flux > 0.0 &&          /* Has a flux measurement */
	 data[star].aper[meas].fluxvar > 0.0 &&       /* And a reliable error */
	 mefinfo->stars[star].sigflux[meas] > 0 &&
	 mefinfo->stars[star].sigflux[meas] < rmsclip &&
	 mefinfo->stars[star].refmag > fmin &&
	 mefinfo->stars[star].refmag < fmax &&        /* Right mag range */
	 mefinfo->stars[star].compok) {               /* OK for comp */
	pdx = mefinfo->stars[star].x - cxbar;
	pdy = mefinfo->stars[star].y - cybar;
	dx = mefinfo->stars[star].x - xbar;
	dy = mefinfo->stars[star].y - ybar;
	wt = 1.0 / (mefinfo->stars[star].sigflux[meas] *
		    mefinfo->stars[star].sigflux[meas]);

	pcorr = polyeval(pdx, pdy, coeff, degree);
	corr = polyeval(dx, dy, b, degree);

	val = data[star].aper[meas].flux - mefinfo->stars[star].medflux[meas];
	pval = val - pcorr;

	/* Use the ones within SIGCLIP of the (old) median */
	if(sigoff == 0.0 || fabsf(pval - medoff) < SIGCLIP * sigoff) {
	  medbuf[opt] = val - corr;
	  newchisq += (val - corr)*(val - corr) * wt;
	  opt++;

	  lastvar = data[star].aper[meas].fluxvar;
	
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

    fmedsig(medbuf, opt, &newmed, &newsig);

    /* Last one? */
    if(fabsf(newmed - medoff) < TTHRESH &&
       fabsf(newsig - sigoff) < TTHRESH)
      ilast = 1;

    /* Store result */
    medoff = newmed;
    sigoff = newsig;
    memcpy(coeff, b, sizeof(coeff));
    memcpy(cov, ainv, sizeof(cov));
    cxbar = xbar;
    cybar = ybar;
    chisq = newchisq;
  }

#ifdef DEBUG
  printf("Final: median %.4f sigma %.4f using %ld objects\n", medoff, sigoff, opt);
#endif

  /* Correct covariance using chi^2/dof */
  varscale = (opt > ncoeff ? chisq/(opt-ncoeff) : 1.0);

  for(k = 0; k < ncoeff; k++)
    for(j = 0; j < ncoeff; j++)
      cov[k][j] *= varscale;

  /* Copy out result */
  f->xbar = cxbar;
  f->ybar = cybar;
  memcpy(f->coeff, coeff, sizeof(f->coeff));
  memcpy(f->cov, cov, sizeof(f->cov));
  f->degree = degree;
  f->medoff = medoff;
  f->sigoff = sigoff;
  f->sigm = (opt > 1 ? (opt > ncoeff ? sigoff / sqrt(opt-ncoeff) : 0.0) : sqrtf(lastvar));
  f->npt = opt;

  /* Print coeffs */
  //printf("Coefficients %9.5f %9.5f x %9.5f y %9.5f x**2 %9.5f y**2 %9.5f x*y\n", coeff[0]*1000, coeff[1]*1000, coeff[3]*1000, coeff[2]*1000000, coeff[5]*1000000, coeff[4]*1000000);

  return(0);

 error:
  return(1);
}

int systematic_apply (struct lc_point *data, struct lc_mef *mefinfo, long frame, long meas,
		      float *medbuf, struct systematic_fit *sysbuf, char *errstr) {
  long star;
  float dx, dy, corr, var;

  /* Apply fit */
  for(star = 0; star < mefinfo->nstars; star++) {
    if(data[star].aper[meas].flux > 0.0) {
      dx = mefinfo->stars[star].x - sysbuf[frame].xbar;
      dy = mefinfo->stars[star].y - sysbuf[frame].ybar;

      corr = polyeval(dx, dy, sysbuf[frame].coeff, sysbuf[frame].degree);
      var = polyvar(dx, dy, sysbuf[frame].cov, sysbuf[frame].degree);

      //printf("%f %f %f %f\n", corr, sqrt(var), sysbuf[frame].sigoff, sysbuf[frame].sigm);

      data[star].aper[meas].flux -= corr;

      if(data[star].aper[meas].fluxvar > 0)
	data[star].aper[meas].fluxvarcom = data[star].aper[meas].fluxvar + var;
      else
	data[star].aper[meas].fluxvarcom = 0.0;
    }
  }

  return(0);
}

