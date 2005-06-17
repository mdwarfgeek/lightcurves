#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>

#include <cpgplot.h>

#include "lightcurves.h"

#include "cvtunit.h"
#include "floatmath.h"
#include "util.h"

/* Sampling for theoretical curve */
#define THEOSTEP  0.01  /* mag */

int do_plots (struct lc_mef *meflist, int nmefs,
	      float medsat, float medlim, float sysbodge, char *errstr) {
  float magmin, magmax, mag, rms, chi, photons, loge, area, tmp;
  long star, pt;
  int mef;

  float *medbuf1 = (float *) NULL, *medbuf2, *medbuf3, *medbuf4;
  float gain, rcore, sigma, avzp;

  float *theox = (float *) NULL, *theoy, *theoys;
  long t, ntheo;

  float tmpx[2], tmpy[2];

  char title[1024];
  float ptmin, ptmax, tpt, ty;

  /* Open PGPLOT */
  cpgopen("?");
  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);

  /* Allocate workspace for median parameters */
  medbuf1 = (float *) malloc(4 * nmefs * sizeof(float));
  if(!medbuf1) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  medbuf2 = medbuf1 + nmefs;
  medbuf3 = medbuf1 + 2 * nmefs;
  medbuf4 = medbuf1 + 3 * nmefs;

  /* RMS plot */
  magmin = medsat;
  magmax = medlim;

  cpgenv(magmin, magmax, -0.1, 2.9, 0, 20);
  cpglab("Magnitude", "RMS (millimag)", "");

  for(mef = 0; mef < nmefs; mef++) {
    for(star = 0; star < meflist[mef].nstars; star++) {
      /* Use only the ones classified as stars */
      if(meflist[mef].stars[star].medflux[0] > 0.0 &&
	 meflist[mef].stars[star].rms > 0.0 &&
	 meflist[mef].stars[star].cls == -1) {
	mag = meflist[mef].zp - meflist[mef].stars[star].medflux[0];
	rms = 3.0 + log10f(meflist[mef].stars[star].rms);

	cpgpt(1, &mag, &rms, 1);
      }
    }

    medbuf1[mef] = meflist[mef].refgain;
    medbuf2[mef] = meflist[mef].refrcore;
    medbuf3[mef] = meflist[mef].avsigma;
    medbuf4[mef] = meflist[mef].zp;
  }

  medsig(medbuf1, nmefs, &gain, (float *) NULL);
  medsig(medbuf2, nmefs, &rcore, (float *) NULL);
  medsig(medbuf3, nmefs, &sigma, (float *) NULL);
  medsig(medbuf4, nmefs, &avzp, (float *) NULL);

  free((void *) medbuf1);
  medbuf1 = (float *) NULL;

  /* Generate theoretical curve */
  ntheo = ceil((magmax - magmin) / THEOSTEP);

  theox = (float *) malloc(3 * ntheo * sizeof(float));
  if(!theox) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  theoy = theox + ntheo;
  theoys = theox + 2 * ntheo;

  loge = log10f(expf(1.0));
  area = M_PI * rcore * rcore;

  for(t = 0; t < ntheo; t++) {
    mag = t * 0.01 + magmin;

    photons = powf(10.0, 0.4 * (avzp - mag)) * gain;
    tmp = 2.5 * loge * sqrtf(photons +
			     area * (gain * gain * sigma * sigma)) / photons;

    theox[t] = mag;
    theoy[t] = 3.0 + log10f(tmp);
    theoys[t] = 3.0 + 0.5 * log10f(tmp*tmp + sysbodge*sysbodge);
  }

  cpgsls(2);
  cpgsci(2);
  cpgline(ntheo, theox, theoy);
  cpgsci(1);

  if(sysbodge > 0.0) {
    cpgsls(4);
    cpgsci(2);
    cpgline(ntheo, theox, theoys);
    cpgsci(1);
  }

  free((void *) theox);
  theox = (float *) NULL;

  tmpx[0] = medsat;
  tmpx[1] = medsat;
  tmpy[0] = -0.1;
  tmpy[1] = 2.9;
  cpgline(2, tmpx, tmpy);

  tmpx[0] = medlim;
  tmpx[1] = medlim;
  tmpy[0] = -0.1;
  tmpy[1] = 2.9;
  cpgline(2, tmpx, tmpy);

  cpgsls(1);

  /* Chi squared plot */
  cpgenv(magmin, magmax, -1.0, 4.0, 0, 20);
  cpglab("Magnitude", "\\gx\\u2\\d\\dred\\u", "");

  for(mef = 0; mef < nmefs; mef++)
    for(star = 0; star < meflist[mef].nstars; star++) {
      /* Use only the ones classified as stars */
      if(meflist[mef].stars[star].medflux[0] > 0.0 &&
	 meflist[mef].stars[star].chisq > 0.0 &&
	 meflist[mef].stars[star].nchisq > 1 &&
	 meflist[mef].stars[star].cls == -1) {
	mag = meflist[mef].zp - meflist[mef].stars[star].medflux[0];
	chi = log10f(meflist[mef].stars[star].chisq /
		     (meflist[mef].stars[star].nchisq - 1));
	
	cpgpt(1, &mag, &chi, 1);
      }
    }

  cpgsls(1);

  /* Frame corrections and RMS */
  cpgsubp(1, 2);

  for(mef = 0; mef < nmefs; mef++) {
    snprintf(title, sizeof(title), "MEF %d", mef+1);

    ptmin = 0.0;
    ptmax = meflist[mef].nf+1;

    cpgsch(1.4);

    cpgenv(ptmin, ptmax, -0.1, 0.1, 0, 0);
    cpglab("Sample", "Frame offset", title);

    cpgsch(2.0);
    for(pt = 0; pt < meflist[mef].nf; pt++) {
      tpt = pt + 1;
      ty = meflist[mef].frames[pt].offset;

      cpgpt(1, &tpt, &ty, 17);
    }

    cpgsci(2);
    for(pt = 2; pt < meflist[mef].nf; pt++) 
      if(meflist[mef].frames[pt].mjd - meflist[mef].frames[pt-1].mjd > 0.5) {
	tmpx[0] = pt + 0.5;
	tmpx[1] = tmpx[0];
	tmpy[0] = -0.1;
	tmpy[1] = 0.1;
	cpgline(2, tmpx, tmpy);
      }
    cpgsci(1);

    cpgsch(1.4);

    cpgenv(ptmin, ptmax, -0.01, 0.10, 0, 0);
    cpglab("Sample", "Frame RMS", "");

    cpgsch(2.0);
    for(pt = 0; pt < meflist[mef].nf; pt++) {
      tpt = pt + 1;
      ty = meflist[mef].frames[pt].rms;

      cpgpt(1, &tpt, &ty, 17);
    }

    cpgsci(2);
    for(pt = 2; pt < meflist[mef].nf; pt++) 
      if(meflist[mef].frames[pt].mjd - meflist[mef].frames[pt-1].mjd > 0.5) {
	tmpx[0] = pt + 0.5;
	tmpx[1] = tmpx[0];
	tmpy[0] = -0.01;
	tmpy[1] = 0.10;
	cpgline(2, tmpx, tmpy);
      }
    cpgsci(1);

    cpgsls(2);
    tmpx[0] = ptmin;
    tmpx[1] = ptmax;
    tmpy[0] = 0.0;
    tmpy[1] = 0.0;
    cpgline(2, tmpx, tmpy);
    cpgsls(1);

    cpgsch(1.0);
  }

  cpgclos();

  return(0);

 error:
  if(medbuf1)
    free((void *) medbuf1);
  if(theox)
    free((void *) theox);

  return(1);
}

int plot_corr (float *beforehist, float *beforewthist,
	       float *corrhist, float *corrwthist,
	       float *afterhist, float *afterwthist,
	       long nbinx, long nbiny, float binsize,
	       float xmin, float xmax, float ymin, float ymax,
	       float medoff, float sigoff, char *errstr) {
  long bin, nbin;
  float vmin = 0.0, vmax = 0.0, tr[6];

  /* Normalise histograms */
  nbin = nbinx * nbiny;

  for(bin = 0; bin < nbin; bin++) {
    if(beforewthist[bin] > 0.0)
      beforehist[bin] /= beforewthist[bin];
    else
      beforehist[bin] = medoff;

    if(corrwthist[bin] > 0.0)
      corrhist[bin] /= corrwthist[bin];
    else
      corrhist[bin] = medoff;

    if(afterwthist[bin] > 0.0)
      afterhist[bin] /= afterwthist[bin];
    else
      afterhist[bin] = medoff;

    if(bin == 0 || corrhist[bin] < vmin)
      vmin = corrhist[bin];
    if(bin == 0 || corrhist[bin] > vmax)
      vmax = corrhist[bin];
  }

  /* Calculate transformation matrix */
  tr[0] = -binsize / 2.0;
  tr[1] = binsize;
  tr[2] = 0.0;
  tr[3] = -binsize / 2.0;
  tr[4] = 0.0;
  tr[5] = binsize;

  /* Plot */
  cpgenv(xmin, xmax, ymin, ymax, 0, 0);
  cpglab("x", "y", "Before correction");
  cpggray(beforehist, nbinx, nbiny, 1, nbinx, 1, nbiny, vmax, vmin, tr);

  cpgenv(xmin, xmax, ymin, ymax, 0, 0);
  cpglab("x", "y", "Polynomial correction");
  cpggray(corrhist, nbinx, nbiny, 1, nbinx, 1, nbiny, vmax, vmin, tr);

  cpgenv(xmin, xmax, ymin, ymax, 0, 0);
  cpglab("x", "y", "After correction");
  cpggray(afterhist, nbinx, nbiny, 1, nbinx, 1, nbiny, vmax, vmin, tr);

  cpgpage();

  return(0);
}


