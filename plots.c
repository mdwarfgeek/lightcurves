#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>

#include <cpgplot.h>

#include "lightcurves.h"

#include "cvtunit.h"
#include "util.h"

/* Sampling for theoretical curve */
#define THEOSTEP  0.01  /* mag */

int do_plots (struct lc_mef *meflist, int nmefs,
	      float medsat, float medlim, float umlim, float lmlim, char *errstr) {
  float magmin, magmax, lrmin, lrmax;
  float mag, rms, chi, photons, skyvar, area, tmp, tmpp, tmps;
  long star, pt;
  int mef;

  float *medbuf1 = (float *) NULL, *medbuf2, *medbuf3, *medbuf4;
  float *medbuf5, *medbuf6, *medbuf7, *medbuf8, *medbuf9;
  float gain, rcore, sigma, skyfiterr, avzp, avapcor, avextinc, avsigm, avscint;

  float *theox = (float *) NULL, *theop, *theos, *theoy, *theoys;
  long t, ntheo;

  float tmpx[2], tmpy[2];

  char title[1024], xlab[1024];
  float ptmin, ptmax, tpt, ty;

  /* Open PGPLOT */
  cpgopen("?");
  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);
  cpgsch(1.4);

  /* Allocate workspace for median parameters */
  medbuf1 = (float *) malloc(9 * nmefs * sizeof(float));
  if(!medbuf1) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  medbuf2 = medbuf1 + nmefs;
  medbuf3 = medbuf1 + 2 * nmefs;
  medbuf4 = medbuf1 + 3 * nmefs;
  medbuf5 = medbuf1 + 4 * nmefs;
  medbuf6 = medbuf1 + 5 * nmefs;
  medbuf7 = medbuf1 + 6 * nmefs;
  medbuf8 = medbuf1 + 7 * nmefs;
  medbuf9 = medbuf1 + 8 * nmefs;

  for(mef = 0; mef < nmefs; mef++) {
    medbuf1[mef] = meflist[mef].refgain;
    medbuf2[mef] = meflist[mef].refrcore;
    medbuf3[mef] = meflist[mef].avsigma;
    medbuf4[mef] = meflist[mef].avskyfit;
    medbuf5[mef] = meflist[mef].zp;
    medbuf6[mef] = meflist[mef].avapcor;
    medbuf7[mef] = meflist[mef].avextinc;
    medbuf8[mef] = meflist[mef].avsigm;
    medbuf9[mef] = meflist[mef].avscint;
  }

  medsig(medbuf1, nmefs, &gain, (float *) NULL);
  medsig(medbuf2, nmefs, &rcore, (float *) NULL);
  medsig(medbuf3, nmefs, &sigma, (float *) NULL);
  medsig(medbuf4, nmefs, &skyfiterr, (float *) NULL);
  medsig(medbuf5, nmefs, &avzp, (float *) NULL);
  medsig(medbuf6, nmefs, &avapcor, (float *) NULL);
  medsig(medbuf7, nmefs, &avextinc, (float *) NULL);
  medsig(medbuf8, nmefs, &avsigm, (float *) NULL);
  medsig(medbuf9, nmefs, &avscint, (float *) NULL);

  free((void *) medbuf1);
  medbuf1 = (float *) NULL;

  /* Plot x range */
  magmin = MIN(medsat, umlim);
  magmax = medlim;

  /* Apply average aperture and extinction correction to ZP */
  avzp -= 2.5*log10(avapcor*avextinc);

  /* Generate theoretical curve */
  ntheo = ceil((magmax - magmin) / THEOSTEP);

  theox = (float *) malloc(5 * ntheo * sizeof(float));
  if(!theox) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  theop = theox + ntheo;
  theos = theox + 2 * ntheo;
  theoy = theox + 3 * ntheo;
  theoys = theox + 4 * ntheo;

  area = M_PI * rcore * rcore;

  if(verbose)
    printf("avzp=%g gain=%g sigma=%g skyfiterr=%g area=%g avextinc=%g avsigm=%g avscint=%g\n",
	   avzp, gain, sigma, skyfiterr, area, avextinc, avsigm, avscint);

  for(t = 0; t < ntheo; t++) {
    mag = t * THEOSTEP + magmin;

    photons = powf(10.0, 0.4 * (avzp - mag)) * gain;
    skyvar = gain*gain*area * (sigma*sigma +
                               area*skyfiterr*skyfiterr);

    tmp = 2.5 * M_LOG10E * sqrtf(photons + skyvar) / photons;
    tmpp = 2.5 * M_LOG10E * sqrtf(photons) / photons;
    tmps = 2.5 * M_LOG10E * sqrtf(skyvar) / photons;

    theox[t] = mag;
    theop[t] = 3.0 + log10f(tmpp);
    theos[t] = 3.0 + log10f(tmps);
    theoy[t] = 3.0 + log10f(tmp);
    theoys[t] = 3.0 + 0.5 * log10f(tmp*tmp + avsigm*avsigm + avscint*avscint);
  }

  /* log(RMS) range */
  lrmin = theoys[0];
  lrmax = theoys[0];

  for(t = 0; t < ntheo; t++) {
    if(theoys[t] < lrmin)
      lrmin = theoys[t];
    if(theoys[t] > lrmax)
      lrmax = theoys[t];
  }

  lrmin -= 0.05*(lrmax-lrmin);
  lrmax += 0.05*(lrmax-lrmin);

  lrmin = floor(lrmin);
  if(lrmin > 0)
    lrmin = 0;

  lrmax = ceil(lrmax);
  if(lrmax < lrmin + 3)
    lrmax = lrmin + 3;

  /* RMS plot */
  snprintf(xlab, sizeof(xlab), "%s magnitude", meflist[0].filter);

  cpgenv(magmin, magmax, lrmin, lrmax-1.0e-3, 0, 20);
  cpglab(xlab, "RMS (millimag)", "");

  for(mef = 0; mef < nmefs; mef++)
    for(star = 0; star < meflist[mef].nstars; star++) {
      /* Use only the ones classified as stars */
      if(meflist[mef].stars[star].med > 0.0 &&
	 meflist[mef].stars[star].rms > 0.0) {
	mag = meflist[mef].zp - meflist[mef].stars[star].med;
	rms = 3.0 + log10f(meflist[mef].stars[star].rms);

	cpgsci(1+mef);

	if(meflist[mef].stars[star].cls == -1 &&
	   meflist[mef].stars[star].bflag == 0 &&
	   meflist[mef].stars[star].cflag == 0) {  /* BODGE */
	  cpgpt(1, &mag, &rms, 1);
	}
	else if(meflist[mef].stars[star].cls == 9) {
	  cpgpt(1, &mag, &rms, 1);
	  cpgsci(2);
	  cpgpt(1, &mag, &rms, 22);
	}

	cpgsci(1);
      }
    }

  cpgsci(2);
  cpgslw(4);
  cpgsls(1);
  cpgline(ntheo, theox, theoys);
  cpgsls(2);
  cpgline(ntheo, theox, theop);
  cpgsls(3);
  cpgline(ntheo, theox, theos);

  cpgsls(4);
  tmpx[0] = magmin;
  tmpx[1] = magmax;
  tmpy[0] = 3.0 + log10f(avsigm);
  tmpy[1] = tmpy[0];
  cpgline(2, tmpx, tmpy);
  cpgsls(1);

  if(avscint > 0.0) {
    cpgsls(5);
    tmpx[0] = magmin;
    tmpx[1] = magmax;
    tmpy[0] = 3.0 + log10f(avscint);
    tmpy[1] = tmpy[0];
    cpgline(2, tmpx, tmpy);
    cpgsls(1);
  }
  cpgslw(1);
  cpgsci(1);

  free((void *) theox);
  theox = (float *) NULL;

  tmpx[0] = medsat;
  tmpx[1] = medsat;
  tmpy[0] = lrmin;
  tmpy[1] = lrmax;
  cpgline(2, tmpx, tmpy);

  cpgsls(2);

  tmpx[0] = umlim;
  tmpx[1] = umlim;
  tmpy[0] = lrmin;
  tmpy[1] = lrmax;
  cpgline(2, tmpx, tmpy);

  tmpx[0] = lmlim;
  tmpx[1] = lmlim;
  tmpy[0] = lrmin;
  tmpy[1] = lrmax;
  cpgline(2, tmpx, tmpy);

  cpgsls(1);

  tmpx[0] = medlim;
  tmpx[1] = medlim;
  tmpy[0] = lrmin;
  tmpy[1] = lrmax;
  cpgline(2, tmpx, tmpy);

  /* Chi squared plot */
  cpgenv(magmin, magmax, -1.0, 4.0, 0, 20);
  cpglab(xlab, "\\gx\\u2\\d\\dred\\u", "");

  for(mef = 0; mef < nmefs; mef++)
    for(star = 0; star < meflist[mef].nstars; star++) {
      /* Use only the ones classified as stars */
      if(meflist[mef].stars[star].med > 0.0 &&
	 meflist[mef].stars[star].chisq > 0.0 &&
	 meflist[mef].stars[star].nchisq > 1) {
	mag = meflist[mef].zp - meflist[mef].stars[star].med;
	chi = log10f(meflist[mef].stars[star].chisq /
		     (meflist[mef].stars[star].nchisq - 1));


	cpgsci(1+mef);

	if(meflist[mef].stars[star].cls == -1)
	  cpgpt(1, &mag, &chi, 1);
	else if(meflist[mef].stars[star].cls == 9) {
	  cpgpt(1, &mag, &chi, 1);
	  cpgsci(2);
	  cpgpt(1, &mag, &chi, 22);
	}

	cpgsci(1);
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


