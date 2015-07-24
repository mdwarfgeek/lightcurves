#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>

#include "lightcurves.h"

#ifdef HJD
#include "hjd.h"
#endif

#include "cvtunit.h"
#include "util.h"

int output_init (struct lc_output *op,
                 struct lc_mef *mefinfo,
                 struct dtai_table *dtab,
                 struct iers_table *itab,
                 struct jpleph_table *jtab,
                 struct jpleph_table *ttab,
                 char *errstr) {
  long pt;

  int utcdn;
  double mjd, utcdf, dtt;

  /* Init these */
  op->obs = (struct observer *) NULL;
#ifdef HJD
  op->epos = (double *) NULL;
#endif

  /* Figure out apertures */
  if(mefinfo->aperture) {  /* only requested aperture */
    op->ap1 = mefinfo->aperture-1;
    op->ap2 = mefinfo->aperture;
    op->napcol = 1;
    op->nlapcol = 1;
  }
  else {
    op->ap1 = 0;
    op->ap2 = NFLUX;
    op->napcol = op->ap2-op->ap1 + 1;

    if(mefinfo->apselmode & APSEL_ALL)  /* write all */
      op->nlapcol = op->napcol;
    else  /* only chosen */
      op->nlapcol = 1;
  }

  /* Allocate buffers */
  op->obs = (struct observer *) malloc(mefinfo->nf * sizeof(struct observer));
  if(!op->obs) {
    report_syserr(errstr, "malloc");
    goto error;
  }

#ifdef HJD
  op->epos = (double *) malloc(3 * mefinfo->nf * sizeof(double));
  if(!op->epos) {
    report_syserr(errstr, "malloc");
    goto error;
  }
#endif

  for(pt = 0; pt < mefinfo->nf; pt++) {
    mjd = mefinfo->mjdref + mefinfo->frames[pt].mjd;

    if(mefinfo->frames[pt].doairm && jtab) {
      /* Create observer structure */
      observer_init(op->obs+pt,
                    mefinfo->frames[pt].longitude,
                    mefinfo->frames[pt].latitude,
                    mefinfo->frames[pt].height);
      
      memcpy(op->obs[pt].refco, mefinfo->frames[pt].refco, sizeof(op->obs[pt].refco));
    }
    else
      observer_geoc(op->obs+pt);  /* geocentric corrections only */

    if(dtab) {
      /* Look up delta(TT) */
      utcdn = floor(mjd);
      utcdf = mjd - utcdn;
      
      dtt = dtai(dtab, utcdn, utcdf) + DTT;
    }
    else
      dtt = DTT;  /* assumes TAI=UTC */
    
    if(jtab) {
      /* Compute time-dependent quantities */
      observer_update(op->obs+pt, jtab, ttab, itab, mjd, dtt, OBSERVER_UPDATE_ALL);
    }

#ifdef HJD
    /* Calculate Earth's heliocentric position at this MJD */
    getearth(mjd, op->epos + 3*pt);
#endif    
  }

  return(0);

 error:
  if(op->obs)
    free((void *) op->obs);
#ifdef HJD
  if(op->epos)
    free((void *) op->epos);
#endif

  return(-1);
}

void output_free (struct lc_output *op) {
  if(op->obs) {
    free((void *) op->obs);
    op->obs = (struct observer *) NULL;
  }
#ifdef HJD
  if(op->epos) {
    free((void *) op->epos);
    op->epos = (double *) NULL;
  }
#endif
}

void output_prepare (struct lc_output *op,
                     struct lc_mef *mefinfo,
                     struct lc_point *lcbuf,
                     struct jpleph_table *jtab,
                     long star, long soff, long stride,
                     double *bjdbuf,
#ifdef HJD
                     double *hjdbuf,
#endif
                     float *fluxbuf, float *fluxerrbuf,
                     double *xlcbuf, double *ylcbuf,
                     float *airbuf, float *habuf,
                     unsigned char *wtbuf,
                     float *locskybuf, float *peakbuf,
                     unsigned char *flagbuf,
                     long *satflag, float *chisq, long *nchisq) {
  long pt;

  int ap;
  long saoff;

  unsigned char flags;

  double s[3], pr, seq[3];

  double airmass;
  float var, fitvar, scvar, tmp;
  
  for(pt = 0; pt < mefinfo->nf; pt++) {
    /* Calculate airmass, HA, BJD */
    if(jtab) {
      /* Apply space motion */
      source_place(op->obs+pt,
                   &(mefinfo->stars[star].src),
                   jtab, TR_MOTION, s, NULL, &pr);
      
      /* Modified BJD(TDB) */
      bjdbuf[soff+pt] = op->obs[pt].tt + (op->obs[pt].dtdb +
                                          bary_delay(op->obs+pt, s, pr)) / DAY;
    }
    else
      bjdbuf[soff+pt] = -999.0;  /* can't compute */
    
    /* Compute airmass */
    if(mefinfo->frames[pt].doairm && jtab) {
      /* Observed place - parallax left out, we don't have it anyway */
      observer_ast2obs(op->obs+pt, s, NULL, pr, TR_TO_OBS_AZ);
      
      /* Airmass */
      airmass = v_airmass(s);
      airbuf[soff+pt] = airmass;
      
      /* Scintillation variance */
      scvar = mefinfo->frames[pt].scvarconst * airmass*airmass*airmass;
      
      /* Hour angle */
      mt_x_v(op->obs[pt].phm, s, seq);
      habuf[soff+pt] = -atan2(seq[1], seq[0]);
    }
    else {
      airbuf[soff+pt] = -999.0;  /* flag unusability */
      habuf[soff+pt] = -999.0;
      
      scvar = 0;
    }
    
#ifdef HJD
    /* Calculate HJD (as UTC) as it was originally done */
    hjdbuf[soff+pt] = mefinfo->mjdref + mefinfo->frames[pt].mjd +
                      hjdcorr(op->epos + 3*pt,
                              mefinfo->stars[star].ra,
                              mefinfo->stars[star].dec);
#endif
    
    /* Chosen aperture and flags */
    flags = 0;
    
    if(lcbuf[pt].aper[mefinfo->stars[star].iap].flux != 0.0) {
      fluxbuf[soff+pt] = mefinfo->zp
                       - lcbuf[pt].aper[mefinfo->stars[star].iap].flux;
      
      if(abs(lcbuf[pt].aper[mefinfo->stars[star].iap].flux > 20)) {
        printf("Warning: daft-looking flux for star %ld point %ld: %.2g\n",
               star+1, pt+1, lcbuf[pt].aper[mefinfo->stars[star].iap].flux);
      }
      
      /* Unset the all saturated flag if not saturated */
      if(lcbuf[pt].conf)
        flags |= FLAG_CONF;
      if(lcbuf[pt].satur) {
        flags |= FLAG_SATUR;
        (*satflag)++;
      }
      
      if(lcbuf[pt].aper[mefinfo->stars[star].iap].fluxvar > 0.0) {
        fitvar = systematic_var_star_frame(lcbuf, mefinfo,
                                           pt, star,
                                           mefinfo->stars[star].iap);

        var = lcbuf[pt].aper[mefinfo->stars[star].iap].fluxvar
            + fitvar
            + scvar;
        fluxerrbuf[soff+pt] = sqrtf(var);
        
        if(mefinfo->stars[star].med > 0.0) {
          tmp = lcbuf[pt].aper[mefinfo->stars[star].iap].flux
              - mefinfo->stars[star].med;
          
          (*chisq) += tmp*tmp / var;
          (*nchisq)++;
        }
      }
      else
        fluxerrbuf[soff+pt] = -999.0;
    }
    else {
      fluxbuf[soff+pt] = -999.0;
      fluxerrbuf[soff+pt] = -999.0;
      flags |= FLAG_NODP;
    }
    
#ifndef NOXY
    xlcbuf[soff+pt] = lcbuf[pt].x;
    ylcbuf[soff+pt] = lcbuf[pt].y;
#else
    xlcbuf[soff+pt] = -999.0;
    ylcbuf[soff+pt] = -999.0;
#endif
    wtbuf[soff+pt] = (lcbuf[pt].comp >> mefinfo->stars[star].iap) & 0x01;
    locskybuf[soff+pt] = lcbuf[pt].sky;
    peakbuf[soff+pt] = lcbuf[pt].peak;
    
    flagbuf[soff+pt] = flags;
    
    /* Do other apertures */
    if(!mefinfo->aperture && mefinfo->apselmode & APSEL_ALL)
      for(ap = op->ap1; ap < op->ap2; ap++) {
        saoff = soff + (ap-op->ap1+1)*stride;
        
        /* Fill in buffer */
        if(lcbuf[pt].aper[ap].flux != 0.0) {
          fluxbuf[saoff+pt] = mefinfo->zp - lcbuf[pt].aper[ap].flux;
          if(lcbuf[pt].aper[ap].fluxvar > 0.0) {
            fitvar = systematic_var_star_frame(lcbuf, mefinfo,
                                               pt, star, ap);

            fluxerrbuf[saoff+pt] = sqrtf(lcbuf[pt].aper[ap].fluxvar +
                                         fitvar + scvar);
          }
          else
            fluxerrbuf[saoff+pt] = -999.0;
        }
        else {
          fluxbuf[saoff+pt] = -999.0;
          fluxerrbuf[saoff+pt] = -999.0;
        }
        
        wtbuf[saoff+pt] = (lcbuf[pt].comp >> ap) & 0x01;
      }
  }
}
