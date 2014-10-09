#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "lightcurves.h"

#include "util.h"

#define NITER 5
#define NSIG  3

int xytoxy (struct lc_point *data, struct lc_mef *mefinfo,
	    double tr[6], char *errstr) {
  double *wtbuf = (double *) NULL;
  float *tmpbufx = (float *) NULL, *tmpbufy;
  unsigned char *rejbuf = (unsigned char *) NULL;

  long star, opt, nuse, nrej;
  double dx, dy, err;
  int iter;

  float averrx, averry, averr;

  /* Init transformation */
  tr[0] = 1.0;
  tr[1] = 0.0;
  tr[2] = 0.0;
  tr[3] = 1.0;
  tr[4] = 0.0;
  tr[5] = 0.0;

  /* Allocate workspace */
  wtbuf = (double *) malloc(mefinfo->nstars * sizeof(double));
  tmpbufx = (float *) malloc(2 * mefinfo->nstars * sizeof(float));
  rejbuf = (unsigned char *) malloc(mefinfo->nstars * sizeof(unsigned char));
  if(!wtbuf || !tmpbufx || !rejbuf) {
    report_syserr(errstr, "malloc");
    goto error;
  }

  tmpbufy = tmpbufx + mefinfo->nstars;

  /* Decide initial weights and rejection */
  opt = 0;

  for(star = 0; star < mefinfo->nstars; star++) {
    if(data[star].aper[0].flux > 0.0 &&          /* Has a flux measurement */
       data[star].aper[0].fluxvar > 0.0 &&       /* And a reliable error */
       mefinfo->stars[star].sigflux[0] > 0 &&
       data[star].aper[0].flux < mefinfo->sysulim &&  /* Not saturated */
       mefinfo->stars[star].cls == -1) {   /* Is classified as stellar */
      wtbuf[star] = 1.0;
      rejbuf[star] = 0;
      opt++;
    }
    else {
      wtbuf[star] = 0.0;
      rejbuf[star] = 1;
    }
  }

  averr = 0;

  for(iter = 0; iter < NITER; iter++) {
    nuse = 0;

    if(iter > 0) {
      for(star = 0; star < mefinfo->nstars; star++)
	if(!rejbuf[star]) {
	  dx = dplate_tr_x(data[star].x, data[star].y, tr)
             - mefinfo->stars[star].x;
          dy = dplate_tr_y(data[star].x, data[star].y, tr)
             - mefinfo->stars[star].y;
	  err = sqrt(dx*dx+dy*dx);
	
	  if(err > NSIG * averr) {
	    wtbuf[star] = 0;
	    rejbuf[star] = 1;
	  }
	  else
	    nuse++;
	}
    }
    else {
      for(star = 0; star < mefinfo->nstars; star++)
	if(!rejbuf[star])
	  nuse++;
    }    

    if(nuse >= 3)
      dplate(data, offsetof(struct lc_point, x), sizeof(struct lc_point),
             data, offsetof(struct lc_point, y), sizeof(struct lc_point),
             mefinfo->stars, offsetof(struct lc_star, x), sizeof(struct lc_star),
             mefinfo->stars, offsetof(struct lc_star, y), sizeof(struct lc_star),
             wtbuf, 0, sizeof(double),
             mefinfo->nstars, 6, tr);
    else
      break;

    nuse = 0;
    nrej = 0;
    for(star = 0; star < mefinfo->nstars; star++) {
      if(rejbuf[star])
	nrej++;
      else {
        dx = dplate_tr_x(data[star].x, data[star].y, tr)
           - mefinfo->stars[star].x;
        dy = dplate_tr_y(data[star].x, data[star].y, tr)
           - mefinfo->stars[star].y;
	
	tmpbufx[nuse] = fabsf(dx);
	tmpbufy[nuse] = fabsf(dy);

	nuse++;
      }
    }
    
    medsig(tmpbufx, nuse, &averrx, (float *) NULL);
    medsig(tmpbufx, nuse, &averry, (float *) NULL);
    
    averrx *= 1.48;  /* converts to rms equivalent */
    averry *= 1.48;

    averr = sqrt(averrx*averrx + averry*averry);

    if(verbose > 2) {
      printf("Iteration %d Av error %.4f (%.4f %.4f) nfit %ld (rej %ld)\n",
	     iter+1, averr, averrx, averry, nuse, nrej);
    }
  }
  
  if(verbose > 2)
    printf("Transform constants: %9.6f %9.6f %9.6f\n"
	   "                     %9.6f %9.6f %9.6f\n",
	   tr[0], tr[1], tr[2],
	   tr[3], tr[4], tr[5]);
  
  free((void *) wtbuf);
  wtbuf = (double *) NULL;
  free((void *) tmpbufx);
  tmpbufx = (float *) NULL;
  free((void *) rejbuf);
  rejbuf = (unsigned char *) NULL;

  return(0);

 error:
  if(wtbuf)
    free((void *) wtbuf);
  if(tmpbufx)
    free((void *) tmpbufx);
  if(rejbuf)
    free((void *) rejbuf);

  return(1);
}
