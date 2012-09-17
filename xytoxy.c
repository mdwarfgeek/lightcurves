#include <stdio.h>
#include <stdlib.h>

#include "lightcurves.h"

#include "floatmath.h"
#include "util.h"

#define NITER 5
#define NSIG  3

static void cplate (struct lc_point *combuf, struct lc_star *refbuf, double *wtbuf,
		    long npt, double tr[6]) {
  double xmatx[50][50], vect[50];

  double sx1sq = 0.0, sx1 = 0.0, sy1sq = 0.0, sy1 = 0.0;
  double sx1y1 = 0.0, sx1x2 = 0.0, sy1x2 = 0.0;
  double sx2 = 0.0, sy1y2 = 0.0, sx1y2 = 0.0, sy2 = 0.0;
  double xbar = 0.0, ybar = 0.0, xref = 0.0, yref = 0.0, sw = 0.0;
  double wt, xx1, yy1, xx2, yy2;

  long i;

  /* find <x> and <y> for comparison and reference files */
  for(i = 0; i < npt; i++) {
    wt = (wtbuf ? wtbuf[i] : 1.0);

    xbar += combuf[i].x * wt;
    ybar += combuf[i].y * wt;
    xref += refbuf[i].x * wt;
    yref += refbuf[i].y * wt;
    sw += wt;
  }

  xbar /= sw;
  ybar /= sw;
  xref /= sw;
  yref /= sw;

  /* accumulate sums */
  sw = 0;

  for(i = 0; i < npt; i++) {
    wt = (wtbuf ? wtbuf[i] : 1.0);

    xx1 = combuf[i].x-xbar;           /* zero-mean for more stable solution */
    xx2 = refbuf[i].x-xref;
    yy1 = combuf[i].y-ybar;
    yy2 = refbuf[i].y-yref;
    sx1sq += xx1*xx1 * wt;
    sx1 += xx1 * wt;
    sy1sq += yy1*yy1 * wt;
    sy1 += yy1 * wt;
    sx1y1 += xx1*yy1 * wt;
    sx1x2 += xx1*xx2 * wt;
    sy1x2 += yy1*xx2 * wt;
    sx2 += xx2 * wt;
    sy1y2 += yy1*yy2 * wt;
    sx1y2 += xx1*yy2 * wt;
    sy2 += yy2 * wt;
    sw += wt;
  }

  /* fill in matrix-vector for a b c */
  xmatx[0][0] = sx1sq;
  xmatx[1][0] = sx1y1;
  xmatx[2][0] = sx1;
  xmatx[0][1] = sx1y1;
  xmatx[1][1] = sy1sq;
  xmatx[2][1] = sy1;
  xmatx[0][2] = sx1;
  xmatx[1][2] = sy1;
  xmatx[2][2] = sw;
  vect[0] = sx1x2;
  vect[1] = sy1x2;
  vect[2] = sx2;

  dsolve(xmatx, vect, 3);

  tr[0] = vect[0];
  tr[1] = vect[1];
  tr[2] = vect[2] - xbar*tr[0] - ybar*tr[1] + xref;    /* put back in original system */

  if(vect[0] == 0 && vect[1] == 0 && vect[2] == 0) {
    /* Failed */
    tr[0] = 1.0;
    tr[1] = 0.0;
    tr[2] = xbar - xref;
  }

  /* fill in matrix-vector for d e f */
  xmatx[0][0] = sy1sq;
  xmatx[1][0] = sx1y1;
  xmatx[2][0] = sy1;
  xmatx[0][1] = sx1y1;
  xmatx[1][1] = sx1sq;
  xmatx[2][1] = sx1;
  xmatx[0][2] = sy1;
  xmatx[1][2] = sx1;
  xmatx[2][2] = sw;
  vect[0] = sy1y2;
  vect[1] = sx1y2;
  vect[2] = sy2;
  
  dsolve(xmatx, vect, 3);

  tr[3] = vect[0];
  tr[4] = vect[1];
  tr[5] = vect[2] - xbar*tr[4] - ybar*tr[3] + yref;

  if(vect[0] == 0 && vect[1] == 0 && vect[2] == 0) {
    /* Failed */
    tr[3] = 1.0;
    tr[4] = 0.0;
    tr[5] = ybar - yref;
  }
}

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
    if(data[star].flux > 0.0 &&          /* Has a flux measurement */
       data[star].fluxerr > 0.0 &&       /* And a reliable error */
       mefinfo->stars[star].sigflux[0] > 0 &&
       data[star].flux < mefinfo->sysulim &&  /* Not saturated */
       mefinfo->stars[star].cls == -1) { /* Is classified as stellar */
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
	  dx = tr[0]*data[star].x+tr[1]*data[star].y+tr[2] - mefinfo->stars[star].x;
	  dy = tr[3]*data[star].y+tr[4]*data[star].x+tr[5] - mefinfo->stars[star].y;
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

    if(nuse > 6)
      cplate(data, mefinfo->stars, wtbuf, mefinfo->nstars, tr);
    else
      break;

    nuse = 0;
    nrej = 0;
    for(star = 0; star < mefinfo->nstars; star++) {
      if(rejbuf[star])
	nrej++;
      else {
	dx = tr[0]*data[star].x+tr[1]*data[star].y+tr[2] - mefinfo->stars[star].x;
	dy = tr[3]*data[star].y+tr[4]*data[star].x+tr[5] - mefinfo->stars[star].y;
	
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
