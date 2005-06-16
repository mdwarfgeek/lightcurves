#include <stdio.h>
#include <stdlib.h>

#include "hjd.h"

#include "sla.h"
#include "floatmath.h"
#include "util.h"

/* Conversion factor AU to light seconds */
#define AUSEC 499.0047837

void getearth (double mjd, double ep[3]) {
  double dtt, tt, bev[3], bep[3], ev[3];

  /* Convert MJD to TT */
  dtt = slaDtt(mjd);
  tt = mjd + dtt;

  /* Compute Earth position and velocity using TT as an approximation
   * to the correct TDB.
   */
  slaEvp(tt, 2000.0, bev, bep, ev, ep);
}

double hjdcorr (double ep[3], double ra, double dec) {
  double sv[3], cosdec, pls;

  /* Calculate direction cosines of star */
  cosdec = cos(dec);
  sv[0] = cos(ra) * cosdec;
  sv[1] = sin(ra) * cosdec;
  sv[2] = sin(dec);

  /* Take scalar product with Earth position */
  pls = AUSEC * (ep[0] * sv[0] + ep[1] * sv[1] + ep[2] * sv[2]);

  /* Convert to days and return */
  return(pls / 86400);
}

