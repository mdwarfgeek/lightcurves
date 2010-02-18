#ifndef __SLA_H__
#define __SLA_H__

#ifdef CSLALIB
#include <slalib.h>
#include <slamac.h>
#else
void slaAddet (double rm, double dm, double eq, double *rc, double *dc);
double slaAirmas (double zd);
void slaAoppa (double date, double dut, double elongm, double phim, double hm,
	       double xp, double yp, double tdk, double pmb, double rh,
	       double wl, double tlr, double aoprms[14]);
void slaAoppat (double date, double aoprms[14]);
void slaAopqk (double rap, double dap, double aoprms[14],
	       double *aob, double *zob, double *hob, double *dob, double *rob);
void slaCldj (int iy, int im, int id, double *djm, int *j);
void slaCr2af (int ndp, float angle, char *sign, int idmsf[4]);
void slaCr2tf (int ndp, float angle, char *sign, int ihmsf[4]);
void slaCtf2d (int ihour, int imin, float sec, float *days, int *j);
double slaDtt (double dju);
void slaE2h (float ha, float dec, float phi, float *az, float *el);
void slaEvp (double date, double deqx,
	     double dvb[3], double dpb[3], double dvh[3], double dph[3]);
double slaEpb (double date);
double slaEpj (double date);
void slaFk45z (double r1950, double d1950, double bepoch, double *r2000, double *d2000);
void slaMap (double rm, double dm, double pr, double pd, double px, double rv,
	     double eq, double date, double *ra, double *da);
void slaMappa (double eq, double date, double amprms[21]);
void slaMapqkz (double rm, double dm, double amprms[21], double *ra, double *da);
void slaPm (double r0, double d0, double pr, double pd, double px, double rv,
	    double ep0, double ep1, double *r1, double *d1);
void slaPreces (char *system, double ep0, double ep1, double *ra, double *dc);
float slaRange (float a);
float slaRanorm (float a);
void slaRdplan (double date, int np, double elong, double phi, double *ra, double *dec,
		double *diam);
void slaSubet (double rc, double dc, double eq, double *rm, double *dm);
#endif

#endif  /* __SLA_H__ */
