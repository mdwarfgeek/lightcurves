#include <stdlib.h>
#include <string.h>

#include "sla.h"

#ifndef CSLALIB

#ifdef GNUFortran
#include <g2c.h>
#else
typedef unsigned int ftnlen;
#endif

/* Prototypes */
void sla_addet_ (double *, double *, double *, double *, double *);
double sla_airmas_ (double *);
void sla_aoppa_ (double *, double *, double *, double *, double *,
		 double *, double *, double *, double *, double *,
		 double *, double *, double *);
void sla_aoppat_ (double *, double *);
void sla_aopqk_ (double *, double *, double *,
		 double *, double *, double *, double *, double *);
void sla_cldj_ (int *, int *, int *, double *, int *);
void sla_cr2af_ (int *, float *, char *, int *, ftnlen);
void sla_cr2tf_ (int *, float *, char *, int *, ftnlen);
void sla_ctf2d_ (int *, int *, float *, float *, int *);
double sla_dtt_ (double *);
void sla_e2h_ (float *, float *, float *, float *, float *);
void sla_evp_ (double *, double *, double *, double *, double *, double *);
double sla_epb_ (double *);
double sla_epj_ (double *);
void sla_fk45z_ (double *, double *, double *, double *, double *);
void sla_map_ (double *, double *, double *, double *, double *,
	       double *, double *, double *, double *, double *);
void sla_mappa_ (double *, double *, double *);
void sla_mapqkz_ (double *, double *, double *, double *, double *);
void sla_pm_ (double *, double *, double *, double *, double *,
	      double *, double *, double *, double *, double *);
void sla_preces_ (char *, double *, double *, double *, double *, ftnlen);
float sla_range_ (float *);
float sla_ranorm_ (float *);
void sla_rdplan_ (double *, int *, double *, double *,  double *, double *, double *);
void sla_subet_ (double *, double *, double *, double *, double *);

/* Wrappers */

void slaAddet (double rm, double dm, double eq, double *rc, double *dc) {
  sla_addet_ (&rm, &dm, &eq, rc, dc);
}

double slaAirmas (double zd) {
  return sla_airmas_(&zd);
}

void slaAoppa (double date, double dut, double elongm, double phim, double hm,
	       double xp, double yp, double tdk, double pmb, double rh,
	       double wl, double tlr, double aoprms[14]) {
  sla_aoppa_(&date, &dut, &elongm, &phim, &hm,
	     &xp, &yp, &tdk, &pmb, &rh,
	     &wl, &tlr, aoprms);
}

void slaAoppat (double date, double aoprms[14]) {
  sla_aoppat_(&date, aoprms);
}

void slaAopqk (double rap, double dap, double aoprms[14],
	       double *aob, double *zob, double *hob, double *dob, double *rob) {
  sla_aopqk_(&rap, &dap, aoprms, aob, zob, hob, dob, rob);
}

void slaCldj (int iy, int im, int id, double *djm, int *j) {
  sla_cldj_(&iy, &im, &id, djm, j);
}

void slaCr2af (int ndp, float angle, char *sign, int idmsf[4]) {
  sla_cr2af_(&ndp, &angle, sign, idmsf, 1);
}

void slaCr2tf (int ndp, float angle, char *sign, int ihmsf[4]) {
  sla_cr2tf_(&ndp, &angle, sign, ihmsf, 1);
}

void slaCtf2d (int ihour, int imin, float sec, float *days, int *j) {
  sla_ctf2d_(&ihour, &imin, &sec, days, j);
}

double slaDtt (double dju) {
  return sla_dtt_ (&dju);
}

void slaE2h (float ha, float dec, float phi, float *az, float *el) {
  sla_e2h_(&ha, &dec, &phi, az, el);
}

void slaEvp (double date, double deqx,
	     double dvb[3], double dpb[3], double dvh[3], double dph[3]) {
  sla_evp_(&date, &deqx, dvb, dpb, dvh, dph);
}

double slaEpb (double date) {
  return sla_epb_ (&date);
}

double slaEpj (double date) {
  return sla_epj_ (&date);
}

void slaFk45z (double r1950, double d1950, double bepoch, double *r2000, double *d2000) {
  sla_fk45z_ (&r1950, &d1950, &bepoch, r2000, d2000);
}

void slaMap (double rm, double dm, double pr, double pd, double px, double rv,
	     double eq, double date, double *ra, double *da) {
  sla_map_ (&rm, &dm, &pr, &pd, &px, &rv, &eq, &date, ra, da);
}

void slaMappa (double eq, double date, double amprms[21]) {
  sla_mappa_(&eq, &date, amprms);
}

void slaMapqkz (double rm, double dm, double amprms[21], double *ra, double *da) {
  sla_mapqkz_(&rm, &dm, amprms, ra, da);
}

void slaPm (double r0, double d0, double pr, double pd, double px, double rv,
	    double ep0, double ep1, double *r1, double *d1) {
  sla_pm_ (&r0, &d0, &pr, &pd, &px, &rv, &ep0, &ep1, r1, d1);
}

void slaPreces (char *system, double ep0, double ep1, double *ra, double *dc) {
  sla_preces_ (system, &ep0, &ep1, ra, dc, strlen(system));
}

float slaRange (float a) {
  return sla_range_ (&a);
}

float slaRanorm (float a) {
  return sla_ranorm_ (&a);
}

void slaRdplan (double date, int np, double elong, double phi, double *ra, double *dec,
		double *diam) {
  sla_rdplan_ (&date, &np, &elong, &phi, ra, dec, diam);
}

void slaSubet (double rc, double dc, double eq, double *rm, double *dm) {
  sla_subet_ (&rc, &dc, &eq, rm, dm);
}

#endif
