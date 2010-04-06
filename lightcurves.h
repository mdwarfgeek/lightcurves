#ifndef __LIGHTCURVES_H__
#define __LIGHTCURVES_H__

#include <fitsio.h>

/* Number of flux measures */
#define NFLUX  3

/* Flags column in lightcurve file */
#define FLAG_NODP  0x01  /* No data point (off chip) */
#define FLAG_CONF  0x02  /* Aperture contains bad pixels */
#define FLAG_SATUR 0x04  /* Saturated */

/* Default to use all objects within 4 mags below saturation for fit */
#define USEMAG  4

struct lc_point {
  float x;
  float y;
  float flux;
  float fluxerr;
  float fluxerrcom;  /* combined uncertainty including fit */
  float airmass;
  float ha;
  float wt;  /* weight given in computing polynomial corr. */
  unsigned char satur : 1;
  unsigned char conf : 1;
};

struct lc_star {
  /* Reference frame info */
  long ptr;
  float x;
  float y;
  float ra;
  float dec;
  int cls;
  int bflag;
  long cflag;

  struct lc_point ref[NFLUX];

  /* Reference magnitude to which to tie the photometry */
  float refmag;

  /* Median flux and sigma in each aperture */
  float medflux[NFLUX];
  float sigflux[NFLUX];

  /* Meridian flip correction if enabled */
  float merid[NFLUX];

  /* Which aperture did we use? */
  float apradius;

  /* Median flux, RMS and chisq for chosen aperture */
  float med;
  float rms;
  float chisq;
  long nchisq;
};

struct lc_frame {
  /* Frame MJD */
  double mjd;

  /* Seeing */
  float seeing;

  /* Field angle */
  float fang;
  long iang;

  /* Frame offset and RMS after correction */
  float offset;
  float rms;

  /* Delta(mag) from median - ie. zero order coefficient from fit */
  float extinc;

  /* Uncertainty contribution from systematics fit */
  float sigm;

  /* MEarth-specific: frame grouping information */
  long split_iexp;
  long split_nexp;

  /* MEarth-specific: weather information */
  float tamb;
  float humid;
  float press;
  float skytemp;
};

struct lc_mef {
  /* Star info */
  struct lc_star *stars;
  long nstars;

  /* Degree of polynomial fit to be applied, zero for none */
  int degree;

  /* Overall reference frame info */
  float zp;
  float satmag;
  float reffang;
  int havefang;
  float refexp;
  float refsigma;
  float refflim;

  float refgain;  /* for theoretical curve */
  float refrcore;

  float apcor[NFLUX];
  float percorr;

  /* Saturation clip level used for systematics fit */
  float satclip[NFLUX];

  /* Filter */
  char filter[64];

  /* Average of sigma for theoretical curve */
  float avsigma;

  /* Average of sky fit error for theoretical curve */
  float avskyfit;

  /* Average of aperture correction for theoretical curve */
  float avapcor;

  /* Average of extinction for theoretical curve */
  float avextinc;

  /* Average of sigma_m for theoretical curve */
  float avsigm;

  /* Scintillation for theoretical curve */
  float avscint;

  /* Reference MJD */
  double mjdref;

  /* Individual frame info */
  struct lc_frame *frames;
  long nf;

  /* Systematics fitting upper mag limit */
  float sysulim;
  float sysllim;

  /* Use this aperture in output (0 to select automatically) */
  long aperture;

  /* Account for meridian offsets? */
  long domerid;

  /* Number of updates made so far - for setting update count flags */
  long nupdate;
};

/* Structure holding disk buffer information */
struct buffer_info {
  /* File info */
  char filename[1024];
  int fd;

  /* Buffer sizing */
  long nobj;
  long nmeas;
  int nflux;

  /* Memory-mapped array */
  unsigned char *buf;
};

/* Structure holding systematics fit results */
struct systematic_fit {
  float xbar;
  float ybar;
  double coeff[50];
  int degree;
};

/* Intrapixel sensitivity map */
struct intra {
  float *map;
  long nbin;
  float binsize;
};

/* Globals */
extern float flux_apers[NFLUX];
extern int verbose;

/* Disk buffer management: buffer.c */
int buffer_init (struct buffer_info *b, char *errstr);
int buffer_alloc (struct buffer_info *b, long nobj, long nmeas, int nflux, char *errstr);
void buffer_close (struct buffer_info *b);

int buffer_fetch_frame (struct buffer_info *b, struct lc_point *buf,
			long noff, long nelem,
			long iframe, int iflux, char *errstr);
int buffer_fetch_object (struct buffer_info *b, struct lc_point *buf,
			 long noff, long nelem,
			 long ipoint, int iflux, char *errstr);

int buffer_put_frame (struct buffer_info *b, struct lc_point *buf,
		      long noff, long nelem,
		      long iframe, int iflux, char *errstr);
int buffer_put_object (struct buffer_info *b, struct lc_point *buf,
		       long noff, long nelem,
		       long ipoint, int iflux, char *errstr);

/* Aperture combination: chooseap.c */
int chooseap (struct buffer_info *buf, struct lc_mef *mefinfo,
	      struct lc_point *ptbuf, float *medbuf, char *errstr);

/* Intrapixel sensitivity correction: intra.c */
int read_intra (char *filename, struct intra *intralist, int nmefs,
		char *errstr);

float calc_intra (float x, float y, struct intra *corr);

/* Main routine: lightcurves.c */
int lightcurves (struct buffer_info *buf, struct lc_mef *mefinfo,
		 int norenorm, char *errstr);
int lightcurves_append (struct buffer_info *buf, struct lc_mef *mefinfo,
			char *errstr);

/* Diagnostic plots: plots.c */
int do_plots (struct lc_mef *meflist, int nmefs,
	      float medsat, float medlim, float umlim, float lmlim, char *errstr);

int plot_corr (float *beforehist, float *beforewthist,
	       float *corrhist, float *corrwthist,
	       float *afterhist, float *afterwthist,
	       long nbinx, long nbiny, float binsize,
	       float xmin, float xmax, float ymin, float ymax,
	       float medoff, float sigoff, char *errstr);

/* Systematics correction: systematic.c */
int systematic_fit (struct lc_point *data, struct lc_mef *mefinfo, long frame, long meas,
		    float *medbuf, int degree, struct systematic_fit *f,
		    float *med_r, float *rms_r, float *sigm_r, long *npt_r,
		    char *errstr);
int systematic_apply (struct lc_point *data, struct lc_mef *mefinfo, long frame, long meas,
		      float *medbuf, struct systematic_fit *f,
		      float sigm, char *errstr);

/* Catalogue and list driven file reading: readfits.c */
int read_lc (fitsfile *fits, struct lc_mef *mefinfo,
	     char *errstr);
int read_ref (fitsfile *fits, struct lc_mef *mefinfo,
	      int diffmode, float satlev,
	      char *errstr);
int read_cat (char *catfile, int iframe, int mef, struct lc_mef *mefinfo,
	      struct buffer_info *buf,
	      int dointra, struct intra *icorr,
	      int diffmode, float satlev,
	      char *errstr);

/* Utility functions: dsolve.c, linear.c, medsig.c, sortfloat.c */
void dsolve (double a[50][50], double b[50], int m);
int hanning (float xbuf[], int npt, char *errstr);
void medsig (float *a, long n, float *median_r, float *sigma_r);

void sortfloat (float *ia, long n);
void sortlong (long *ia, long n);

#endif  /* __LIGHTCURVES_H__ */
