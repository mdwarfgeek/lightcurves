#ifndef __LIGHTCURVES_H__
#define __LIGHTCURVES_H__

/* Number of flux measures */
#define NFLUX  3

struct lc_point {
  float flux;
  float fluxerr;
  unsigned char satur;
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

  /* Median flux and sigma in each aperture */
  float medflux[NFLUX];
  float sigflux[NFLUX];

  /* Which aperture did we use? */
  float apradius;

  /* RMS and chisq */
  float rms;
  float chisq;
  long nchisq;
};

struct lc_frame {
  /* Frame MJD */
  double mjd;

  /* Frame offset and RMS after correction */
  float offset;
  float rms;
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
  float refexp;
  float refsigma;
  float refflim;

  float refgain;  /* for theoretical curve */
  float refrcore;

  float apcor[NFLUX];
  float percorr;

  /* Saturation clip level used for systematics fit */
  float satclip[NFLUX];

  /* Average of sigma for theoretical curve */
  float avsigma;

  /* Reference MJD */
  double mjdref;

  /* Individual frame info */
  struct lc_frame *frames;
  long nf;
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
		 int noapsel, int dopca, char *errstr);

/* PCA-like systematics removal: pcasys.c */
int pcasys (struct buffer_info *buf, struct lc_point *ptbuf, struct lc_mef *mefinfo,
	    long meas, float *medbuf, char *errstr);

/* Diagnostic plots: plots.c */
int do_plots (struct lc_mef *meflist, int nmefs,
	      float medsat, float medlim, float sysbodge, char *errstr);

int plot_corr (float *beforehist, float *beforewthist,
	       float *corrhist, float *corrwthist,
	       float *afterhist, float *afterwthist,
	       long nbinx, long nbiny, float binsize,
	       float xmin, float xmax, float ymin, float ymax,
	       float medoff, float sigoff, char *errstr);

/* Systematics correction: systematic.c */
int systematic_fit (struct lc_point *data, struct lc_mef *mefinfo, long frame, long meas,
		    float *medbuf, struct systematic_fit *f,
		    float *med_r, float *rms_r, char *errstr);
int systematic_apply (struct lc_point *data, struct lc_mef *mefinfo, long frame, long meas,
		      float *medbuf, struct systematic_fit *f, char *errstr);

/* Utility functions: dsolve.c, linear.c, medsig.c, sortfloat.c */
void dsolve (double a[50][50], double b[50], int m);
int hanning (float xbuf[], int npt, char *errstr);
void medsig (float *a, long n, float *median_r, float *sigma_r);
void sortfloat (float *ia, long n);

#endif  /* __LIGHTCURVES_H__ */
