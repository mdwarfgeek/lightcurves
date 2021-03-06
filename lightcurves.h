#ifndef LIGHTCURVES_H
#define LIGHTCURVES_H

#include <fitsio.h>

#ifdef _WIN32
/* CFITSIO's TBYTE define conflicts with a typedef on Win32.
   We don't use it, so get rid of it here. */
#undef TBYTE
#include <windows.h>
#endif

#include "lfa.h"

#ifndef NFLUX
/* Default number of flux measures */
#define NFLUX  3
#endif

#ifndef REFAP
/* Default flux measure to use for normalization */
#define REFAP  0
#endif

/* Flags column in lightcurve file */
#define FLAG_NODP  0x01  /* No data point (off chip) */
#define FLAG_CONF  0x02  /* Aperture contains bad pixels */
#define FLAG_SATUR 0x04  /* Saturated */

/* Default to use all objects within 4 mags below saturation for fit */
#define USEMAG  4

/* Number of coefficients in polynomial */
#define POLY_NCOEFF(degree) ((((degree) + 1) * ((degree) + 2)) / 2)

/* Intrapixel sensitivity map */
struct intra {
  float *map;
  long nbin;
  float binsize;
};

/* "Instrument version" information */
struct instvers {
  int iver;
  long date;
};

struct lc_point {
#ifndef NOXY
  double x;
  double y;
#endif

  struct {
    float flux;
    float fluxvar;
  } aper[NFLUX];

  float sky;   /* local sky */
  float peak;  /* peak counts including sky */

  /* Flags */
  unsigned int comp : 16;
  unsigned int satur : 1;
  unsigned int conf : 1;
};

struct lc_star_segment {
  float corr[NFLUX];
  float medx;
  float sigx;
  float medy;
  float sigy;
};

struct lc_star {
  /* Reference frame info */
  long ptr;
  double x;
  double y;
  double ra;
  double dec;
  float pmra;
  float pmdec;
  int havepm;
  int cls;
  int bflag;
  int cflag;

  struct lc_point ref;

  /* LFA source structure for this star */
  struct source src;

  /* Reference magnitude to which to tie the photometry */
  float refmag;  /* in aperture REFAP */

  /* Information for each piece of the light curve */
  struct lc_star_segment *segs;

  /* Global median and scatter */
  float medflux[NFLUX];
  float sigflux[NFLUX];

  /* Usable as comparison star? */
  int compok;

  /* Used in practice? */
  int used;

  /* Which aperture did we choose? */
  int iap;
  float apradius;

  /* Median flux and RMS for chosen aperture (iap) */
  float med;
  float rms;
  float chisq;
  long nchisq;
};

/* Structure holding systematics fit results */
struct systematic_fit {
  float xbar;
  float ybar;
  double *coeff;
  double *cov;
  int degree;
  int ncoeff;

  float medoff;
  float sigoff;
  float sigm;
  long npt;
};

struct lc_frame {
  /* Frame MJD and exposure time */
  double mjd;
  float exptime;

  /* Observer location and refraction constants for computing
     airmass, ha, bjd, hjd at output stage. */
  double latitude;
  double longitude;
  double height;
  double refco[NREFCO];

  unsigned char doairm;

  /* Scintillation variance constant */
  float scvarconst;

  /* Magnitude zeropoint difference from reference (for absolute cal.) */
  float zpdiff;

  /* Seeing, ellipticity, sky */
  float seeing;
  float ellipt;
  float skylev;
  float skynoise;

  /* Field angle */
  float fang;
  long iang;

  /* Systematic fit results */
  struct systematic_fit sys[NFLUX];

  /* Frame offset and RMS after correction */
  float offset;
  float rms;

  /* Delta(mag) from median - ie. zero order coefficient from fit */
  float extinc;

  /* Uncertainty contribution from systematics fit */
  float sigm;

#ifndef NOXY
  /* Standard 6-coefficient transformation to reference system */
  double tr[6];

  /* Median x,y offsets and uncertainties */
  float xoff;
  float xsig;
  float yoff;
  float ysig;
#endif

  /* MEarth-specific: frame grouping information */
  long split_iexp;
  long split_nexp;

  /* MEarth-specific: weather information */
  float tamb;
  float humid;
  float press;
  float skytemp;

  /* MEarth-specific: realtime status information */
  long rtstat;

  /* MEarth-specific: scheduler information */
  int schpri;
  float schcad;
  char schtype[FLEN_VALUE];

  /* MEarth-specific: "instrument version" */
  struct instvers *instvers;

  /* Kepler-specific: cadence number */
  long cadencenum;

  /* Segment number */
  int iseg;
};

struct lc_segment {
  struct instvers *instvers;
  int iang;

  /* Median and rms x,y position offsets for this MEF, segment */
  float medxoff;
  float sigxoff;
  float medyoff;
  float sigyoff;
};

struct lc_mef {
  /* Star info */
  struct lc_star *stars;
  long nstars;  /* number of stars from table */
  long nrows;  /* original number of rows in table */

  /* Warnings we have emitted */
  unsigned char warned;
#define WARNED_SATLEV   0x01
#define WARNED_EXPTIME  0x02
#define WARNED_RCORE    0x04
#define WARNED_GAIN     0x08
#define WARNED_MAGZPT   0x10
#define WARNED_TIME     0x20

  /* Use theoretical sky level calc? */
  int theosky;

  /* Degree of polynomial fit to be applied, zero for none */
  int degree;
  int ncoeff;

  /* Overall reference frame info */
  float zp;
  float satmag;
  float reffang;
  float refexp;
  float refextinct;
  float refairmass;
  float refmagzpt;
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

  /* Systematics fit buffers */
  double *sysbuf_coeff;
  double *sysbuf_cov;
  int nsysalloc;

  /* Systematics fitting upper mag limit */
  float sysulim;
  float sysllim;

  /* Use this aperture in output (0 to select automatically) */
  long aperture;

  /* Aperture selection output modes */
  long apselmode;
#define APSEL_SEL 0x01
#define APSEL_ALL 0x02  /* not supported on own yet, todo */

  /* Account for meridian offsets? */
  long domerid;

  /* Segment info */
  struct lc_segment *segs;
  long nseg;
};

struct lc_output {
  struct observer *obs;
#ifdef HJD
  double *epos;
#endif

  /* Apertures we processed */
  int ap1;
  int ap2;

  int napcol;  /* number of columns for medflux etc. */
  int nlapcol;  /* number of columns for flux etc. */
};

/* Structure holding disk buffer information */
struct buffer_info {
  /* File info */
#ifdef _WIN32
  char filename[MAX_PATH];
  HANDLE hf;
  HANDLE hm;
#else
  char filename[1024];
#endif
  int fd;

  /* Buffer sizing */
  long nobj;
  long nmeas;

  /* Memory-mapped array */
  unsigned char *buf;
};

struct input_file {
  char *filename;
  int ient;  /* index within container (-1 = normal file) */
  int iarg;  /* argument this input file came from */
  char *arg;  /* pointer to argument text */
};

/* Globals */
extern float flux_apers[];
extern int verbose;

/* Disk buffer management: buffer.c */
int buffer_init (struct buffer_info *b, char *errstr);
int buffer_alloc (struct buffer_info *b, long nobj, long nmeas, char *errstr);
void buffer_close (struct buffer_info *b);

int buffer_fetch_frame (struct buffer_info *b, struct lc_point *buf,
			long noff, long nelem,
			long iframe, char *errstr);
int buffer_fetch_object (struct buffer_info *b, struct lc_point *buf,
			 long noff, long nelem,
			 long ipoint, char *errstr);

int buffer_put_frame (struct buffer_info *b, struct lc_point *buf,
		      long noff, long nelem,
		      long iframe, char *errstr);
int buffer_put_object (struct buffer_info *b, struct lc_point *buf,
		       long noff, long nelem,
		       long ipoint, char *errstr);

/* Aperture combination: chooseap.c */
int chooseap (struct buffer_info *buf, struct lc_mef *mefinfo,
	      struct lc_point *ptbuf, float *medbuf, int doapcor, char *errstr);

/* "Instrument version": instvers.c */
int read_instvers (char *filename, struct instvers **instverslist_r, int *ninstvers_r,
		   char *errstr);

/* Intrapixel sensitivity correction: intra.c */
int read_intra (char *filename, struct intra *intralist, int nmefs,
		char *errstr);

float calc_intra (float x, float y, struct intra *corr);

/* Main routine: lightcurves.c */
int lightcurves (struct buffer_info *buf, struct lc_mef *mefinfo,
		 int norenorm, int noastrom, int niter, char *errstr);
int lightcurves_append (struct buffer_info *buf, struct lc_mef *mefinfo,
			int noastrom, char *errstr);

/* Common output processing: output.c */
int output_init (struct lc_output *op,
                 struct lc_mef *mefinfo,
                 struct dtai_table *dtab,
                 struct iers_table *itab,
                 struct jpleph_table *jtab,
                 struct jpleph_table *ttab,
                 char *errstr);
void output_free (struct lc_output *op);
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
                     long *satflag, float *chisq, long *nchisq);

#ifdef PLOTS
/* Diagnostic plots: plots.c */
int do_plots (struct lc_mef *meflist, int nmefs,
	      float medsat, float medlim, float umlim, float lmlim,
              int *outcls, int wantoutcls, char *errstr);
#endif

int plot_corr (float *beforehist, float *beforewthist,
	       float *corrhist, float *corrwthist,
	       float *afterhist, float *afterwthist,
	       long nbinx, long nbiny, float binsize,
	       float xmin, float xmax, float ymin, float ymax,
	       float medoff, float sigoff, char *errstr);

/* Systematics correction: systematic.c */
void systematic_select (struct lc_star *stars, long nstars, float fmin, float fmax);
int systematic_fit (struct lc_point *data, struct lc_mef *mefinfo, long frame, long meas,
		    float *medbuf, int degree, struct systematic_fit *f, char *errstr);

int systematic_apply_frame (struct lc_point *data, struct lc_mef *mefinfo,
                            long frame, long meas, char *errstr);
int systematic_apply_star (struct lc_point *data, struct lc_mef *mefinfo,
                           long star, long meas, char *errstr);

float systematic_var_star_frame (struct lc_point *data,
                                 struct lc_mef *mefinfo,
                                 long frame, long star,
                                 long meas);

#ifndef NOXY
/* Plate solutions: xytoxy.c */
int xytoxy (struct lc_point *data, struct lc_mef *mefinfo,
	    double tr[6], char *errstr);
#endif

/* Catalogue and list driven file reading: readfits.c */
int read_lc (fitsfile *fits, struct lc_mef *mefinfo,
	     char *errstr);
int read_ref (fitsfile *fits, struct lc_mef *mefinfo,
	      int diffmode, float satlev,
	      float sysllim, float sysulim,
	      int *outcls, int wantoutcls,
	      char *errstr);
int read_cat (struct input_file *catfile, int iframe,
              int mef, struct lc_mef *mefinfo,
	      struct buffer_info *buf,
	      int dointra, struct intra *icorr,
	      int doinstvers, struct instvers *instverslist, int ninstvers,
	      int diffmode, float satlev,
	      char *errstr);

/* Utility functions: dsolve.c, dmatinv.c, filelist.c, lightcurves.c */
void dsolve (double *a, double *b, int m);
void dmatinv (double *a, int m);

struct input_file *read_file_list (int argc, char **argv,
                                   int *nf_r, char *errstr);

int checkoutcls (int cls, int *outcls, int wantoutcls);

#endif  /* LIGHTCURVES_H */
