#ifndef __HJD_H__
#define __HJD_H__

/* Get Earth position */
void getearth (double mjd, double ep[3]);

/* Get correction from MJD to HJD */
double hjdcorr (double ep[3], double ra, double dec);

#endif  /* __HJD_H__ */
