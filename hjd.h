#ifndef HJD_H
#define HJD_H

/* Get Earth position */
void getearth (double mjd, double ep[3]);

/* Get correction from MJD to HJD */
double hjdcorr (double ep[3], double ra, double dec);

#endif  /* HJD_H */
