#include <sys/types.h>

#if defined(HAVE_MMAP) && !defined(_WIN32)
#include <sys/mman.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>

#include "lightcurves.h"

#include "util.h"

static int buffer_read (struct buffer_info *b, struct lc_point *rbuf,
			off_t offset, int nread, char *errstr);
static int buffer_write (struct buffer_info *b, struct lc_point *wbuf,
			 off_t offset, int nwrite, char *errstr);

/* Implements a simple disk backed buffer for storing intermediate
 * stages of lightcurve processing.  We only use this for the fluxes
 * since there will be enough room in RAM for everything else.
 *
 * This uses the simple-minded approach - the data is dumped out
 * as simple arrays of struct lc_point into a single temporary file
 * in the current directory.
 *
 */

int buffer_init (struct buffer_info *b, char *errstr) {
#ifdef _WIN32
  /* Obtain temporary file name */
  if(GetTempFileName(TEXT("."), 
                     TEXT("lightcurves"),
                     0,
                     b->filename) == 0) {
    report_err(errstr, "GetTempFileName failed");
    goto error;
  }

  /* Open it */
  b->hf = CreateFile(b->filename,
                     GENERIC_READ | GENERIC_WRITE,
                     0,
                     NULL,
                     CREATE_ALWAYS,
                     FILE_ATTRIBUTE_TEMPORARY | FILE_FLAG_DELETE_ON_CLOSE,
                     NULL);
  if(b->hf == INVALID_HANDLE_VALUE) {
    report_err(errstr, "CreateFile(%s) failed\n", b->filename);
    goto error;
  }

  b->fd = _open_osfhandle((intptr_t) b->hf, 0);
  if(b->fd == -1) {
    report_syserr(errstr, "open");
    goto error;
  }
#else
  int rv;

  /* Generate temporary file */
  strncpy(b->filename, "lightcurves.XXXXXX", sizeof(b->filename));
  b->filename[sizeof(b->filename)-1] = '\0';

  b->fd = mkstemp(b->filename);
  if(b->fd == -1) {
    report_syserr(errstr, "mkstemp");
    goto error;
  }

  /* Set it to unlink when closed */
  rv = unlink(b->filename);
  if(rv < 0) {
    report_syserr(errstr, "unlink: %s", b->filename);
    goto error;
  }
#endif

  b->buf = (unsigned char *) NULL;

  return(0);

 error:
  return(1);
}

int buffer_alloc (struct buffer_info *b, long nobj, long nmeas, char *errstr) {
  off_t totsize;
  int rv;

#ifdef HAVE_MMAP
  /* Unmap the buffer */
  if(b->buf) {
#ifdef _WIN32
    if(!UnmapViewOfFile((LPCVOID) b->buf)) {
      report_err(errstr, "UnmapViewOfFile failed");
      goto error;
    }

    if(!CloseHandle(b->hm)) {
      report_err(errstr, "CloseHandle failed");
      goto error;
    }
#else
    rv = munmap((void *) b->buf,
		((off_t) b->nobj) * ((off_t) b->nmeas) *
		sizeof(struct lc_point));
    if(rv < 0) {
      report_syserr(errstr, "munmap");
      goto error;
    }
#endif

    b->buf = (unsigned char *) NULL;
  }
#endif

  /* Calculate total required file size */
  totsize = ((off_t) nobj) * ((off_t) nmeas) * sizeof(struct lc_point);

  /* Truncate the file to that size */
  rv = ftruncate(b->fd, totsize);
  if(rv < 0) {
    report_syserr(errstr, "ftruncate");
    goto error;
  }

  if(totsize > 0) {
#ifdef HAVE_MMAP
    /* Map the buffer */
#ifdef _WIN32
    b->hm = CreateFileMapping(b->hf,
                              NULL,
                              PAGE_READWRITE,
                              0,
                              0,
                              NULL);
    if(!b->hm) {
      report_err(errstr, "CreateFileMapping failed");
      goto error;
    }

    b->buf = (unsigned char *) MapViewOfFile(b->hm,
                                             FILE_MAP_ALL_ACCESS,
                                             0,
                                             0,
                                             totsize);
    if(!b->buf) {
      report_err(errstr, "MapViewOfFile failed");
      goto error;
    }
#else
    b->buf = (unsigned char *) mmap((void *) NULL, totsize,
                                    PROT_READ | PROT_WRITE, MAP_SHARED, b->fd, 0);
    if(b->buf == ((unsigned char *) MAP_FAILED)) {
      report_syserr(errstr, "mmap");
      goto error;
    }
#endif
#endif
  }

  b->nobj = nobj;
  b->nmeas = nmeas;

  return(0);

 error:
  return(1);
}

void buffer_close (struct buffer_info *b) {

#ifdef HAVE_MMAP
  /* Unmap the buffer */
  if(b->buf) {
#ifdef _WIN32
    UnmapViewOfFile((LPCVOID) b->buf);
    CloseHandle(b->hm);
#else
    munmap((void *) b->buf,
           ((off_t) b->nobj) * ((off_t) b->nmeas) *
           sizeof(struct lc_point));
#endif

    b->buf = (unsigned char *) NULL;
  }
#endif

  /* Close the file */
  close(b->fd);
}

int buffer_fetch_frame (struct buffer_info *b, struct lc_point *buf,
			long noff, long nelem,
			long iframe, char *errstr) {
  off_t offset;

  offset = ((off_t) iframe) * ((off_t) b->nobj) +
           ((off_t) noff);

  if(buffer_read(b, buf, offset, nelem, errstr))
    goto error;

  return(0);

 error:
  return(1);
}

int buffer_fetch_object (struct buffer_info *b, struct lc_point *buf,
			 long noff, long nelem,
			 long ipoint, char *errstr) {
  off_t offset;
  long iframe;

  for(iframe = 0; iframe < nelem; iframe++) {
    offset = ((off_t) (iframe+noff)) * ((off_t) b->nobj) + ((off_t) ipoint);

    if(buffer_read(b, buf + iframe, offset, 1, errstr))
      goto error;
  }

  return(0);

 error:
  return(1);
}

int buffer_put_frame (struct buffer_info *b, struct lc_point *buf,
		      long noff, long nelem,
		      long iframe, char *errstr) {
  off_t offset;

  offset = ((off_t) iframe) * ((off_t) b->nobj) +
           ((off_t) noff);

  if(buffer_write(b, buf, offset, nelem, errstr))
    goto error;

  return(0);

 error:
  return(1);
}

int buffer_put_object (struct buffer_info *b, struct lc_point *buf,
		       long noff, long nelem,
		       long ipoint, char *errstr) {
  off_t offset;
  long iframe;

  for(iframe = 0; iframe < nelem; iframe++) {
    offset = ((off_t) (iframe+noff)) * ((off_t) b->nobj) + ((off_t) ipoint);

    if(buffer_write(b, buf + iframe, offset, 1, errstr))
      goto error;
  }

  return(0);

 error:
  return(1);
}

#ifdef HAVE_MMAP
static int buffer_read (struct buffer_info *b, struct lc_point *rbuf,
			off_t offset, int nread, char *errstr) {
  /* Convert to bytes */
  offset *= sizeof(struct lc_point);
  nread *= sizeof(struct lc_point);

  /* Copy in */
  memcpy(rbuf, b->buf + offset, nread);

  return(0);
}
#else
/* This could do with a more intelligent, buffered, implementation */
static int buffer_read (struct buffer_info *b, struct lc_point *rbuf,
			off_t offset, int nread, char *errstr) {
  char *obuf;
  int nbread, nbdone, nbleft;

  int rv;
  off_t rvs;

  obuf = (char *) rbuf;

  /* Convert offset to bytes */
  offset *= sizeof(struct lc_point);

  /* Seek to the right place */
  rvs = lseek(b->fd, offset, SEEK_SET);
  if(rvs == (off_t) -1) {
    report_syserr(errstr, "lseek");
    goto error;
  }

  /* How many bytes to do? */
  nbread = nread * sizeof(struct lc_point);
  nbdone = 0;
  
  while(nbdone < nbread) {
    nbleft = nbread - nbdone;

    /* Read */
  do_read:
    rv = read(b->fd, obuf + nbdone, nbleft);
    if(rv < 0) {
      if(errno == EINTR)
	goto do_read;
      else {
	report_syserr(errstr, "read");
	goto error;
      }
    }
    else if(rv == 0) {
      report_err(errstr, "unexpected EOF reading buffer");
      goto error;
    }
    
    nbdone += rv;
  }

  return(0);

 error:
  return(1);
}
#endif

#ifdef HAVE_MMAP
static int buffer_write (struct buffer_info *b, struct lc_point *wbuf,
			 off_t offset, int nwrite, char *errstr) {
  /* Convert to bytes */
  offset *= sizeof(struct lc_point);
  nwrite *= sizeof(struct lc_point);

  /* Copy out */
  memcpy(b->buf + offset, wbuf, nwrite);

  return(0);
}
#else
static int buffer_write (struct buffer_info *b, struct lc_point *wbuf,
			 off_t offset, int nwrite, char *errstr) {
  char *ibuf;
  int nbwrite, nbwritten, nbleft;

  int rv;
  off_t rvs;

  ibuf = (char *) wbuf;

  /* Convert offset to bytes */
  offset *= sizeof(struct lc_point);

  /* Seek to the right place */
  rvs = lseek(b->fd, offset, SEEK_SET);
  if(rvs == (off_t) -1) {
    report_syserr(errstr, "lseek");
    goto error;
  }

  /* How many bytes to do? */
  nbwrite = nwrite * sizeof(struct lc_point);
  nbwritten = 0;
  
  while(nbwritten < nbwrite) {
    nbleft = nbwrite - nbwritten;

    /* Write */
  do_write:
    rv = write(b->fd, ibuf + nbwritten, nbleft);
    if(rv < 0) {
      if(errno == EINTR)
	goto do_write;
      else {
	report_syserr(errstr, "write");
	goto error;
      }
    }

    nbwritten += rv;
  }

  return(0);

 error:
  return(1);
}
#endif
