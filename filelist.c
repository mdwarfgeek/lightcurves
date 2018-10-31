#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifndef _WIN32
/* No glob on Win32, but we can get it done for us by linking with
   setargv.obj.  On my MinGW install this (or some equivalent) seems
   to be the default anyway. */
#include <glob.h>
#endif

#ifdef ZIPSUPPORT
#include <zip.h>
#endif

#include "lightcurves.h"
#include "util.h"

static int is_zip_file (char *filename, char *errstr);

struct input_file *read_file_list (int argc, char **argv,
                                   int *nf_r, char *errstr) {
  struct input_file *list = (struct input_file *) NULL;

  FILE *fp = (FILE *) NULL;
  char line[16384], *av, *p, *s;

  int a, f, nf = 0;

#ifndef _WIN32
  glob_t gbuf;
  int rv, op;

  /* Init */
  gbuf.gl_pathv = (char **) NULL;
#endif

  int is_zip;

#ifdef ZIPSUPPORT
  zip_t *z = (zip_t *) NULL;
  zip_error_t ze;
  struct zip_stat sb;

  int zerrno;
  int ent, nent;
  const char *name;
  int namelen;

  zip_uint8_t opsys;
  zip_uint32_t extattr;
#endif

  /* Loop over arguments */
  for(a = 0; a < argc; a++) {
    av = argv[a];

    if(*av == '@') {  /* @list form */
      av++;  /* skip past the @ */

      is_zip = is_zip_file(av, errstr);
      if(is_zip < 0)
        goto error;

      if(is_zip) {
#ifdef ZIPSUPPORT
        z = zip_open(av, 0, &zerrno);
        if(!z) {
          zip_error_init_with_code(&ze, zerrno);
          report_err(errstr, "zip_open: %s: %s",
                     av, zip_error_strerror(&ze));
          zip_error_fini(&ze);
          goto error;
        }
        
        nent = zip_get_num_entries(z, 0);
        if(nent < 0) {
          report_err(errstr, "zip_get_num_entries: %s",
                     zip_strerror(z));
          goto error;
        }
        
        for(ent = 0; ent < nent; ent++) {
          if(zip_stat_index(z, ent, 0, &sb)) {
            report_err(errstr, "zip_stat_index(%d): %s\n",
                       ent, zip_strerror(z));
            goto error;
          }
          
          if(zip_file_get_external_attributes(z, ent, 0, &opsys, &extattr)) {
            report_err(errstr, "zip_get_external_attributes(%d): %s\n",
                       ent, zip_strerror(z));
            goto error;
          }
          
          name = sb.name;
          namelen = strlen(name);
          
          if((opsys != ZIP_OPSYS_AMIGA && extattr & 0x10) ||
             (namelen > 0 && name[namelen-1] == '/'))
            /* Is directory */
            continue;
          
          /* Otherwise, add to list */
          list = (struct input_file *) realloc(list, (nf+1) * sizeof(struct input_file));
          if(!list) {
            report_syserr(errstr, "realloc");
            goto error;
          }
          
          s = strdup(name);
          if(!s) {
            report_syserr(errstr, "strdup");
            goto error;
          }
          
          list[nf].filename = s;
          list[nf].ient = ent;
          list[nf].iarg = a;
          list[nf].arg = av;
          nf++;
        }
        
        if(zip_close(z)) {
          report_err(errstr, "zip_close: %s",
                     zip_strerror(z));
          goto error;
        }

        z = (zip_t *) NULL;
#else
        report_err(errstr, "not compiled with zip support");
        goto error;
#endif
      }
      else {
        fp = fopen(av, "r");
        if(!fp) {
          report_syserr(errstr, "open: %s", av);
          goto error;
        }
        
        while(fgets(line, sizeof(line), fp)) {
          p = sstrip(line);
          
          if(*p) {
            list = (struct input_file *) realloc(list, (nf+1) * sizeof(struct input_file));
            if(!list) {
              report_syserr(errstr, "realloc");
              goto error;
            }
            
            s = strdup(p);
            if(!s) {
              report_syserr(errstr, "strdup");
              goto error;
            }
            
            list[nf].filename = s;
            list[nf].ient = -1;
            list[nf].iarg = a;
            list[nf].arg = av;
            nf++;
          }
        }
        
        if(ferror(fp)) {
          report_syserr(errstr, "read: %s", av);
          goto error;
        }
        
        fclose(fp);
        fp = (FILE *) NULL;
      }
    }
    else {  /* glob or single filename */
#ifndef _WIN32
      rv = glob(av, 0, NULL, &gbuf);
      if(rv == 0) {  /* succeeded */
        /* Allocate block */
        list = (struct input_file *) realloc(list,
                                             (nf+gbuf.gl_pathc) * sizeof(struct input_file));
        if(!list) {
          report_syserr(errstr, "realloc");
          goto error;
        }

        /* Zero out the new pointers */
        memset(list+nf, 0, gbuf.gl_pathc * sizeof(struct input_file));

        /* Record so we know to free them */
        op = nf;

        nf += gbuf.gl_pathc;
        
        for(f = 0; f < gbuf.gl_pathc; f++) {
          s = strdup(gbuf.gl_pathv[f]);
          if(!s) {
            report_syserr(errstr, "strdup");
            goto error;
          }

          list[op].filename = s;
          list[op].ient = -1;
          list[op].iarg = a;
          list[op].arg = av;
          op++;
        }

        globfree(&gbuf);
        gbuf.gl_pathv = (char **) NULL;
      }
      else if(rv == GLOB_NOMATCH) {  /* no match */
#endif  /* _WIN32 */

        /* Assume it's a single entry. */
        list = (struct input_file *) realloc(list, (nf+1) * sizeof(struct input_file));
        if(!list) {
          report_syserr(errstr, "realloc");
          goto error;
        }
        
        s = strdup(av);
        if(!s) {
          report_syserr(errstr, "strdup");
          goto error;
        }
        
        list[nf].filename = s;
        list[nf].ient = -1;
        list[nf].iarg = a;
        list[nf].arg = av;
        nf++;
#ifndef _WIN32
      }
      else {
        report_err(errstr, "glob error: %s", av);
        goto error;
      }
#endif
    }
  }

  *nf_r = nf;

  return(list);

 error:
  if(fp)
    fclose(fp);

#ifndef _WIN32
  if(gbuf.gl_pathv)
    globfree(&gbuf);
#endif

#ifdef ZIPSUPPORT
  if(z) {
    zip_close(z);
    z = (zip_t *) NULL;
  }
#endif

  if(list) {
    for(f = 0; f < nf; f++)
      if(list[f].filename)
        free((void *) list[f].filename);

    free((void *) list);
  }

  return((struct input_file *) NULL);
}

#define NMAGIC 4

static int is_zip_file (char *filename, char *errstr) {
  FILE *fp;
  char buf[NMAGIC];
  int nread, rv;

  fp = fopen(filename, "rb");
  if(!fp) {
    report_syserr(errstr, "open: %s", filename);
    goto error;
  }

  /* Somewhat simplified method for detection of a zip file.  This tests
     for the signatures:
     PK\x03\x04  local file header (normal file)
     PK\x05\x06  end of central directory (empty file)
     PK\x07\x08  multi-volume archive
     which should cover most normal zip files, but strictly speaking we
     should detect zip format by seeking to end of file and backing up
     looking for the end of central directory record. */
  nread = fread(buf, 1, NMAGIC, fp);

  if(nread == NMAGIC &&
     (buf[0] == 'P' && buf[1] == 'K' &&
      ((buf[2] == '\x03' && buf[3] == '\x04') ||
       (buf[2] == '\x05' && buf[3] == '\x06') ||
       (buf[2] == '\x07' && buf[3] == '\x08'))))
    rv = 1;
  else if(ferror(fp))
    goto error;
  else
    rv = 0;

  fclose(fp);

  return(rv);

 error:

  return(-1);
}

