lightcurves
===========

This program combines list-driven photometry files from the CASU
pipeline toolkit or CASUTOOLS into light curves, optionally performing
differential photometry.  Output is written to FITS binary tables.

Most of the documentation is in the form of a man page.  To view this
without installing it:

  nroff -man -Tascii lightcurves.man | more

GNU "less" produces slightly nicer formatting, if available.

Dependencies
------------

* A C99 compiler on a POSIX / UNIX-like platform.

The use of the mmap system call can be disabled if needed, but there
are other dependencies and assumptions specific to this type of
system, such as the use of getopt and isatty, which are POSIX
functions.  I may be motivated to fix this if it is a serious problem.

The program has been successfully compiled and run on Windows using
MinGW.  Native Windows API calls are used in buffer.c to perform
memory mapping if HAVE_MMAP is defined (this is the default in the
Makefile).  Full functionality has not been tested thoroughly so some
issues may remain.

* My C subroutine library from "lib" on github.

The Makefiles currently assume it's located in a directory "../lib".
Please refer to the README.md file there for notes on how to install
the IERS data files and JPL ephemerides.  These are needed to compute
Barycentric Julian Dates (BJDs) in the output file.  The program will
run without them, but the BJD values will not be written.

* CFITSIO, version 3.x
http://heasarc.gsfc.nasa.gov/fitsio/

The Makefile uses pkg-config to obtain suitable compiler and linker
flags for CFITSIO.  On some systems, if CFITSIO is installed in a
non-standard location, the PKG_CONFIG_PATH environment variable may
need to be adjusted to allow pkg-config to find the cfitsio.pc file
installed by CFITSIO.

* PGPLOT >= 5.1.1 (optional, plots disabled if unavailable)
http://www.astro.caltech.edu/~tjp/pgplot/

To disable use of PGPLOT, comment out the appropriate lines in the
Makefile.  This will result in no diagnostic plots being produced, but
this is not a serious problem, most of the quantities needed to
produce them appear in the output files.  The most important of these
is an RMS versus magnitude diagram, and can be produced by plotting
the "rms" versus "medflux" columns from the output file using your
favourite method for plotting columns of FITS binary tables.

The Fortran compiler (FC variable in "make") is used to link the
binaries in order to bring in any runtime dependencies for PGPLOT.
This can be changed to use the C compiler if PGPLOT is not being
linked.  Please refer to the appropriate commented section of the
Makefile.

CASUTOOLS (publicly available) or the CASU pipeline toolkit (not
publicly available) are needed to produce suitable input files, but
these packages are not directly needed to compile, link, or run
anything in this package.

For CASUTOOLS, see:
http://casu.ast.cam.ac.uk/surveys-projects/software-release

Building
--------

Type "make" to build both binaries.  There is no install target, but
this can be done if desired by copying the binaries "lightcurves" and
"update" to a suitable location, and "lightcurves.man" can be placed
in one of the man1 directories on the man path so the standard "man"
utility can find it.  If installing these binaries on the system, it
would be advisable to rename "update" to a more descriptive name, less
likely to clash with other packages, such as "lightcurves-update".

The standard Makefile should work unmodified if:

(a) pkg-config can find cfitsio.pc (see above).
(b) the PGPLOT_DIR environment variable points to your PGPLOT
    installation, or PGPLOT is installed in one of the system include
    and library paths.
(c) the FC variable in "make" points to the correct Fortran compiler.

Invariably, one of these things is not true for most real systems,
especially in astronomy departments.

Item (c) deserves a special note.  Most "make" programs I have come
across recently default to FC=f77, which is almost always the wrong
choice on modern systems (in fact, it often does not even exist!).
On most Linux and other systems using the GNU compiler collection, you
will need to set FC=gfortran.  This can be done on the "make" command
line or with an environment variable.

