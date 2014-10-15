# Makefile for lightcurves
#

### Constants: edit to suit your system ####

# CFITSIO include and library paths
CFITSIO_INC=`pkg-config cfitsio --cflags`
CFITSIO_LIBS=`pkg-config cfitsio --libs`

# PGPLOT is is now optional.  Comment out these lines to disable plots
# and remove the dependency.
PGPLOT_DIR?=/usr/local/pgplot
PGPLOT_INC=-I$(PGPLOT_DIR) -DPLOTS
PGPLOT_LIBS=-L$(PGPLOT_DIR) -lcpgplot -lpgplot -L/usr/X11R6/lib -lX11
PGPLOT_SRCS=plots.c

# SLALIB is now needed only for HJD (backwards compatibility) support.
# Uncomment these lines and adjust accordingly if you need it.
#SLA_INC=-I/usr/local/include -DHJD
#SLA_LIBS=-L/usr/local/lib -lsla
#SLA_SRCS=hjd.c sla.c

# C compiler
#CC=gcc

# Fortran compiler
#FC=gfortran

# Linking is done with the Fortran compiler in the standard configuration
# because PGPLOT is a Fortran library, and this seems to be the most
# portable way to bring in all the runtime libraries.  However, if plots
# are disabled above, the C compiler can be used instead by changing this
# line to CC (e.g. if there is no Fortran compiler available).
LINK=$(FC)

# Optimization flags

# DEBUG:
#OPT=-g

# OPT:
OPT=-g -O3 -ffast-math

# Compiler flags
CFLAGS=-std=gnu99 $(OPT) -Wall -I../lib $(CFITSIO_INC) $(PGPLOT_INC) $(SLA_INC) -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_MMAP

# Linker flags
LIBS=$(CFITSIO_LIBS) $(PGPLOT_LIBS) $(SLA_LIBS) -lm

#### End constants section ####

COMMON_SRCS=buffer.c chooseap.c instvers.c intra.c lightcurves.c \
	readfits.c systematic.c xytoxy.c \
	dsolve.c dmatinv.c filelist.c

LIGHTCURVES_SRCS=main.c $(COMMON_SRCS) $(SLA_SRCS) $(PGPLOT_SRCS)
LIGHTCURVES_OBJS=${LIGHTCURVES_SRCS:%.c=%.o}

UPDATE_SRCS=update.c $(COMMON_SRCS) $(SLA_SRCS) $(PGPLOT_SRCS)
UPDATE_OBJS=${UPDATE_SRCS:%.c=%.o}

TESTBUF_SRCS=testbuf.c buffer.c
TESTBUF_OBJS=${TESTBUF_SRCS:%.c=%.o}

DEPEND_SRCS=main.c update.c testbuf.c $(COMMON_SRCS) $(SLA_SRCS) $(PGPLOT_SRCS)

LIB_OBJS=../lib/fitsutil.o ../lib/liblfa.a

# Rules for building C

.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

all: lightcurves update

depend:
	$(CC) $(CFLAGS) -E -MM $(DEPEND_SRCS) > .depend

../lib/liblfa.a:
	(cd ../lib && $(MAKE) liblfa.a)

../lib/fitsutil.o:
	(cd ../lib && $(MAKE) fitsutil.o CFITSIO_INC=$(CFITSIO_INC))

# Rules for lightcurves

lightcurves: $(LIGHTCURVES_OBJS) $(LIB_OBJS)
	$(LINK) -o $@ $(LIGHTCURVES_OBJS) $(LIB_OBJS) $(LIBS)

update: $(UPDATE_OBJS) $(LIB_OBJS)
	$(LINK) -o $@ $(UPDATE_OBJS) $(LIB_OBJS) $(LIBS)

testbuf: $(TESTBUF_OBJS) $(LIB_OBJS)
	$(LINK) -o $@ $(TESTBUF_OBJS) $(LIB_OBJS) $(LIBS)

clean:
	rm -f $(LIGHTCURVES_OBJS) lightcurves
	rm -f $(UPDATE_OBJS) update
	rm -f $(TESTBUF_OBJS) testbuf
	rm -f *.core *~
