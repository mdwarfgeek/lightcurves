# Makefile for lightcurves
#

### Constants: edit to suit your system ####

# CFITSIO include and library paths
CFITSIO_INC:=`pkg-config cfitsio --cflags`
CFITSIO_LIB:=`pkg-config cfitsio --libs`

# PGPLOT and X11 directories
PGPLOT_INC=-I$(PGPLOT_DIR)
PGPLOT_LIB=-L$(PGPLOT_DIR) -lcpgplot -lpgplot
X11_LIB=-L/usr/X11R6/lib -lX11

# SLALIB is now needed only for HJD (backwards compatibility) support.
# Uncomment these lines and adjust accordingly if you need it.
#SLA_INC=-I/usr/local/include -DHJD -DCSLALIB
#SLA_LIB=-L/usr/local/lib -lsla
#SLA_SRCS=hjd.c sla.c

# C compiler
#CC=gcc

# Opt/debug flags
#OPT=-g  # NO OPT FOR DEBUG
OPT=-g -O3 -ffast-math

# Compiler flags
CFLAGS=-std=gnu99 $(OPT) -Wall -I../lib $(CFITSIO_INC) $(PGPLOT_INC) $(SLA_INC) -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_MMAP

# Linker flags
LIBS=$(CFITSIO_LIB) $(PGPLOT_LIB) $(SLA_LIB) $(X11_LIB) -lX11 -lg2c -lm

#### End constants section ####

COMMON_SRCS=buffer.c chooseap.c instvers.c intra.c lightcurves.c \
	plots.c readfits.c systematic.c xytoxy.c \
	dsolve.c dmatinv.c filelist.c medsig.c sortfloat.c sortlong.c

LIGHTCURVES_SRCS=main.c $(COMMON_SRCS) $(SLA_SRCS)
LIGHTCURVES_OBJS=${LIGHTCURVES_SRCS:%.c=%.o}

UPDATE_SRCS=update.c $(COMMON_SRCS) $(SLA_SRCS)
UPDATE_OBJS=${UPDATE_SRCS:%.c=%.o}

TESTBUF_SRCS=testbuf.c buffer.c
TESTBUF_OBJS=${TESTBUF_SRCS:%.c=%.o}

LIB_OBJS=../lib/cvtunit.o ../lib/fitsutil.o ../lib/util.o ../lib/liblfa.a

# Rules for building C

.SUFFIXES: .c

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

all: lightcurves update

../lib/liblfa.a:
	(cd ../lib && $(MAKE) liblfa.a)

../lib/cvtunit.o:
	(cd ../lib && $(MAKE) cvtunit.o)

../lib/fitsutil.o:
	(cd ../lib && $(MAKE) fitsutil.o)

../lib/util.o:
	(cd ../lib && $(MAKE) util.o)

# Rules for lightcurves

lightcurves: $(LIGHTCURVES_OBJS) $(LIB_OBJS)
	$(CC) -o $@ $(LIGHTCURVES_OBJS) $(LIB_OBJS) $(LIBS)

update: $(UPDATE_OBJS) $(LIB_OBJS)
	$(CC) -o $@ $(UPDATE_OBJS) $(LIB_OBJS) $(LIBS)

testbuf: $(TESTBUF_OBJS) $(LIB_OBJS)
	$(CC) -o $@ $(TESTBUF_OBJS) $(LIB_OBJS) $(LIBS)

clean:
	rm -f $(LIGHTCURVES_OBJS) lightcurves
	rm -f $(UPDATE_OBJS) update
	rm -f $(TESTBUF_OBJS) testbuf
	rm -f *.core *~
