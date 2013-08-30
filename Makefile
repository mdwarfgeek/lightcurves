# Makefile for lightcurves
#

### Constants: edit to suit your system ####

# CFITSIO include and library paths
CFITSIO_INC=/usr/local/include
CFITSIO_LIB=/usr/local/lib

# PGPLOT and X11 directories
PGPLOT_INC=/usr/local/pgplot
PGPLOT_LIB=/usr/local/pgplot
X11_LIB=/usr/X11R6/lib

# SLALIB include and library paths
SLA_INC=/usr/local/include
SLA_LIB=/usr/local/lib

# C compiler
#CC=gcc

# Opt/debug flags
#OPT=-g  # NO OPT FOR DEBUG
OPT=-g -O3 -ffast-math

# Compiler flags
CFLAGS=$(OPT) -Wall -I$(CFITSIO_INC) -I$(PGPLOT_INC) -I$(SLA_INC) -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DCSLALIB -DHAVE_MMAP

# Linker flags
LIBS=-L$(CFITSIO_LIB) -lcfitsio -L$(PGPLOT_LIB) -lcpgplot -lpgplot -L$(SLA_LIB) -lsla -lg2c -L$(X11_LIB) -lX11 -lm

#### End constants section ####

LIGHTCURVES_SRCS=main.c buffer.c chooseap.c instvers.c intra.c lightcurves.c plots.c systematic.c xytoxy.c \
	readfits.c hjd.c sla.c \
	dsolve.c dmatinv.c hanning.c medsig.c sortfloat.c sortlong.c cvtunit.c util.c
LIGHTCURVES_OBJS=${LIGHTCURVES_SRCS:%.c=%.o}

UPDATE_SRCS=update.c buffer.c chooseap.c instvers.c intra.c lightcurves.c plots.c systematic.c xytoxy.c \
	readfits.c hjd.c sla.c \
	dsolve.c dmatinv.c hanning.c medsig.c sortfloat.c sortlong.c cvtunit.c util.c
UPDATE_OBJS=${UPDATE_SRCS:%.c=%.o}

TESTBUF_SRCS=testbuf.c buffer.c cvtunit.c util.c
TESTBUF_OBJS=${TESTBUF_SRCS:%.c=%.o}

# Rules for building C

.SUFFIXES: .c

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

all: lightcurves update

# Rules for lightcurves

lightcurves: $(LIGHTCURVES_OBJS)
	$(CC) -o $@ $(LIGHTCURVES_OBJS) $(LIBS)

update: $(UPDATE_OBJS)
	$(CC) -o $@ $(UPDATE_OBJS) $(LIBS)

testbuf: $(TESTBUF_OBJS)
	$(CC) -o $@ $(TESTBUF_OBJS) $(LIBS)

clean:
	rm -f $(LIGHTCURVES_OBJS) lightcurves
	rm -f $(UPDATE_OBJS) update
	rm -f $(TESTBUF_OBJS) testbuf
	rm -f *.core *~
