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

# Compiler flags
CFLAGS=-g -O -Wall -I$(CFITSIO_INC) -I$(PGPLOT_INC) -I$(SLA_INC) -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DCSLALIB -DHAVE_MMAP
#CFLAGS=-xO4 -g -I$(CFITSIO_INC) -I$(PGPLOT_INC) -I$(SLA_INC) -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_MMAP -D__EXTENSIONS__ -D_LONGLONG_TYPE

# Linker flags
LIBS=-L$(CFITSIO_LIB) -lcfitsio -L$(PGPLOT_LIB) -lcpgplot -lpgplot -L$(SLA_LIB) -lsla -lg2c -L$(X11_LIB) -lX11 -lm
#LIBS=-L$(CFITSIO_LIB) -lcfitsio -L$(PGPLOT_LIB) -lcpgplot -lpgplot -L$(SLA_LIB) -lsla -L/usr/local/lib -lpng -lz -L$(X11_LIB) -lX11 -L/opt/SUNWspro/lib -R/opt/SUNWspro/lib -lF77 -lM77 -lsunmath -lm -lnsl -lsocket

#### End constants section ####

LIGHTCURVES_SRCS=main.c buffer.c chooseap.c intra.c lightcurves.c plots.c systematic.c \
	readfits.c hjd.c sla.c \
	dsolve.c hanning.c medsig.c sortfloat.c sortlong.c cvtunit.c util.c
LIGHTCURVES_OBJS=${LIGHTCURVES_SRCS:%.c=%.o}

TESTBUF_SRCS=testbuf.c buffer.c cvtunit.c util.c
TESTBUF_OBJS=${TESTBUF_SRCS:%.c=%.o}

# Rules for building C

.SUFFIXES: .c

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

all: lightcurves testbuf

# Rules for lightcurves

lightcurves: $(LIGHTCURVES_OBJS)
	$(CC) -o $@ $(LIGHTCURVES_OBJS) $(LIBS)

testbuf: $(TESTBUF_OBJS)
	$(CC) -o $@ $(TESTBUF_OBJS) $(LIBS)

clean:
	rm -f $(LIGHTCURVES_OBJS) lightcurves
	rm -f $(TESTBUF_OBJS) testbuf
	rm -f *.core *~
