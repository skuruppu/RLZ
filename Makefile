# Makefile for compiling the RLZ program.
# Author: Shanika Kuruppu
# Date: 5/08/2011

CC = g++

CFLAGS = -Wall -O9 -g

INCLUDES = -I../lib/cds/include

LFLAGS = -L../lib/cds/lib

LIBS = -lcds -ldivsufsort64

SRCS = RLZ.cpp Bits.cpp main.cpp RLZ_index.cpp

OBJS = RLZ.o Bits.o main.o

IDXOBJS = RLZ_index.o Bits.o

all: rlz

rlz: $(OBJS) 
	$(CC) $(CFLAGS) -o rlz $(OBJS) $(INCLUDES) $(LFLAGS) $(LIBS)

rlzindex: $(IDXOBJS)
	$(CC) $(CFLAGS) -o $@ $(IDXOBJS) $(INCLUDES) $(LFLAGS) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

clean:
	rm -f $(OBJS)

clobber: clean
	rm -f rlz

depend: $(SRCS)
	$(CC) $(INCLUDES) -MM $(SRCS) > depend

include depend
