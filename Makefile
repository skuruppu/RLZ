# Makefile for compiling the RLZ program.
# Author: Shanika Kuruppu
# Date: 5/08/2011

CC = g++

CFLAGS = -Wall -O9 -g

LIBS = -lcds -ldivsufsort64

SRCS = RLZ.cpp Bits.cpp main.cpp RLZ_index.cpp

OBJS = RLZ.o Bits.o main.o lib_wrapper/wrapper.o

IDXOBJS = RLZ_index.o Bits.o

all: rlz rlzindex

rlz: wrapper $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)

rlzindex: $(IDXOBJS)
	$(CC) $(CFLAGS) -o $@ $(IDXOBJS) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< 

wrapper:
	$(MAKE) -C lib_wrapper

clean:
	rm -f $(OBJS) $(IDXOBJS)

clobber: clean
	rm -f rlz rlzindex

depend: $(SRCS)
	$(CC) -MM $(SRCS) > depend

include depend
