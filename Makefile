# Makefile for compiling the RLZ program.
# Author: Shanika Kuruppu
# Date: 5/08/2011

CC = g++

CFLAGS = -Wall -O9

INCLUDES = -I/home/data/kuruppu/src/lib/cds/include

LFLAGS = -L/home/data/kuruppu/src/lib/cds/lib

LIBS = -lcds -ldivsufsort64

SRCS = RLZ.cpp Bits.cpp main.cpp

OBJS = RLZ.o Bits.o main.o

all: rlz

rlz: $(OBJS) 
	$(CC) $(CFLAGS) -o rlz $(OBJS) $(INCLUDES) $(LFLAGS) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

clean:
	rm -f $(OBJS)

clobber: clean
	rm -f rlz

depend: $(SRCS)
	$(CC) $(INCLUDES) -MM $(SRCS) > depend

include depend
