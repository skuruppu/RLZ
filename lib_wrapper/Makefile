# Makefile for compiling the data structure library wrapper.
# Author: Shanika Kuruppu
# Date: 4/05/2012

CC = g++

CFLAGS = -Wall -O9 -g

SRCS = cds_wrapper.cpp

OBJS = wrapper.o cds_wrapper.o

all: $(OBJS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJS)

depend: $(SRCS)
	$(CC) -MM $(SRCS) > depend

include depend

.PHONY: all
