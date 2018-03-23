CC=g++

CFLAGS=-std=c++11 -w

DIR=/Users/alok/Setups/bamtools-master
IDIR=$(DIR)/include
LDIR=$(DIR)/lib
LIBS=-lbamtools -ltbb

all: indel

indel:
	$(CC) $(CFLAGS) -I$(IDIR) -L$(LDIR) src/Cluster.cpp -o indel $(LIBS)

clean:
	rm indel
