CC := gcc
CFLAGS := -g -I$(PWD)/htslib/htslib
LDFLAGS := htslib/libhts.a -lz -lm -lpthread

.PHONY: htslib clean

all: htslib count_in_bins

htslib:
	cd htslib && $(MAKE)
count_in_bins: runtime.o time.o errors.o memalloc.o files.o slinklist.o sequences.o hashtable.o count_in_bins.o

runtime.o: runtime.c runtime.h
time.o: runtime.o time.c time.h asserts.h
errors.o : errors.c errors.h asserts.h 
memalloc.o: errors.o memalloc.c memalloc.h constants.h
files.o: memalloc.o files.c files.h asserts.h 
hashtable.o: memalloc.o slinklist.o hashtable.c hashtable.h asserts.h
slinklist.o: memalloc.o slinklist.c slinklist.h asserts.h
sequences.o: memalloc.o files.o sequences.c sequences.h asserts.h constants.h
count_in_bins.o: runtime.o time.o errors.o memalloc.o files.o slinklist.o sequences.o hashtable.o count_in_bins.c

clean:
	cd htslib && $(MAKE) clean
	rm *.o
	rm count_in_bins
