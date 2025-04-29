CC = mpicc
CFLAGS = -Wall -O3 -g
BINS = quicksort

all: $(BINS)

quicksort: quicksort.h quicksort.c pivot.h pivot.c
	$(CC) $(CFLAGS) -o $@ quicksort.c pivot.c

clean:
	$(RM) $(BINS)

