CC = gcc
CFLAGS = -Wall -g 

all: visu.o a.out

visu.o: visu.h visu.c
	$(CC) $(CFLAGS) visu.c -c

a.out: test1.c visu.o
	$(CC) $(CFLAGS) test1.c visu.o -o $@

clean:
	rm -rf ./a.out *.o 

