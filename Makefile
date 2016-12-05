CC=gcc

all:
	$(CC) -o membrane membrane.c -O3 -lm
