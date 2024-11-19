CC = gcc
CFLAGS = -Wall -Wextra -O2
OBJS = test.o poly_mul.o data.o

all: test

test: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o test

test.o: test.c data.h
	$(CC) $(CFLAGS) -c test.c

poly_mul.o: poly_mul.c poly_mul.h data.h
	$(CC) $(CFLAGS) -c poly_mul.c

data.o: data.c data.h
	$(CC) $(CFLAGS) -c data.c

clean:
	rm -f *.o test
