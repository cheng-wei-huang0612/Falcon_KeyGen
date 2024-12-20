CC = gcc
CFLAGS = -O2 -Wall -Wextra -lm
TARGET = main

all: data.c $(TARGET)

data.c data.h: datagen.py
	python3 datagen.py

$(TARGET): main.o data.o
	$(CC) $(CFLAGS) -o $(TARGET) main.o data.o

main.o: main.c data.h
	$(CC) $(CFLAGS) -c main.c

data.o: data.c data.h
	$(CC) $(CFLAGS) -c data.c

clean:
	rm -f *.o $(TARGET) data.c data.h
