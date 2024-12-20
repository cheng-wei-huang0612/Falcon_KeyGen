CC = gcc
CFLAGS = -O2 -Wall -Wextra
TARGET = main


all: $(TARGET)

$(TARGET): main.o data.o
	$(CC) $(CFLAGS) -o $(TARGET) main.o data.o

main.o: main.c data.h
	$(CC) $(CFLAGS) -c main.c

data.o: data.c data.h
	$(CC) $(CFLAGS) -c data.c

data: data.c data.h

clean:
	rm -f *.o $(TARGET)
