# Compiler and flags
CC = gcc
CFLAGS = -Wall -O2 -g
LDFLAGS =

# Executable name
TARGET = program

# Source files and object files
SRCS = main.c butterfly.c data.c utils.c
OBJS = $(SRCS:.c=.o)

# Include directories
INCLUDES = -I.

# Default target
all: $(TARGET)

# Build the executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

# Compile source files to object files
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

# Clean up generated files
clean:
	rm -f $(OBJS) $(TARGET)

# Run the program
run: $(TARGET)
	./$(TARGET)

# Debugging options
debug: CFLAGS += -DDEBUG -g
debug: clean all
