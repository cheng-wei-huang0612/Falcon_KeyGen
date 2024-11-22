# Compiler and flags
CC = gcc
CFLAGS = -O2 -Wall

# Target executable
TARGET = program

# Source files (excluding subtest.c)
SRCS = main.c data.c utils.c
# Object files
OBJS = $(SRCS:.c=.o)

# Build target
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

# Rule for object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean
