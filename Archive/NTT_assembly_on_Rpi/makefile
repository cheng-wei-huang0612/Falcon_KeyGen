cat makefile
# Compiler and flags
CC = gcc
AS = as
CFLAGS = -Wall -O2
TARGET = program

# Source and object files
C_SRCS = main.c data.c hal-cortexa.c
ASM_SRCS = butterfly_ntt.s
OBJS = $(C_SRCS:.c=.o) $(ASM_SRCS:.s=.o)

# Build target
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET)

# Compile C source files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Assemble assembly source files
%.o: %.s
	$(AS) -g $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)
