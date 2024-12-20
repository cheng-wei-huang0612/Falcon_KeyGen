# Compiler and flags
CC = gcc
CFLAGS = -Wall -O2
ASM = as
ASMFLAGS = -g

# Target executable
TARGET = program

# Source files
C_SOURCES = main.c data.c
ASM_SOURCES = butterfly_ntt.s

# Object files (generated from C and assembly files)
OBJECTS = $(C_SOURCES:.c=.o) $(ASM_SOURCES:.s=.o)

# Default target to build the executable
all: $(TARGET)

# Rule to link object files and create the executable
$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET)

# Rule to compile C source files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to assemble butterfly_ntt.s
butterfly_ntt.o: butterfly_ntt.s
	$(ASM) $(ASMFLAGS) butterfly_ntt.s -o butterfly_ntt.o

# Clean up build files
clean:
	rm -f $(OBJECTS) $(TARGET)

# Print debug information
debug: CFLAGS += -g
debug: clean all
