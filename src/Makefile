# Makefile to compile the following C++ files:
# - UMDcompressor.cpp to create the executable 'umfgen'
# - UMFParser.cpp to create the executable 'umfparser-live'
# - UMFQuery.cpp to create the executable 'umfparser'
# - UMFChecker.cpp to create the executable 'umfchecker'

# Compiler and flags
CXX = g++
CXXFLAGS = -g -w -O3 -I.
LDFLAGS = 

# Source files and corresponding executables
SOURCES = UMDcompressor.cpp UMFParser.cpp UMFQuery.cpp UMFChecker.cpp
EXECUTABLES = umfgen umfparser-live umfparser umfchecker

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Default target to build all executables
all: $(EXECUTABLES)

# Rule to build each executable from its corresponding source file
umfgen: UMDcompressor.o
	$(CXX) $(LDFLAGS) -o $@ $^
umfparser-live: UMFParser.o
	$(CXX) $(LDFLAGS) -o $@ $^
umfparser: UMFQuery.o
	$(CXX) $(LDFLAGS) -o $@ $^
umfchecker: UMFChecker.o
	$(CXX) $(LDFLAGS) -o $@ $^

# Rule to compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean target to remove executables and object files
clean:
	rm -f $(EXECUTABLES) $(OBJECTS)

.PHONY: all clean