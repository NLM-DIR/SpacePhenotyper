CPPFLAGS=-std=c++23 -O3 -Wall -Wextra -pthread
DEBUGFLAGS =-D_GLIBCXX_DEBUG

all: objects
	g++ main.cpp DataMatrix.o -o main $(CPPFLAGS) $(DEBUGFLAGS)

objects:
	g++ DataMatrix.cpp -c -o DataMatrix.o $(CPPFLAGS) $(DEBUGFLAGS)
