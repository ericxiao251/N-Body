CC=mpiCC
CFLAGS=-c -g -Wall -std=c++11
LDFLAGS=-lm
OBJECTS=$(SOURCES:.cpp=.o)

SOURCES=main.cpp image.c util.c cyclic_dist.c physics.c

EXECUTABLE=project.x

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o *.x
