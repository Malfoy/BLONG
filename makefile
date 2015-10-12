CC=g++
CFLAGS= -Wall -O3 -std=c++11 -march=native -pthread -pedantic -Wextra
LDFLAGS=-pthread

ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O4 -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif

EXEC=blong

all: $(EXEC)

Utils.o: Utils.cpp Utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

nw.o: nw.cpp nw.h
	$(CC) -o $@ -c $< $(CFLAGS)

graph.o: graph.cpp graph.h
	$(CC) -o $@ -c $< $(CFLAGS)

binSeq.o: binSeq.cpp binSeq.h
	$(CC) -o $@ -c $< $(CFLAGS)

MappingSupervisor.o: MappingSupervisor.cpp MappingSupervisor.h Utils.h nw.h
	$(CC) -o $@ -c $< $(CFLAGS)

main.o: main.cpp MappingSupervisor.h Utils.h graph.h binSeq.h
	$(CC) -o $@ -c $< $(CFLAGS)

blong: main.o MappingSupervisor.o Utils.o graph.o binSeq.o nw.o
	$(CC) -o $@ $^ $(LDFLAGS)



clean:
	rm -rf *.o
	rm -rf $(EXEC)

rebuild: clean $(EXEC)

