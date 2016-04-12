# CC=g++
CC=/usr/bin/g++
CFLAGS= -Wall -Ofast -std=c++11 -march=native -pthread -Wextra -flto
LDFLAGS=-pthread

ifeq ($(gprof),1)
CFLAGS= -Wall -Ofast -std=c++11 -march=native -pthread -Wextra -pg
LDFLAGS=-pg -pthread
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif

EXEC=blong

all: $(EXEC)

Utils.o: Utils.cpp Utils.h
	$(CC) -o $@ -c $< $(CFLAGS)
	
graph.o: graph.cpp graph.h
	$(CC) -o $@ -c $< $(CFLAGS)

MappingSupervisor.o: MappingSupervisor.cpp MappingSupervisor.h Utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

main.o: main.cpp MappingSupervisor.h Utils.h graph.h
	$(CC) -o $@ -c $< $(CFLAGS)

blong: main.o MappingSupervisor.o Utils.o graph.o
	$(CC) -o $@ $^ $(LDFLAGS)



clean:
	rm -rf *.o
	rm -rf $(EXEC)

rebuild: clean $(EXEC)
