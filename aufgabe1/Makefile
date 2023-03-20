CC = g++
CFLAGS = -Wall -Wextra -O3 -march=native

all: ilp randomized

ilp:
	$(CC) $(CFLAGS) aufgabe1_ilp.cpp -lhighs -lz -o aufgabe1_ilp

randomized:
	$(CC) $(CFLAGS) aufgabe1_randomized.cpp -o aufgabe1_randomized

clean:
	rm *.o *.gch aufgabe1_ilp aufgabe1_randomized