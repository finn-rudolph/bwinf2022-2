CC = g++
CFLAGS = -Wall -Wextra -O3 -march=native

all: ip randomized

ip:
	$(CC) $(CFLAGS) aufgabe1_ip.cpp -lhighs -lz -o aufgabe1_ip

randomized:
	$(CC) $(CFLAGS) aufgabe1_randomized.cpp -o aufgabe1_randomized

clean:
	rm *.o *.gch aufgabe1_ilp aufgabe1_randomized