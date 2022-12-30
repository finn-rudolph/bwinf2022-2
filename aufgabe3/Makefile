CC = g++
CFLAGS = -Wall -g -O3 -march=native
a_objects = aufgabe3_a.o util.o
b_objects = aufgabe3_b.o util.o

all: a b

a: $(a_objects)
	$(CC) $(CFLAGS) $(a_objects) -o aufgabe3_a

b: $(b_objects)
	$(CC) $(CFLAGS) $(b_objects) -o aufgabe3_b

%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) -c $^

clean:
	rm *.o *.gch aufgabe3_a aufgabe3_b