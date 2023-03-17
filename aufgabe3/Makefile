CC = g++
CFLAGS = -Wall -Wextra -std=c++20 -O3 -march=native -DNDEBUG
obj_a = aufgabe3_a.o util.o
obj_b = aufgabe3_b.o util.o

all: a b

a: $(obj_a)
	$(CC) $(CFLAGS) $(obj_a) -o aufgabe3_a

b: $(obj_b)
	$(CC) $(CFLAGS) $(obj_b) -o aufgabe3_b

%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) -c $^

clean:
	rm *.o *.gch aufgabe3_a aufgabe3_b