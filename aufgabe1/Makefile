CC = g++
CFLAGS = -Wall -g
obj = aufgabe1.o

all: $(obj)
	$(CC) $(CFLAGS) $(obj) -lhighs -lz -o aufgabe1

%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) -c $^

clean:
	rm *.o *.gch aufgabe1