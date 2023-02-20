CC = gcc
CFLAGS = -Wall -Wextra -g -fsanitize=address,undefined -D_GLIBCXX_DEBUG
obj = aufgabe1.o

all: $(obj)
	$(CC) $(CFLAGS) $(obj) -lglpk -lgmp -lm -o aufgabe1

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $^

clean:
	rm *.o *.gch aufgabe1