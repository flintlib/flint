CFLAGS=-pedantic -Wall -g -O2

test: test.c
	gcc -o test test.c $(CFLAGS) -I/home/fredrik/src/arb -I/home/fredrik/src/arb/build/include -L/home/fredrik/src/arb -L/home/fredrik/src/flint2 -larb -lflint -lmpfr -lgmp

check:
	valgrind ./test

clean:
	rm test

