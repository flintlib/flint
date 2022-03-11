CFLAGS=-pedantic -Wall -g -O2 -I. -I/home/fredrik/src/arb -I/home/fredrik/src/arb/build/include -L/home/fredrik/src/arb -L/home/fredrik/src/flint2 -larb -lflint -lmpfr -lgmp

test: test.c
	gcc -o test test.c gr/*.c $(CFLAGS)

check:
	./test

valgrind:
	valgrind ./test

clean:
	rm test

