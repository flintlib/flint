CFLAGS=-ansi -pedantic -Wall -g -O2 -I. -I/home/fredrik/src/calcium -I/home/fredrik/src/arb -I/home/fredrik/src/arb/build/include -L/home/fredrik/src/calcium -L/home/fredrik/src/arb -L/home/fredrik/src/flint2 -lcalcium -larb -lflint -lmpfr -lgmp

test: test.c
	gcc -o test test.c gr/*.c gr_vec/*.c gr_mat/*.c gr_poly/*.c $(CFLAGS)

check:
	./test

valgrind:
	valgrind ./test

clean:
	rm test

