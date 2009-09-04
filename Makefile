CC=gcc
CFLAGS=-std=c99 -O2
LIBS=-L$(CURDIR) -L/home/wbhart/mpir-trunk/.libs -lmpir -lm
INCS=-I$(CURDIR) -I/home/wbhart/mpir-trunk/
export

SOURCES = $(wildcard *.c)

HEADERS = $(wildcard *.h)

OBJS = $(patsubst %.c, %.o, $(SOURCES))

LOBJS = $(patsubst %.c, %.lo, $(SOURCES))

LIB_SOURCES = $(SOURCES) $(foreach dir, $(BUILD_DIRS), $(wildcard $(dir)/*.c))

LIB_OBJS = $(patsubst %.c, %.lo, $(LIB_SOURCES))

TEST_SOURCES = $(wildcard test/*.c)

PROF_SOURCES = $(wildcard profile/*.c)

TESTS = $(patsubst %.c, %, $(TEST_SOURCES))

PROFS = $(patsubst %.c, %, $(PROF_SOURCES))

all: $(OBJS) recursive 

clean:
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir) clean;)
	rm -f $(OBJS) $(LOBJS) $(TESTS) 

profile: all profiler.o
	$(foreach prog, $(PROFS), $(CC) -O2 -std=c99 $(INCS) $(prog).c profiler.o -o $(prog) $(LIBS) -lflint;)
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir) profile;)

recursive:
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir);) 

check: all $(LOBJS) library
	$(foreach prog, $(TESTS), $(CC) $(CFLAGS) $(INCS) $(prog).c -o $(prog) $(LIBS) -lflint;)
	$(foreach prog, $(TESTS), $(prog);)
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir) check;)

library: library-recursive $(LIB_OBJS)
	$(CC) -fPIC -shared $(LIB_OBJS) -o libflint.so

library-recursive:
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir) library;) 

.PHONY: profile library library-recursive recursive clean check check-recursive all

%.lo: %.c
	$(CC) -fPIC $(CFLAGS) $(INCS) -c $< -o $@

%.o: %.c
	$(CC) -fPIC $(CFLAGS) $(INCS) -c $< -o $@

BUILD_DIRS = ulong_extras fmpz

