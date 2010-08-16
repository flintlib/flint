CFLAGS=-std=c99 -O2 
LIBS=-L$(CURDIR) -L/home/wbhart/lib/ -lflint -lmpir -lmpfr -lm
LIBS2=-L$(CURDIR) -L/home/wbhart/lib/ -lmpir -lmpfr -lm
INCS=-I$(CURDIR) -I/home/wbhart/include
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
	$(foreach prog, $(PROFS), $(CC) -O2 -std=c99 $(INCS) $(prog).c profiler.o -o $(prog) $(LIBS);)
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir) profile;)

recursive:
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir);) 

check: all $(LOBJS) library
ifndef MOD
	$(foreach prog, $(TESTS), $(CC) $(CFLAGS) $(INCS) $(prog).c -o $(prog) $(LIBS);)
	$(foreach prog, $(TESTS), $(prog);)
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir) check;)
else
	$(foreach dir, $(MOD), $(MAKE) -C $(dir) check;)
endif

library: library-recursive $(LIB_OBJS)
	$(CC) -fPIC -shared $(LIB_OBJS) $(LIBS2) -o libflint.so

library-recursive:
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir) library;) 

.PHONY: profile library library-recursive recursive clean check check-recursive all

%.lo: %.c
	$(CC) -fPIC $(CFLAGS) $(INCS) -c $< -o $@

%.o: %.c
	$(CC) -fPIC $(CFLAGS) $(INCS) -c $< -o $@

BUILD_DIRS = ulong_extras fmpz fmpz_vec fmpz_poly fmpq_poly fmpz_mpoly \
   fmpz_mat mpfr_vec mpfr_mat mpfr_poly LLL nmod_vec nmod_poly nmod_mpoly \
   arith


