SHELL=/bin/sh
LD_PRE=LD_LIBRARY_PATH=..

CC=gcc

WARNS=-ansi -pedantic -Wall -Wextra -Werror -Wno-unused-parameter -Wno-unused
CFLAGS=$(WARNS) -mpopcnt -g -O0
CFLAGS_PROF=$(WARNS) -mpopcnt -funroll-loops -O3
LIBS=-lflint -lmpfr -lgmp -lm

HEADERS=../*.h
SOURCES=$(wildcard ./*.c)
TEST_SOURCES=$(wildcard test/t-*.c)
PROF_SOURCES=$(wildcard profile/p-*.c)
TESTS=$(patsubst test/t-%.c,%,$(TEST_SOURCES))
PROFS=$(patsubst profile/p-%.c,%,$(PROF_SOURCES))
TEST_RUNS=$(patsubst %,t-%,$(TESTS))
PROF_RUNS=$(patsubst %,p-%,$(PROFS))

BD=build

TESTS_EXTRA_SRC=../test_helpers.c
PROFS_EXTRA_SRC=../profiler.c

TEST_EXES=$(patsubst %,$(BD)/t-%,$(TESTS))
PROF_EXES=$(patsubst %,$(BD)/p-%,$(PROFS))

all: tests profile 

$(TEST_RUNS): %: build/%
	@echo
	@echo "RUNNING TEST $@"
	$(LD_PRE) $<

tests: $(TEST_RUNS)

$(TEST_EXES): $(BD)/%: test/%.c $(SOURCES) $(HEADERS) $(TESTS_EXTRA_SRC)
	@echo
	@echo "COMPILING $@"
	$(CC) $(CFLAGS) -I.. $< $(SOURCES) $(TESTS_EXTRA_SRC) -L.. $(LIBS) -o $@

$(PROF_RUNS): %: build/%
	@echo
	@echo "RUNNING PROFILE $@"
	$(LD_PRE) $<

profile: $(PROF_RUNS)

$(PROF_EXES): $(BD)/%: profile/%.c $(SOURCES) $(HEADERS) $(PROFS_EXTRA_SRC)
	@echo
	@echo "COMPILING $@"
	$(CC) $(CFLAGS_PROF) -I.. $< $(SOURCES) $(PROFS_EXTRA_SRC) -L.. $(LIBS) -o $@

clean:
	rm -f $(BD)/*

.PHONY: clean all tests profile $(TEST_RUNS) $(PROF_RUNS)
