SOURCES = $(wildcard *.c)

OBJS = $(patsubst %.c, %.o, $(SOURCES))

LIB_OBJS = $(patsubst %.c, %.lo, $(SOURCES))

TEST_SOURCES = $(wildcard test/*.c)

PROF_SOURCES = $(wildcard profile/*.c)

TESTS = $(patsubst %.c, %, $(TEST_SOURCES))

PROFS = $(patsubst %.c, %, $(PROF_SOURCES))

all: $(OBJS)

library: $(LIB_OBJS)

profile:
	$(foreach prog, $(PROFS), $(CC) -O2 -std=c99 $(INCS) $(prog).c ../profiler.o -o $(BUILD_DIR)/$(prog) $(LIBS) -lflint;)
        
%.o: %.c
	$(CC) $(CFLAGS) -c $(INCS) $< -o $(BUILD_DIR)/$@

%.lo: %.c
	$(CC) -fPIC $(CFLAGS) $(INCS) -c $< -o $(BUILD_DIR)/$@

clean:
	rm -rf $(BUILD_DIR)	

check: library
	$(foreach prog, $(TESTS), $(CC) $(CFLAGS) $(INCS) $(prog).c -o $(BUILD_DIR)/$(prog) $(LIBS) -lflint;)
	$(foreach prog, $(TESTS), $(BUILD_DIR)/$(prog);)

.PHONY: profile clean check all
