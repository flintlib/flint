SOURCES = $(wildcard *.c)

OBJS = $(patsubst %.c, $(BUILD_DIR)/$(MOD_DIR)_%.o, $(SOURCES))

LIB_OBJS = $(patsubst %.c, $(BUILD_DIR)/%.lo, $(SOURCES))

TEST_SOURCES = $(wildcard test/*.c)

PROF_SOURCES = $(wildcard profile/*.c)

TUNE_SOURCES = $(wildcard tune/*.c)

TESTS = $(patsubst %.c, $(BUILD_DIR)/%, $(TEST_SOURCES))

PROFS = $(patsubst %.c, %, $(PROF_SOURCES))

TUNE = $(patsubst %.c, %, $(TUNE_SOURCES))

all: $(OBJS)

library: $(LIB_OBJS)

profile:
	$(foreach prog, $(PROFS), $(CC) -O2 -std=c99 $(INCS) $(prog).c ../profiler.o -o $(BUILD_DIR)/$(prog) $(LIBS);)
        
tune: $(TUNE_SOURCES)
	$(foreach prog, $(TUNE), $(CC) -O2 -std=c99 $(INCS) $(prog).c -o $(BUILD_DIR)/$(prog) $(LIBS);)

$(BUILD_DIR)/$(MOD_DIR)_%.o: %.c
	$(CC) $(CFLAGS) -c $(INCS) $< -o $@

$(BUILD_DIR)/%.lo: %.c
	$(CC) -fPIC $(CFLAGS) $(INCS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)	

check: library $(TESTS)
	$(foreach prog, $(TESTS), $(prog);)

$(BUILD_DIR)/test/%: test/%.c
	$(CC) $(CFLAGS) $(INCS) $< -o $@ $(LIBS)

.PHONY: profile clean check all
