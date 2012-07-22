SOURCES = $(wildcard *.c)

OBJS = $(patsubst %.c, $(BUILD_DIR)/$(MOD_DIR)_%.o, $(SOURCES))

LIB_OBJS = $(patsubst %.c, $(BUILD_DIR)/%.lo, $(SOURCES))

TEST_SOURCES = $(wildcard test/*.c)

PROF_SOURCES = $(wildcard profile/*.c)

TUNE_SOURCES = $(wildcard tune/*.c)

TESTS = $(patsubst %.c, $(BUILD_DIR)/%, $(TEST_SOURCES))

TESTS_RUN = $(foreach file, $(TESTS), $(file)_RUN)

PROFS = $(patsubst %.c, %, $(PROF_SOURCES))

TUNE = $(patsubst %.c, %, $(TUNE_SOURCES))

.SECONDARY:

all: $(OBJS)

library: $(LIB_OBJS)

profile: $(PROF_SOURCES)
	$(foreach prog, $(PROFS), $(CC) -O2 -std=c99 $(INCS) $(prog).c ../profiler.o -o $(BUILD_DIR)/$(prog) $(LIBS) || exit $$?;)
        
tune: $(TUNE_SOURCES)
	$(foreach prog, $(TUNE), $(CC) -O2 -std=c99 $(INCS) $(prog).c -o $(BUILD_DIR)/$(prog) $(LIBS) || exit $$?;)

$(BUILD_DIR)/$(MOD_DIR)_%.o: %.c
	$(CC) $(CFLAGS) -c $(INCS) $< -o $@

$(BUILD_DIR)/%.lo: %.c
	$(CC) -fPIC $(CFLAGS) $(INCS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)	

check: $(TESTS) $(TESTS_RUN)

%_RUN: %
	@$<

$(BUILD_DIR)/test/%: test/%.c
	$(CC) $(CFLAGS) $(INCS) $< -o $@ $(LIBS)

.PHONY: profile tune clean check all %_RUN
