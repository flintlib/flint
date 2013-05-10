SOURCES = $(wildcard *.c)

OBJS = $(patsubst %.c, $(BUILD_DIR)/$(MOD_DIR)_%.o, $(SOURCES))

LOBJS = $(patsubst %.c, $(BUILD_DIR)/%.lo, $(SOURCES))
MOD_LOBJ = $(BUILD_DIR)/../$(MOD_DIR).lo

TEST_SOURCES = $(wildcard test/*.c)

PROF_SOURCES = $(wildcard profile/*.c)

TUNE_SOURCES = $(wildcard tune/*.c)

TESTS = $(patsubst %.c, $(BUILD_DIR)/%, $(TEST_SOURCES))

TESTS_RUN = $(patsubst %, %_RUN, $(TESTS))

PROFS = $(patsubst %.c, %, $(PROF_SOURCES))

TUNE = $(patsubst %.c, %, $(TUNE_SOURCES))

all: shared static

shared: $(MOD_LOBJ)

static: $(OBJS)

profile: $(PROF_SOURCES)
	$(foreach prog, $(PROFS), $(CC) $(ABI_FLAG) -std=c99 -O2 -g $(INCS) $(prog).c ../profiler.o -o $(BUILD_DIR)/$(prog) $(LIBS) || exit $$?;)

tune: $(TUNE_SOURCES)
	$(foreach prog, $(TUNE), $(CC) $(CFLAGS) $(INCS) $(prog).c -o $(BUILD_DIR)/$(prog) $(LIBS) || exit $$?;)

$(BUILD_DIR)/$(MOD_DIR)_%.o: %.c
	$(CC) $(CFLAGS) -c $(INCS) $< -o $@

$(MOD_LOBJ): $(LOBJS)
	$(CC) $(ABI_FLAG) -Wl,-r $^ -o $@ -nostdlib

$(BUILD_DIR)/%.lo: %.c
	$(CC) $(PIC_FLAG) $(CFLAGS) $(INCS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(MOD_LOBJ)

check: $(TESTS) $(TESTS_RUN)

$(BUILD_DIR)/test/%: test/%.c
	$(CC) $(CFLAGS) $(INCS) $< ../test_helpers.o -o $@ $(LIBS)

%_RUN: %
	@$<

.PHONY: profile tune clean check all shared static %_RUN
