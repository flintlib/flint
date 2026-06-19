# Convenience wrapper for the Meson build.
#
# Override BUILD to use a different build directory, for example:
#   make check BUILD=build-debug

BUILD ?= build
MESON ?= $(if $(wildcard .venv/bin/meson),.venv/bin/meson,meson)
MESON_SETUP_ARGS ?=
MESON_COMPILE_ARGS ?=
MESON_TEST_ARGS ?=
GDB ?= gdb
VALGRIND ?= valgrind
VALGRIND_ARGS ?= --track-origins=yes --leak-check=full --show-reachable=yes
DESTDIR ?=

.PHONY: all setup test tests check debug valgrind examples checkexamples profile tune coverage coverage-html install uninstall clean distclean

all: setup
	$(MESON) compile -C $(BUILD) $(MESON_COMPILE_ARGS)

setup:
	@if [ ! -f "$(BUILD)/build.ninja" ]; then \
		$(MESON) setup "$(BUILD)" $(MESON_SETUP_ARGS); \
	fi

test: check

tests: setup
	$(MESON) compile -C $(BUILD) tests $(MESON_COMPILE_ARGS)

check: setup
	$(MESON) test -C $(BUILD) $(MESON_TEST_ARGS) $(MOD)

debug: setup
ifndef MOD
	$(error Use make debug MOD=module ARGS='test arguments')
endif
ifndef ARGS
	$(error Use make debug MOD=module ARGS='test arguments')
endif
	$(MESON) test -C $(BUILD) --gdb --gdb-path $(GDB) --interactive $(MESON_TEST_ARGS) --test-args '$(ARGS)' $(MOD)

valgrind: setup
	$(MESON) test -C $(BUILD) --wrapper '$(VALGRIND) $(VALGRIND_ARGS)' $(MESON_TEST_ARGS) $(MOD)

examples: setup
	$(MESON) compile -C $(BUILD) examples $(MESON_COMPILE_ARGS)

checkexamples: setup
	$(MESON) compile -C $(BUILD) examples $(MESON_COMPILE_ARGS)
	@for src in examples/*.c; do \
		example=$$(basename "$$src" .c); \
		./dev/check_examples.sh "$$example" "$(BUILD)/examples"; \
	done

profile: setup
	$(MESON) compile -C $(BUILD) profile $(MESON_COMPILE_ARGS)

tune: setup
	$(MESON) compile -C $(BUILD) tune $(MESON_COMPILE_ARGS)

coverage:
	$(MESON) setup "$(BUILD)" --reconfigure -Db_coverage=true $(MESON_SETUP_ARGS); \
	$(MESON) test -C $(BUILD) $(MESON_TEST_ARGS) $(MOD)
	$(MESON) coverage -C $(BUILD)

coverage-html:
	$(MESON) setup "$(BUILD)" --reconfigure -Db_coverage=true $(MESON_SETUP_ARGS); \
	$(MESON) test -C $(BUILD) $(MESON_TEST_ARGS) $(MOD)
	$(MESON) coverage -C $(BUILD) --html

install: setup
	$(MESON) install -C $(BUILD) $(if $(DESTDIR),--destdir $(DESTDIR),)

uninstall:
	$(MESON) compile -C $(BUILD) uninstall

clean:
	@if [ -f "$(BUILD)/build.ninja" ]; then \
		$(MESON) compile -C $(BUILD) --clean; \
	fi

distclean:
	rm -rf "$(BUILD)"
