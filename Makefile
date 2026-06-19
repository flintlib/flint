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

.PHONY: all setup library shared static test tests check debug valgrind examples checkexamples profile tune coverage coverage_html install uninstall clean distclean print-%

ifneq ($(NJOBS),)
number_generator = $(words $2) $(if $(word $1,$2),,$(call number_generator,$1,w $2))
THREAD_LIST := $(call number_generator,$(NJOBS),w)
CHECK_THREAD_TARGETS := $(foreach mod,$(MOD),$(foreach num,$(THREAD_LIST),check-$(mod)-thread-$(num)))

define check_mod_build_rule
.PHONY: check-$(1)-build
check-$(1)-build: setup
	$$(MESON) compile -C $$(BUILD) src/$(1)/test/main $$(MESON_COMPILE_ARGS)
endef

define check_thread_rule
.PHONY: check-$(1)-thread-$(2)
check-$(1)-thread-$(2): check-$(1)-build
	@exe="$$(BUILD)/src/$(1)/test/main"; \
	if [ ! -f "$$$$exe" ]; then exe="$$$$exe.exe"; fi; \
	echo "$$$$exe --numthreads=$$(NJOBS) --thread=$(2)"; \
	LD_LIBRARY_PATH="$$(BUILD):$$$$LD_LIBRARY_PATH" DYLD_LIBRARY_PATH="$$(BUILD):$$$$DYLD_LIBRARY_PATH" "$$$$exe" --numthreads=$$(NJOBS) --thread=$(2)
endef

$(foreach mod,$(MOD),$(eval $(call check_mod_build_rule,$(mod))))
$(foreach mod,$(MOD),$(foreach num,$(THREAD_LIST),$(eval $(call check_thread_rule,$(mod),$(num)))))
endif

all: library

library: setup
	$(MESON) compile -C $(BUILD) $(MESON_COMPILE_ARGS)

shared:
	$(MESON) setup "$(BUILD)" --reconfigure -Ddefault_library=shared $(MESON_SETUP_ARGS)
	$(MESON) compile -C $(BUILD) $(MESON_COMPILE_ARGS)

static:
	$(MESON) setup "$(BUILD)" --reconfigure -Ddefault_library=static $(MESON_SETUP_ARGS)
	$(MESON) compile -C $(BUILD) $(MESON_COMPILE_ARGS)

setup:
	@if [ ! -f "$(BUILD)/build.ninja" ]; then \
		$(MESON) setup "$(BUILD)" $(MESON_SETUP_ARGS); \
	fi

test: check

tests: setup
	$(MESON) compile -C $(BUILD) tests $(MESON_COMPILE_ARGS)

check: setup
ifneq ($(NJOBS),)
ifdef MOD
ifndef ARGS
check: $(CHECK_THREAD_TARGETS)
endif
endif
endif

ifdef PYTHON
	$(MESON) compile -C $(BUILD) $(MESON_COMPILE_ARGS)
	@LD_LIBRARY_PATH="$(BUILD):$${LD_LIBRARY_PATH}" DYLD_LIBRARY_PATH="$(BUILD):$${DYLD_LIBRARY_PATH}" python3 src/python/flint_ctypes.py
	@echo ''
	@echo 'All Python tests passed.'
else
ifdef ARGS
ifneq ($(words $(sort $(MOD))),1)
	$(error Can only check one module with arguments)
endif
	$(MESON) test -C $(BUILD) $(MESON_TEST_ARGS) --test-args='$(ARGS)' $(MOD)
else
ifdef NJOBS
ifdef MOD
	@echo ''
	@echo 'All tests passed for $(sort $(MOD)).'
else
	$(MESON) test -C $(BUILD) $(MESON_TEST_ARGS)
endif
else
	$(MESON) test -C $(BUILD) $(MESON_TEST_ARGS) $(MOD)
endif
endif
endif

debug: setup
ifndef MOD
	$(error Use make debug MOD=module ARGS='test arguments')
endif
ifndef ARGS
	$(error Use make debug MOD=module ARGS='test arguments')
endif
	$(MESON) test -C $(BUILD) --gdb --gdb-path $(GDB) --interactive $(MESON_TEST_ARGS) --test-args='$(ARGS)' $(MOD)

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
ifdef MOD
	$(MESON) compile -C $(BUILD) $(MOD)-profile $(MESON_COMPILE_ARGS)
else
	$(MESON) compile -C $(BUILD) profile $(MESON_COMPILE_ARGS)
endif

tune: setup
	$(MESON) compile -C $(BUILD) tune $(MESON_COMPILE_ARGS)

coverage:
	$(MESON) setup "$(BUILD)" --reconfigure -Db_coverage=true $(MESON_SETUP_ARGS); \
	$(MESON) test -C $(BUILD) $(MESON_TEST_ARGS) $(MOD)
	$(MESON) coverage -C $(BUILD)

coverage_html:
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

print-%:
	@echo "$*=$($*)"
