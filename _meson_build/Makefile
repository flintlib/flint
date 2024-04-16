BUILDDIR = build.meson

all: setup
	meson compile -C $(BUILDDIR)

setup: $(BUILDDIR)/.build-dir-created

$(BUILDDIR)/.build-dir-created:
	meson setup $(BUILDDIR)
	touch $(BUILDDIR)/.build-dir-created

clean:
	rm -rf $(BUILDDIR)

install: compile
	meson install -C $(BUILDDIR)

test: setup
	meson configure $(BUILDDIR) -Dtests=enabled
	meson test -C $(BUILDDIR)

coverage: setup
	# Coverage tool e.g. gcovr needs to be installed *before* meson configure:
	#    pip install gcovr
	meson configure $(BUILDDIR) -Dtests=enabled -Db_coverage=true
	meson test -C $(BUILDDIR)
	cd $(BUILDDIR) && ninja coverage-html
	# open $(BUILDDIR)/meson-logs/coveragereport/index.html
