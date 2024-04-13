BUILDDIR = build2

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
