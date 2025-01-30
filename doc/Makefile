# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    = -j=auto
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = flint
SOURCEDIR     = source
ORIGSRCDIR    = origsrc
BUILDDIR      = build

M4 := m4

SPHINX_ACTIONS := html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf latexpdfja text man texinfo info gettext changes xml pseudoxml linkcheck doctest coverage

CONFIG := macros.m4
SRCS := $(wildcard $(ORIGSRCDIR)/*.rst)
OBJS := $(patsubst $(ORIGSRCDIR)/%,$(SOURCEDIR)/%,$(SRCS))

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

$(SOURCEDIR)/%.rst: $(ORIGSRCDIR)/%.rst $(CONFIG)
	@$(M4) $< > $@

preprocess: $(OBJS)

$(SPHINX_ACTIONS): Makefile preprocess
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	rm -f $(OBJS)
	rm -rf $(BUILDDIR)

print-%:
	@echo "$*=$($*)"

.PHONY: help Makefile preprocess
