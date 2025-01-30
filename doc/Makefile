# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    = -j=auto
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = flint
SOURCEDIR     = source
BUILDDIR      = build

M4 := m4

SPHINX_ACTIONS := html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf latexpdfja text man texinfo info gettext changes xml pseudoxml linkcheck doctest coverage clean

MACROS := $(SOURCEDIR)/macros.m4
SRCS := $(wildcard $(SOURCEDIR)/*.rst.pp)
OBJS := $(patsubst %.pp,%,$(SRCS))

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

preprocess: $(OBJS)

%.rst: %.rst.pp $(MACROS)
	@$(M4) -I $(SOURCEDIR) $< > $@

$(SPHINX_ACTIONS): Makefile preprocess
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile preprocess
