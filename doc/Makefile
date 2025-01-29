# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    = -j=auto
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = flint
SOURCEDIR     = source
BUILDDIR      = build

M4 := m4

SRCS := $(wildcard $(SOURCEDIR)/*.rst.pp)
OBJS := $(patsubst %.pp,%,$(SRCS))

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

preprocess: $(OBJS)

%.rst: %.rst.pp
	@$(M4) -I $(SOURCEDIR) $< > $@

.PHONY: help Makefile preprocess

# # Catch-all target: route all unknown targets to Sphinx using the new
# # "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
# %: Makefile preprocess
# 	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
