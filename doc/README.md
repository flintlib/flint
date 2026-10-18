# FLINT documentation

The documentation sources are in `src/`, written in the Markdown dialect
described in `fdoc/README.md`; `fdoc/` is the program that builds them.

From the top-level directory, after `./configure`:

    make doc          # HTML in build/doc/html (needs Python 3)
    make docpdf       # PDF in build/doc/latex/flint.pdf (also needs pdflatex)
    make checkdoc     # check against the headers; compile and run all examples

The online documentation is at https://flintlib.org/doc/.
