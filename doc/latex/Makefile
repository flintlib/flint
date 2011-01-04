CC    ?= gcc
CFLAGS = -ansi -pedantic

doc :
	$(CC) -c $(CFLAGS) doc.c -o doc.o
	$(CC) doc.o -o doc
	./doc
	pdflatex doc.tex
	bibtex doc
	pdflatex doc.tex
	pdflatex doc.tex
	mv doc.tex doc.tex.bak
	mv title.tex title.tex.bak
	rm *.tex
	mv doc.tex.bak doc.tex
	mv title.tex.bak title.tex
	rm doc.o doc
	rm doc.aux doc.log doc.out doc.toc doc.bbl doc.blg

