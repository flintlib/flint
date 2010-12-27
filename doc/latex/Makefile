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

clean :
	rm doc.o doc
	rm doc.aux doc.log doc.out doc.toc

