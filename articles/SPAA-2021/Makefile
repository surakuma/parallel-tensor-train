.PHONY: document

FILE=ptt

FIGURES=

all: pdflatex

document:
	rubber -dv $(FILE).tex

mini: document
	pdfnup --nup 2x1 $(FILE).pdf

pdflatex: $(FIGURES)
	pdflatex $(FILE).tex
	bibtex $(FILE)
	pdflatex $(FILE).tex
	pdflatex $(FILE).tex

clean:
	rubber --clean $(FILE).tex
	rm -f $(FILE).pdf

aspell:
	aspell --lang=en -t=tex  --encoding=utf-8 -c $(FILE).tex

