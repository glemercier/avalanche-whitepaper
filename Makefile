
LATEX=TZ='UTC' pdflatex

BIBTEX=bibtex
BIBFILE=paper
PAPER=paper
TEXFILES = $(wildcard *.tex) $(wildcard sections/*.tex) $(wildcard figures/*.tex) $(wildcard figures/*.pdf)
TEXFILES = $(wildcard *_fr.tex) $(wildcard sections/*_fr.tex) $(wildcard figures/*_fr.tex) $(wildcard figures/*.pdf)

SUBDIRS =
#figs

all: subdirs $(PAPER).pdf $(PAPER)_fr.pdf

subdirs:
	@for dir in $(SUBDIRS); do $(MAKE) -C $$dir; done;

$(PAPER).pdf: $(TEXFILES) $(BIBFILE).bib usenix2019_v3.sty
	$(LATEX) $(PAPER)
	$(BIBTEX) $(BIBFILE)
	$(LATEX) $(PAPER)
	$(LATEX) $(PAPER)

$(PAPER)_fr.pdf: $(TEXFILES_FR) $(BIBFILE).bib usenix2019_v3.sty
	$(LATEX) $(PAPER)_fr
	$(BIBTEX) $(BIBFILE)
	$(LATEX) $(PAPER)_fr
	$(LATEX) $(PAPER)_fr


cleantex:
	rm -f *.aux *.log *~ *.out .DS_Store *.dvi $(PAPER).ps $(PAPER).pdf $(PAPER)_fr.pdf \
		*.lot *.lof *.toc *.blg *.bbl *.ent *.bak *.glo *.gls *.ist

clean: cleantex
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir clean; done;\

