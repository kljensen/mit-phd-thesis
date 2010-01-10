#  Makefile for LaTeX documents
#
TARGET = main.pdf
TOUCH=touch
PS = dvips
PDFLATEX = pdflatex
LATEX = $(PDFLATEX)
BIBTEX = bibtex
PS2PDF= ps2pdf
DELETE =  *.aux *.log *.ps *.dvi *.bbl *.blg *~ $(TARGET) *.brf *.idx *.ilg *.ind *.lof *.lot *.out *.toc
DVIPSFLAGS= -Pcmz -Pamz -Ppdf -G0 -tletter
TOUCH=touch
MAKE=make
MAKEIDX=makeindex

.SUFFIXES: .pdf .tex .ps

all: $(TARGET)
	$(MAKE) $(TARGET)

$(TARGET): 

.ps.pdf:
	$(PS2PDF) $<
	
.dvi.ps:
	$(PS) $(DVIPSFLAGS) -t letter -o $@ $<

.tex.pdf: 
	$(LATEX) $<
	$(LATEX) $<
	$(BIBTEX) $*
#	$(MAKEIDX) $*
	$(LATEX) $<
	$(BIBTEX) $*
	$(LATEX) $<

clean: 
	rm -f $(DELETE)
