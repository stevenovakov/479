FIGURES := optics cavity eom_wave rf_design
FIGURES_PDF := $(foreach figure, $(FIGURES), figures/$(figure).pdf)

.PHONY: tex
tex: final.tex $(FIGURES_PDF)
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<

%.pdf: %.svg
	rsvg-convert -f pdf -o $@ $^
