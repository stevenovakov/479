.PHONY: tex
tex: final.tex
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<

%.pdf: %.svg
	rsvg-convert -f pdf -o $@ $^
