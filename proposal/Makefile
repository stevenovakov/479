.PHONY: tex
tex: proposal479.tex setup.pdf spectrum.pdf
	pdflatex $<
	bibtex $<
	pdflatex $<
	pdflatex $<

%.pdf: %.svg
	rsvg-convert -f pdf -o $@ $^
