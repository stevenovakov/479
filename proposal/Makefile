.PHONY: tex
tex: proposal479.tex setup.pdf spectrum.pdf doppler.pdf
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<

%.pdf: %.svg
	rsvg-convert -f pdf -o $@ $^
