logbookentries := $(shell grep -v '^%' logbook.tex | grep -o 'entries/[A-Za-z0-9_\-]*.[a-z]*' )
entries_inputs = $(shell grep -ho 'input/[A-Za-z0-9_\-]*.[a-z]*' entries/202*.tex)

all: logbook.pdf

clean:
	rm *.aux *.log *.toc *.out

logbook.pdf: logbook.tex $(logbookentries) $(entries_inputs)
	pdflatex -draftmode $<
	bibtex $(basename $<).aux
	pdflatex -draftmode $<
	pdflatex $<
	rm $(basename $<).log $(basename $<).out $(basename $<).toc $(basename $<).aux
	rm $(basename $<).bbl $(basename $<).blg

input:
	mkdir $@

../tasks/%: #Generic recipe to produce outputs from upstream tasks
	$(MAKE) -C $(subst output/,code/,$(dir $@)) ../output/$(notdir $@)
