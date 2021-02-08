TEXT = gpdfa_paper
REFS = gpdfa

all: $(TEXT).pdf

%.pdf: %.Rmd $(REFS).bib options.sty
	Rscript -e "rmarkdown::render('gpdfa_paper.Rmd')"
continuous:
	while true; do make --silent; sleep 1; done
