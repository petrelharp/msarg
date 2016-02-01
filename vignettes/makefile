htmls = $(patsubst %.Rmd,%.html,$(wildcard *.Rmd))

all : $(htmls)

%.html : %.md
	pandoc $< > $@

%.html : %.Rmd
	R --vanilla -e 'library(rmarkdown);render("$<",output_file="$@")'

