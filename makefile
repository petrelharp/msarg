%.html : %.md
	pandoc $< > $@

%.html : %.Rmd
	R --vanilla -e 'library(rmarkdown);render("$<",output_file="$@")'

