actividades_freefem.html: actividades_freefem.md
	pandoc --toc -s -o $@ $<
