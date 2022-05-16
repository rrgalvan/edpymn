all: actividades_freefem.html actividades_fenics.html

actividades_freefem.html: actividades_freefem.md
	pandoc --toc -s -o $@ $<

actividades_fenics.html: actividades_fenics.md
	pandoc --toc -s -o $@ $<
