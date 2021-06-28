CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

BLANK :=
SLASH = \$(BLANK)

default: out.txt EigenEnergies.png

out.txt: main
	./$< > $@

EigenEnergies.png: EigenEnergies.txt
	echo '\
	set term png;\
	set output "$@";\
	set key outside;\
	set fontsize 0.5;\
	set xlabel "x";\
	set ylabel "State";\
	set title "Solutions for a particle in a box";\
	plot \
	"$<" using 1:2 with line linewidth 3 col blue title "n = 0"\
	,"$<" using 1:3 with line linewidth 3 col red title "n = 1"\
	,"$<" using 1:4 with line linewidth 3 col green title "n = 2"\
	,"$<" using 1:2 with line linetype 2 linewidth 3 col orange title "Analytic: n = 0"\
	,"$<" using 1:3 with line linetype 2 linewidth 3 col brown title "Analytic: n = 1"\
	,"$<" using 1:4 with line linetype 2 linewidth 3 col grey title "Analytic: n = 2"\
	' | pyxplot
	
main: main.o Jacobi.o

.PHONEY:clean
clean:
	$(RM) *txt main *png
