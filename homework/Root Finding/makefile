CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt hydrogen.png


out.txt: main
	./$< > $@

main: main.o RKODE.o GramSchmidt.o roots.o



hydrogen.png: hydrogen.txt
	echo '\
	set term png;\
	set output "$@";\
	set key outside;\
	set fontsize 1;\
	set xlabel "r in units of r0";\
	set ylabel "f(r)";\
	set title "Solution to the Hydrogen atom";\
	plot \
	"$<" using 1:3 with lines linewidth 2 color blue title "ODE solution"\
	,"$<" using 1:2 with lines linewidth 2 linetype 2 color red title "Exact solution"\
	' | pyxplot



.PHONEY:clean
clean:
	$(RM) *txt main *png
