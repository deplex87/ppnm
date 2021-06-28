CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt higgs_fit.png

higgs_fit.png: higgs.txt higgs_fit.txt makefile
	echo '\
	set term png;\
	set output "$@";\
	set key outside;\
	set fontsize 1;\
	set xlabel "Energy [GeV]";\
	set ylabel "Cross Section";\
	set title "Fit to BreitWigner function";\
	plot \
	"higgs_fit.txt" using 1:2 with lines color red title "Fit result" \
	,"higgs.txt" using 1:2:3 with yerrorbars color black title "Data points from CERN" \
	' | pyxplot

out.txt: main
	./$< > $@

main: main.o qnewton.o

.PHONEY:clean
clean:
	$(RM) out.txt higgs_fit.txt main *png
