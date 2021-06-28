CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY:default
default: gamma.pyxplot.png

gamma.pyxplot.png: outputdata.txt makefile tabulateddata.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma function";\
	plot \
	 "$<" using 1:2 with line title "math.h gamma"\
	,"$<" using 1:3 with line title "gsl-sf-gamma"\
	,"$<" using 1:4 with line title "mygamma.c"\
	,"tabulateddata.txt" using 1:2 with points pointtype 2 title "Tabulated data"\
	' | tee gamma.pyxplot.gpi | pyxplot

outputdata.txt: main
	./$< > $@

main: main.o mygamma.o
main.o Erf.o: mygamma.h


.PHONEY:clean
clean:
	$(RM) *.o main *.png outputdata.txt gamma.pyxplot.gpi
