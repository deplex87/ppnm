CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY:default
default: erf.pyxplot.png

erf.pyxplot.png: outputdata.txt Makefile data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "math.h erf"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "myerf.c"\
	,"data.txt" using 1:2 with points pointtype 2 title "Tabulated data"\
	' | tee erf.pyxplot.gpi | pyxplot

outputdata.txt: main
	./$< > $@

main: main.o myerf.o
main.o Erf.o: myerf.h


.PHONEY:clean
clean:
	$(RM) *.o main *.png outputdata.txt tmp.gpi
