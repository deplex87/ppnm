CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt fit_graph.png data.txt

fit_graph.png: data.txt makefile
	echo '\
	set term png;\
	set output "$@";\
	set key outside;\
	set fontsize 1;\
	set xlabel "time t (days)";\
	set ylabel "log(activity)";\
	set title "Least square fit";\
	plot \
	"$<" index 0 using 1:3:5 with errorbars title "Data Points" \
	,"$<" index 1 using 1:3 with lines linetype 2 color red title "Least square fit + sigma"\
	,"$<" index 1 using 1:4 with lines linetype 2 color red title "Least square fit - sigma"\
	,"$<" index 1 using 1:2 with lines color black title "Least square fit"\
	' | pyxplot

out.txt:main
	./$< > out.txt

main: main.o GramSchmidt.o LeastSquares.o

.PHONEY:clean
clean:
	$(RM) *txt main *png *.o
