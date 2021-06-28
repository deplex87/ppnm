CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default:out.main.txt gsl_eval.txt LinearSpline.png CubicSpline.png QuadraticSpline.png


out.main.txt: main
	./$< > $@

main: main.o Linear_Spline.o Quadratic_Spline.o Cubic_Spline.o Binary_Search.o

main.o Splines.o quadratic: Linear_Spline.h Quadratic_Spline.h Cubic_Spline.h Binary_Search.h

gsl_eval.txt: gsl_eval
	./$< > $@


LinearSpline.png: out.main.txt Makefile gsl_eval.txt
	echo '\
	set term png;\
	set output "$@";\
	set key outside;\
	set fontsize 1;\
	set linewidth 2;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Linear Spline";\
	plot [0:7.5][] \
	"$<" index 0 using 1:2 with points pointtype 3 title "Data Points sin(x)" \
	,"$<" index 1 using 1:2 with lines color blue title "Linear Spline"\
	,"$<" index 1 using 1:3 with lines color green title "Linear Spline Integral"\
	,"gsl.txt" index 1 using 1:2 with lines linetype 2 color red title "Linear Spline GSL"\
	,"gsl.txt" index 2 using 1:2 with lines linetype 2 color red title "Linear Integral GSL"\
	' | pyxplot

QuadraticSpline.png: out.main.txt Makefile gsl_eval.txt
	echo '\
	set term png;\
	set output "$@";\
	set fontsize 1;\
	set key outside;\
	set linewidth 2;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Quadratic Spline";\
	plot [0:7.5][] \
	"$<" index 0 using 1:2 with points pointtype 3 title "Data Points sin(x)" \
	,"$<" index 2 using 1:2 with lines color blue title "Quadratic Spline"\
	,"$<" index 2 using 1:3 with lines color black title "Quadratic Spline Derivative"\
	,"$<" index 2 using 1:4 with lines color green title "Quadratic Spline Integral"\
	,"$<" index 0 using 1:3 with lines linetype 2 color red title "Analytic Derivative"\
	,"$<" index 0 using 1:4 with lines linetype 2 color red title "Analytic Integral"\
	' | pyxplot


CubicSpline.png:  out.main.txt Makefile gsl_eval.txt
	echo '\
	set term png;\
	set output "$@";\
	set fontsize 1;\
	set linewidth 2;\
	set key outside;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Cubic Spline";\
	plot [0:7.5][] \
	"$<" index 0 using 1:2 with points pointtype 3 title "Data Points sin(x)" \
	,"$<" index 3 using 1:2 with lines color green title "Cubic Spline"\
	,"$<" index 3 using 1:3 with lines color blue title "Cubic Spline Derivative"\
	,"$<" index 3 using 1:4 with lines color black title "Cubic Spline Integral"\
	,"gsl.txt" index 1 using 1:3 with lines linetype 2 color red title "Cubic Spline GSL"\
	,"gsl.txt" index 3 using 1:2 with lines linetype 2 color red title "Cubic Spline GSL Derivative"\
	,"gsl.txt" index 2 using 1:3 with lines linetype 2 color red title "Cubic Spline GSL Integral"\
	' | pyxplot


.PHONEY:clean
clean:
	$(RM) *txt main gsl *png quadratic *.o
