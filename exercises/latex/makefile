CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

exp.pdf: ExponentialFunction.txt ExponentialPyxplot.pdf
	pdflatex $<
	pdflatex $<

out.txt: main
	./$< > out.txt

ExponentialPyxplot.pdf: out.txt makefile
	echo '\
	set term pdf;\
	set width 6;\
	set output "$@";\
	set key outside;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Exponential function";\
	plot \
	 "$<" using 1:2 with line title "Exp from math.h"\
	,"$<" using 1:3 with points pointtype 11 color blue title "Function ex"\
	' | pyxplot

.PHONEY: clean
clean:
	$(RM) ExponentialFunction.pdf ExponentialFunction.log ExponentialFunction.aux main out.txt ExponentialFunction.tex ExponentialPyxplot.pdf
