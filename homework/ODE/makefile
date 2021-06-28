CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt DifferentialEquation.png SIRmodel.png

out.txt: main
	./$< > $@

main: main.c RKODE.c


DifferentialEquation.png: DifferentialEquation.gpi
	pyxplot $<
	
SIRmodel.png: SIRmodel.gpi
	pyxplot $<


.PHONEY:clean
clean:
	$(RM) *txt main *png
