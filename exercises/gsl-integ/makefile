CFLAGS = -Wall -O1 -std=gnu11 -g
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt makefile


out.txt: main
	./$< > $@

main: main.c


.PHONEY:clean
clean:
	$(RM) *txt main
