CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


out.txt: main
	./$< >$@
	cat out.txt


.PHONEY:clean
clean:
	$(RM) main out.txt 
