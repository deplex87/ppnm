CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

out.txt:main
	./$< > result.txt

main: main.o GramSchmidt.o

.PHONEY:clean
clean:
	$(RM) *txt main *.o
