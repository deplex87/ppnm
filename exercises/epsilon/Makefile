CC = gcc
CFLAGS = -Wall -fwrapv

.PHONEY: default
default: out.txt
	cat $<

out.txt: main
	./$< > $@

.PHONEY: clean
clean:
	$(RM) main out.txt
