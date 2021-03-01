CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm

.PHONEY:default

default: q1.out.txt q2.out.txt q3.out.txt
	cat  opg*

q1.out.txt:cmdread
	./cmdread 1 2 3 4 5 > $@

q2.out.txt:stdinput
	echo 5 6 7 8 9 10 | ./stdinput > $@

q3.out.txt:readfile
	./readfile input.txt > $@

.PHONEY:clean
clean:
	$(RM) opg* input* out* stdinput cmdread
