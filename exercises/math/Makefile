CC       =  gcc
CFLAGS   =  -Wall
LDFLAGS  =  -lm

default : out.txt
	cat out.txt

out.txt : math                 
	./$< > $@       

math : math.o       	       		
	$(CC) $< -o $@ $(CFLAGS) $(LDFLAGS)    
	
math.o : math.c      	       
	$(CC) -c $<          

.PHONEY: clean
clean:                 	       
	rm -f out.txt math math.o
