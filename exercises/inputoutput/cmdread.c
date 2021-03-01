#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main( int argc, char** argv){

	if ( argc < 2){
		fprintf(stderr, "No arguments were given.\n");
	}
	else {
		for (int i = 1; i < argc; i++){
			double argX = atof(argv[i]);
			fprintf(stdout, "x= %lg, sin(x)=%lg, cos(x)=%lg \n", argX, sin(argX), cos(argX));
		}
	}

	return 0;
}
