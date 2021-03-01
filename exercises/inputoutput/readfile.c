#include<stdio.h>
#include<math.h>
int main(){
	double x;
	FILE* out=fopen("q3.out.txt","w");
	FILE* input=fopen("inputfile.txt","r");
	while (fscanf(input, "%lg", & x ) == 1 ){
		fprintf(out,"x= %lg, sin(x)=%lg, cos(x)=%lg \n", x,sin(x),cos(x));
        }
	fclose(input);
	fclose(out);
return 0;
}
