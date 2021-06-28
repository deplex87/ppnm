#include<assert.h>

int binsearch(int n, double* x, double z){
	assert(n>1 && x[0]<=z && z<=x[n-1]); // Tjekker om Z faktisk er i intervallet og om vi har mere end 1 datapunkt, hvis ikke giver den fejl.
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
	} 
