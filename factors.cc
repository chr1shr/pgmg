#include <cstdio>
#include <cmath>

const int n=131072;

int main() {
	int i,j,k[n];

	// Calculates the number of factors of each number from 1 to n
	for(i=1;i<n;i++) {
		k[i]=0;
		for(j=1;j*j<=i;j++) {
			if(i%j==0) k[i]++;
		}
	}

	// Print out the results, along with the number of a factors averaged
	// over a +/- 50 window
	int il,ip,s;
	for(i=1;i<n;i++) {

		// Compute the window range, taking into account boundary cases
		il=i-50;if(il<1) il=1;
		ip=i+51;if(ip>n) ip=n;

		// Calculate the average and print out the results.
		for(s=0,j=il;j<ip;j++) s+=k[j];
		printf("%d %d %g\n",i,k[i],(static_cast<double>(s))/(ip-il));
	}
}
