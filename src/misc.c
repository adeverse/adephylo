#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <R_ext/Print.h>

/* interpolate values in a density */
/* x contais n values for which we want y=f(x) */
void predict_density(double *densx, double *densy, int *densn, double *x, double *y, int *n){
	int i, idx;

	for(i=0;i<*n;i++){
		idx=0;
		while(idx < *densn && x[i]>densx[idx]){
			idx++;
		}
		if(idx==0){
			y[i] = densy[idx]/2.0;
			/* printf("\nx: %.5f, idx: %d, x.chosen: %.5f, %.5f = %.5f / 2.0", x[i], idx, densx[idx], y[i], densy[idx]); */
		} else if(idx==*densn){
			y[i] = densy[idx-1]/2.0;
			/* printf("\nx: %.5f, idx: %d, x.chosen: %.5f, %.5f = %.5f / 2.0", x[i], idx, densx[idx-1], y[i], densy[idx-1]); */
		} else {
			y[i] = (densy[idx-1] + densy[idx])/2.0;
			/* printf("\nx: %.5f, idx: %d, x.before: %.5f, x.after: %.5f, %.5f = %.5f + %.5f / 2.0", x[i], idx, densx[idx-1], densx[idx], y[i], densy[idx-1], densy[idx]); */
		}
	}
}
