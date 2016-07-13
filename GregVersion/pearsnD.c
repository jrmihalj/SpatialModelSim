#define TINY 1.0e-20

void pearsn(double x[], double y[], unsigned long n, double *r){

unsigned long j;
float yt, xt, t, df;
float syy=0.0, sxy=0.0, sxx=0.0, ay = 0.0, ax = 0.0;

for(j=1;j<=n;j++){
	ax += x[j];	ay += y[j];
}

ax /= n;
ay /= n;
for(j=1;j<=n;j++){
	xt = x[j]-ax;
	yt = y[j]-ay;
	sxx += xt*xt;
	syy += yt*yt;
	sxy += xt*yt;
}
*r = sxy/(sqrt(sxx*syy));

}
