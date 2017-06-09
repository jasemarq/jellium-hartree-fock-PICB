#include<stdio.h>
#include<cmath>
#include<math.h>
#include<stdlib.h>
#include<complex.h>
#include<time.h>
#include"gauss_legendre.h"

double pi=4*atan(1.0);

double g_pq(double p, double q, double r);
double pq_int(double px, double py, double pz, double qx, double qy, double qz);


int main() {
  int i, n;
  double a, b;
  double *x, *w;
  n = 320;
  x = (double *)malloc(n*sizeof(double));
  w = (double *)malloc(n*sizeof(double));


  n = 320;
  a = 0.;
  b = 1.;

  legendre_compute_glr(n, x, w);
  //legendre_compute_glr ( int n, double x[], double w[] )

  rescale ( a, b, n, x, w );
  for (i=0; i<n; i++) {

    printf("  %12.10f  %12.10f \n",x[i],w[i]);

  }


  printf("  (1|2) is %f\n",pq_int(0,0,2,0,0,2));

return 0;

}


double pq_int(double px, double py, double pz, double qx, double qy, double qz) {

  double dr = 0.005;
  int max = 200;
  double sum = 0.;
  double num, denom;
  double x, y, z;

  for (int i=0; i<max; i++) {
    x = dr*i + 0.005;
    for (int j=0; j<max; j++) {
      y = dr*j + 0.005;
      for (int k=0; k<max; k++) {
        z = dr*k + 0.005;
        num = g_pq(px, qx, x)*g_pq(py, qy, y)*g_pq(pz, qz, z);
        denom = sqrt(x*x+y*y+z*z);
        sum += (num/denom)*dr*dr*dr;
        printf("  sum %f  x %f  y %f  z %f\n",sum, x, y, z);
      }
    }
  }       

  return (8./pi)*sum;

}


double g_pq(double p, double q, double x) {

  int d;
  d = (int)(fabs(p-q));
  double g;
  g = 0.;
  if (p == q && p == 0) {

    g = 1 - x;

  }
  else if ( p == q && p > 0 ) {

    g = (1 - x)*cos(p*pi*x)/2. - sin(p*pi*x)/(2*p*pi);

  }
  else if ( (d % 2)==0) {

    g = (q*sin(q*pi*x) - p*sin(p*pi*x))/((p*p-q*q)*pi);
  } 
  else g = 0.;

  return g;
}
