# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

#ifndef Radial_h
#define Radial_h

void legendre_compute_glr ( int n, double x[], double w[] );
void legendre_compute_glr0 ( int n, double *p, double *pp );
void legendre_compute_glr1 ( int n, double *roots, double *ders );
void legendre_compute_glr2 ( double p, int n, double *roots, double *ders );
void legendre_handle ( int n, double a, double b );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
void rescale ( double a, double b, int n, double x[], double w[] );
double rk2_leg ( double t, double tn, double x, int n );
void timestamp ( void );
double ts_mult ( double *u, double h, int n );
double wtime ( );

#endif
