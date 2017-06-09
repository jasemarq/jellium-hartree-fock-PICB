#ifndef MEMORY_H
#define MEMORY_H

#include<malloc.h>
#include</usr/include/complex.h>
#include<assert.h>
#include<stdio.h>
#include<stdlib.h>

int *VEC_INT(int dim);

double *VEC_DOUBLE(int dim);

char *VEC_CHAR(int dim);

char **MAT_CHAR(int dim1, int dim2);

int **SQMAT_2_INT(int dim);

char ** SQMAT_2_CHAR(int dim);

double **SQMAT_2_DOUBLE(int dim);

double ****SQMAT_4_DOUBLE(int dim);

int ***SQMAT_3_INT(int dim);

int **MAT_INT(int dim1, int dim2);

int ***MAT_3_INT(int dim1, int dim2, int dim3);

double **MAT_DOUBLE(int dim1, int dim2);

double complex *VEC_CDOUBLE(int dim);

void FREE_4(double ****M,int dim1, int dim2, int dim3);

void FREE_4_INT(int ****M, int dim);

void FREE_4_DOUBLE(double ****M,int dim);

void FREE_2_DOUBLE(double **M,int dim);

void FREE_2_CHAR(char **M,int dim);

void FREE_2_INT(int **M,int dim);

void FREE_3_INT(int ***M,int dim);






#endif

