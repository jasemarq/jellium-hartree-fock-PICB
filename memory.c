#include"memory.h"

/**************************************************************

        This header contains functions to allocate 
            memory for matrices and tensors.

**************************************************************/
int *VEC_INT(int dim){
  int *v,i;
  v = (int *)malloc(dim*sizeof(int));
  if (v==NULL) {
     printf("\n\nVEC_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0;
  return v;
}
double *VEC_DOUBLE(int dim){
  int i;
  double *v;
  v = (double *)malloc(dim*sizeof(double));
  if (v==NULL) {
     printf("\n\nVEC_DOUBLE: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0.0;
  return v;
}

double complex *VEC_CDOUBLE(int dim) {
  int i;
  double complex *v;
  v = (double complex *)malloc(dim*sizeof(double complex));
  if (v==NULL) {
     printf("\n\nVEC_CDOUBLE:  Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0.0 + I*0.;
  return v;
}

char *VEC_CHAR(int dim){
  char *v;
  v = (char *)malloc(dim*sizeof(char));
  if (v==NULL) {
     printf("\n\nVEC_CHAR: Memory allocation error\n\n");
     exit(0);
  }
  return v;
}
char **MAT_CHAR(int dim1, int dim2){
  int i,j;
  char **C;
  C = (char**)malloc(dim1*sizeof(char*));
  if (C==NULL) {
     printf("\n\nMAT_CHAR: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim1; i++){
      C[i] = (char*)malloc(dim2*sizeof(char));
      if (C[i]==NULL) {
         printf("\n\nMAT_CHAR: Memory allocation error\n\n");
         exit(0);
      }
      for (j=0; j<dim2; j++)
          C[i][j] = 'n';
  }
return C;
}
int **SQMAT_2_INT(int dim){
  int i,j;
  int **M;
  M = (int **)malloc(dim*sizeof(int *));
  if (M==NULL) {
     printf("\n\nSQMAT_2_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++){
      M[i] = (int *)malloc(dim*sizeof(int));
      if (M[i]==NULL) {
         printf("\n\nSQMAT_2_INT: Memory allocation error\n\n");
         exit(0);
      }
      for (j=0; j<dim; j++)
          M[i][j] = 0;
  }
  return M;
}

char **SQMAT_2_CHAR(int dim){
  int i,j;
  char **M;
  M = (char **)malloc(dim*sizeof(char *));
  if (M==NULL) {
     printf("\n\nSQMAT_2_CHAR: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++){
      M[i] = (char *)malloc(dim*sizeof(char));
      if (M[i]==NULL) {
         printf("\n\nSQMAT_2_CHAR: Memory allocation error\n\n");
         exit(0);
      }
      for (j=0; j<dim; j++)
          M[i][j] = 'n';
  }
  return M;
}
double **SQMAT_2_DOUBLE(int dim){
  int i,j;
  double **M;
  M = (double **)malloc(dim*sizeof(double *));
  if (M==NULL){
     printf("\n\nSQMAT_2_DOUBLE: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++){
      M[i] = (double *)malloc(dim*sizeof(double));
      if (M[i]==NULL) {
         printf("\n\nSQMAT_2_DOUBLE: Memory allocation error\n\n");
         exit(0);
      }
      for (j=0; j<dim; j++)
          M[i][j] = 0.0;
  }
  return M;
}
double ****SQMAT_4_DOUBLE(int dim){
  int i,j,k,l;
  double ****M;
  M = (double ****)malloc(dim*sizeof(double ***));
  if (M==NULL) {
     printf("\n\nSQMAT_4_DOUBLE: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++){
      M[i] = (double ***)malloc(dim*sizeof(double **));
      if (M[i]==NULL) {
         printf("\n\nSQMAT_4_DOUBLE: Memory allocation error\n\n");
         exit(0);
      }
      for (j=0; j<dim; j++){
          M[i][j] = (double **)malloc(dim*sizeof(double *));
          if (M[i][j]==NULL) {
             printf("\n\nSQMAT_4_DOUBLE: Memory allocation error\n\n");
             exit(0);
          }
          for (k=0; k<dim; k++){
              M[i][j][k] = (double *)malloc(dim*sizeof(double));
              if (M[i][j][k]==NULL) {
                 printf("\n\nSQMAT_4_DOUBLE: Memory allocation error\n\n");
                 exit(0);
              }
              for (l=0; l<dim; l++)
                  M[i][j][k][l] = 0.0;
          }
      }
  }
  return M;
}
int ***SQMAT_3_INT(int dim){
  int i,j,k;
  int ***M;
  M = (int ***)malloc(dim*sizeof(int **));
  if (M==NULL) {
     printf("\n\nSQMAT_3_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++){
      M[i] = (int **)malloc(dim*sizeof(int *));
      if (M[i]==NULL) {
         printf("\n\nSQMAT_3_INT: Memory allocation error\n\n");
         exit(0);
      }
      for (j=0; j<dim; j++){
          M[i][j] = (int *)malloc(dim*sizeof(int));
          if (M[i][j]==NULL) {
             printf("\n\nSQMAT_3_INT: Memory allocation error\n\n");
             exit(0);
          }
          for (k=0; k<dim; k++)
              M[i][j][k] = 0;
      }
  }
  return M;
}
int **MAT_INT(int dim1, int dim2){
  int i,j,k;
  double sum=0.0;
  int **M;
  M = (int **)malloc(dim1*sizeof(int *));
  if (M==NULL) {
     printf("\n\nMAT_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim1; i++){
      M[i] = (int *)malloc(dim2*sizeof(int));
      if (M[i]==NULL) {
         printf("\n\nMAT_INT: Memory allocation error\n\n");
         exit(0);
      }
      for (j=0; j<dim2; j++){
          M[i][j] = 0;
      }
  }
  return M;
}
int ***MAT_3_INT(int dim1, int dim2, int dim3){
  int i,j,k;
  int ***M;
  M = (int ***)malloc(dim1*sizeof(int **));
  if (M==NULL) {
     printf("\n\nMAT_3_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim1; i++){
      M[i] = (int **)malloc(dim2*sizeof(int *));
      if (M[i]==NULL) {
         printf("\n\nMAT_3_INT: Memory allocation error\n\n");
         exit(0);
      }
      for (j=0; j<dim2; j++){
          M[i][j] = (int *)malloc(dim3*sizeof(int));
          if (M[i][j]==NULL) {
             printf("\n\nMAT_3_INT: Memory allocation error\n\n");
             exit(0);
          }
          for (k=0; k<dim3; k++)
              M[i][j][k] = 0;
      }
  }
  return M;
}
double **MAT_DOUBLE(int dim1, int dim2){
  int i,j;
  double **M;
  M = (double **)malloc(dim1*sizeof(double *));
  if (M==NULL) {
     printf("\n\nMAT_DOUBLE: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim1; i++){
      M[i] = (double *)malloc(dim2*sizeof(double));
      if (M[i]==NULL) {
         printf("\n\nMAT_DOUBLE: Memory allocation error\n\n");
         exit(0);
      }
      for (j=0; j<dim2; j++)
          M[i][j] = 0.0;
  }
  return M;
}
void FREE_4(double ****M, int dim1,int dim2, int dim3){
  int i,j,k;
  assert(M!=NULL);
  for (i=0; i<dim1; i++){
      for (j=0; j<dim2; j++){
          for (k=0; k<dim3; k++){
              assert(M[i][j][k]!=NULL);
              free(M[i][j][k]);
          }
          assert(M[i][j]!=NULL);
          free(M[i][j]);
      }
      assert(M[i]!=NULL);
      free(M[i]);
  }
  free(M);
}
void FREE_4_INT(int ****M, int dim){
  int i,j,k;
  assert(M!=NULL);
  for (i=0; i<dim; i++){
      for (j=0; j<dim; j++){
          for (k=0; k<dim; k++){
              assert(M[i][j][k]!=NULL);
              free(M[i][j][k]);
          }
          assert(M[i][j]!=NULL);
          free(M[i][j]);
      }
      assert(M[i]!=NULL);
      free(M[i]);
  }
  free(M);
}
void FREE_4_DOUBLE(double ****M, int dim){
  int i,j,k;
  assert(M!=NULL);
  for (i=0; i<dim; i++){
      for (j=0; j<dim; j++){
          for (k=0; k<dim; k++){
              assert(M[i][j][k]!=NULL);
              free(M[i][j][k]);
          }
          assert(M[i][j]!=NULL);
          free(M[i][j]);
      }
      assert(M[i]!=NULL);
      free(M[i]);
  }
  free(M);
}
void FREE_2_DOUBLE(double **M, int dim)
  {
  int i,j;
  assert(M!=NULL);
  for (i=0; i<dim; i++){
      assert(M[i]!=NULL);
      free(M[i]);
  }
  free(M);
}
void FREE_2_CHAR(char **M, int dim){
  int i,j;
  assert(M!=NULL);
  for (i=0; i<dim; i++){
      assert(M[i]!=NULL);
      free(M[i]);
  }
  free(M);
}
void FREE_2_INT(int **M, int dim){
  int i,j;
  assert(M!=NULL);
  for (i=0; i<dim; i++){
      assert(M[i]!=NULL);
      free(M[i]);
  }
  free(M);
}
void FREE_3_INT(int ***M, int dim){
  int i,j;
  assert(M!=NULL);
  for (i=0; i<dim; i++){
      for (j=0; j<dim; j++){
          assert(M[i][j]!=NULL);
          free(M[i][j]);
      }
      assert(M[i]!=NULL);
      free(M[i]);
  }
  free(M);
}
