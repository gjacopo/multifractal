/* ===================================
** test.c
** started on Wed Feb 28 11:29:27 2007 
** ===================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#ifndef I_TCHAR 
#define I_TCHAR 1
typedef char TCHAR;
#endif

#ifndef I_TUCHAR
#define I_TUCHAR 2 
typedef unsigned char TUCHAR;
#endif

#ifndef I_TINT
#define I_TINT 3
typedef int TINT;
#endif

#ifndef I_TUINT
#define I_TUINT 4
typedef unsigned int TUINT;
#endif

#ifndef I_TLONG
#define I_TLONG 5
typedef long TLONG;
#endif

#ifndef I_TULONG
#define I_TULONG 6
typedef unsigned long TULONG;
#endif

#ifndef I_TFLOAT
#define I_TFLOAT 9
typedef float TFLOAT;
#endif

#ifndef I_TDOUBLE
#define I_TDOUBLE 10
typedef double TDOUBLE;
#endif



typedef union array {
  int **i; 
  float **f;
  double **d;
} Array;

typedef struct matrix {
  int itype;
  Array a; 
  int x, y;
}
Matrix;




#ifndef CALLOC3D
#define CALLOC3D(p,type,size) ((p)=(type***)calloc(size,sizeof(type**)))
#endif

#ifndef CALLOC2D
#define CALLOC2D(p,type,size) ((p)=(type**)calloc(size,sizeof(type*)))
#endif


#ifndef CALLOC1D
#define CALLOC1D(p,type,size) ((p)=(type*)calloc(size,sizeof(type)))
#endif


int main() {
  Matrix f;
  f.itype = I_TDOUBLE;
  int x=5, y=4;
  int i, j;
  
  CALLOC2D(f.a.d,TDOUBLE,x);
  for(i=0; i<x;i++) {
    CALLOC1D(f.a.d[i],TDOUBLE,y);
    for(j=0; j<y; j++)
      f.a.d[i][j] = (double) i*y + j+1;
      }
  
  for(i=0; i<x; i++) {
    for(j=0; j<y; j++)
      fprintf(stderr,"\t %f",f.a.d[i][j]); 
    fprintf(stderr,"\n"); 
  } 

  return 0;
}
