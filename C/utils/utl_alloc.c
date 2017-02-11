#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_alloc.h>


/***************************************************************************/
TypeUn *vector( int xdim ){
/***************************************************************************/
  TypeUn *v;
  
  if ((v=(TypeUn*)calloc(xdim,sizeof(TypeUn))) == NULL)
    return NULL;
  
  return v;
}

/***************************************************************************/
double *dvector( int xdim ){
/***************************************************************************/
  double *v;
  
  if ((v=(double*)calloc(xdim,sizeof(double))) == NULL)
    return NULL;
  
  return v;
}

/***************************************************************************/
TypeUn *centvector( int ldim, int rdim){
/***************************************************************************/
  TypeUn *v;
  
  if ((v=(TypeUn*)calloc(ldim+rdim,sizeof(TypeUn))) == NULL)
    return NULL;
  
  return v+ldim;
}

/***************************************************************************/
double *centdvector( int ldim, int rdim){
/***************************************************************************/
  double *v;
  
  if ((v = (double *)calloc(ldim+rdim,sizeof(double))) == NULL)
    return NULL;
  
  return v+ldim;
}

/***************************************************************************/
int *ivector( int xdim ){
/***************************************************************************/
  int *v;
  
  if ((v=(int*)calloc(xdim,sizeof(int))) == NULL)
    return NULL;
  
  return v;
}

/***************************************************************************/
int *iuvector( int xdim ){
/***************************************************************************/
  int *v;
  
  if ((v=(unsigned int*)calloc(xdim,sizeof(unsigned int))) == NULL)
    return NULL;
  
  return v;
}

/***************************************************************************/
char *cvector( int xdim ){
/***************************************************************************/
  char *v;
  
  if ((v=(char*)calloc(xdim,sizeof(char))) == NULL)
    return NULL;
  
  return v;
}

/***************************************************************************/
char *ucvector( int xdim ){
/***************************************************************************/
  char *v;
  
  if ((v=(unsigned char*)calloc(xdim,sizeof(unsigned char))) == NULL)
    return NULL;
  
  return v;
}

/***************************************************************************/
TypeUn ***matrix3D( int zdim, int ydim, int xdim) {
  /***************************************************************************/
  TypeUn ***m;
  int j,k;

  if((m=(TypeUn***) calloc(zdim,sizeof(TypeUn**))) == NULL)
    return NULL;
  
  for(k=0; k<zdim; ++k) 	
    if((m[k]=matrix2D( ydim, xdim)) == NULL)
      return NULL;
  
  return m;
}


/***************************************************************************/
double ***dmatrix3D( int zdim, int ydim, int xdim) {
  /***************************************************************************/
  double ***m;
  int j,k;

  if((m=(double ***) calloc(zdim,sizeof(double **))) == NULL)
    return NULL;
  
  for(k=0;k<zdim;++k) 	
    if((m[k]=dmatrix2D( ydim, xdim)) == NULL)
      return NULL;
  /*
  for(k=0;k<zdim;++k) 	{
  if( (m[k]=(double **) calloc( ydim,sizeof(double *))) == NULL)
  return NULL;
  for(j=0;j<ydim;++j) m[k][j]=(double *) calloc(xdim,sizeof(double));
  }
  */
  
  return m;
}


/***************************************************************************/
int free_matrix3D( TypeUn ***m, int zdim, int ydim ) {
  /***************************************************************************/
  int j,k;

  for(k=0; k<zdim; ++k)	{
    for(j=0; j<ydim; ++j) free(m[k][j]);
    free(m[k]);
  }
  free(m);

  return OK;
}

/***************************************************************************/
int free_dmatrix3D( double ***m, int zdim, int ydim ) {
  /***************************************************************************/
  int j,k;

  for(k=0; k<zdim; ++k)	{
    for(j=0; j<ydim; ++j) free(m[k][j]);
    free(m[k]);
  }
  free(m);

  return OK;
}

/***************************************************************************/
TypeUn ***realloc_matrix3D( int zdim0, int ydim0, int xdim0,
			    int  zdim, int ydim, int xdim, TypeUn ***pointer )
  /*redimensiona_tritensor*/ {
  /***************************************************************************/
  TypeUn ***m;
  int i,j,k;

  m=matrix3D(zdim,ydim,xdim);

  for(k=0;k<zdim;k++)
    for(j=0;j<ydim;j++)
      for(i=0;i<xdim;i++)	{
	if((k<zdim0)&&(j<ydim0)&&(i<xdim0))
	  m[k][j][i]=pointer[k][j][i];
	else m[k][j][i]=0.;
      }
    
  free_matrix3D(pointer,zdim,ydim);

  return m;
}

/***************************************************************************/
double ***realloc_dmatrix3D( int zdim0, int ydim0, int xdim0,
			    int  zdim, int ydim, int xdim, double ***pointer )
  /*redimensiona_tritensor*/ {
  /***************************************************************************/
  double ***m;
  int i,j,k;

  m=dmatrix3D(zdim,ydim,xdim);

  for(k=0;k<zdim;k++)
    for(j=0;j<ydim;j++)
      for(i=0;i<xdim;i++)	{
	if((k<zdim0)&&(j<ydim0)&&(i<xdim0))
	  m[k][j][i]=pointer[k][j][i];
	else m[k][j][i]=0.;
      }
    
  free_dmatrix3D(pointer,zdim,ydim);

  return m;
}

/***************************************************************************/
char ***cmatrix3D( int zdim, int ydim, int xdim ) {
  /***************************************************************************/
  char ***m;
  int j,k;

  m=(char ***) calloc(zdim,sizeof(char **));
  for(k=0;k<zdim;++k) 	{
    m[k]=(char **) calloc( ydim,sizeof(char *));
    for(j=0;j<ydim;++j) m[k][j]=(char *) calloc(xdim,sizeof(char));
  }
	
  return m;
}

/***************************************************************************/
int free_cmatrix3D( char ***m, int zdim, int ydim ) {
  /***************************************************************************/
  int j,k;

  for(k=0; k<zdim; ++k)	{
    for(j=0; j<ydim; ++j) free(m[k][j]);
    free(m[k]);
  }
  free(m);

  return OK;
}


/***************************************************************************/
TypeUn **matrix2D( int ydim, int xdim ) {
  /***************************************************************************/
  TypeUn **m;
  int j;

  m=(TypeUn **) calloc(ydim,sizeof(TypeUn *));
  for(j=0;j<ydim;++j) m[j]=(TypeUn*) calloc( xdim,sizeof(TypeUndouble));
	
  return m;
}
/***************************************************************************/
double **dmatrix2D( int ydim, int xdim ) {
  /***************************************************************************/
  double **m;
  int j;

  m=(double **) calloc(ydim,sizeof(double *));
  for(j=0;j<ydim;++j) m[j]=(double *) calloc( xdim,sizeof(double));
	
  return m;
}

/***************************************************************************/
int free_matrix2D( TypeUn **m, int ydim ) {
  /***************************************************************************/
  int j;

  for(j=0; j<ydim; ++j) free(m[j]);
  free(m);

  return OK;
}

/***************************************************************************/
int free_dmatrix2D( double **m, int ydim ) {
  /***************************************************************************/
  int j;

  for(j=0; j<ydim; ++j) free(m[j]);
  free(m);

  return OK;
}


/***************************************************************************/
char **cmatrix2D( int ydim, int xdim ) {
  /***************************************************************************/
  char **m;
  int j;

  m=(char **) calloc(ydim,sizeof(char *));
  for(j=0;j<ydim;++j) m[j]=(char *) calloc( xdim,sizeof(char));
	
  return m;
}

/***************************************************************************/
int free_cmatrix2D( char **m, int ydim ) {
  /***************************************************************************/
  int j;

  for(j=0; j<ydim; ++j) free(m[j]);
  free(m);
  return OK;
}

/***************************************************************************/
int **imatrix2D( int ydim, int xdim ) {
  /***************************************************************************/
  int **m;
  int j;

  m=(int **) calloc(ydim,sizeof(int *));
  for(j=0;j<ydim;++j) m[j]=(int *) calloc( xdim,sizeof(int));
	
  return m;
}

/***************************************************************************/
int free_imatrix2D( int **m, int ydim ) {
  /***************************************************************************/
  int j;

  for(j=0; j<ydim; ++j) free(m[j]);
  free(m);
  return OK;
}


