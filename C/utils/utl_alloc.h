#ifndef UTL_ALLOC_H
#define UTL_ALLOC_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  TypeUn *vector( int xdim );
  TypeUn *centvector( int ldim, int rdim);
  
  double *dvector( int xdim );
  double *centdvector( int ldim, int rdim);
  int *ivector( int xdim );
  char *cvector( int xdim );
  
  double ***matrix3D( int zdim, int ydim, int xdim);
  int free_matrix3D( double ***m, int zdim, int ydim );
  double ***realloc_matrix3D( int zdim0, int ydim0, int xdim0,
			      int  zdim, int ydim, int xdim, double ***pointer );
  /*redimensiona_tritensor*/;
  
  char ***cmatrix3D( int zdim, int ydim, int xdim );
  int free_cmatrix3D( char ***m, int zdim, int ydim );
  
  double **matrix2D( int ydim, int xdim );
  int free_matrix2D( double **m, int ydim );
  
  char **cmatrix2D( int ydim, int xdim );
  int free_cmatrix2D( char **m, int ydim );
  
  int **imatrix2D( int ydim, int xdim );
  int free_imatrix2D( int **m, int ydim );
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
