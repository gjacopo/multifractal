#ifndef UTL_OPERATOR_H
#define UTL_OPERATOR_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */


  int fill0( int dimx, int dimy, double **matriz, char**m )/*limpia*/;
  int fill1D0( int dimx, double *vector, char*m )/*limpia_lista*/;
  int fill( int dimx, int dimy, double d, double **matriz, char**m );
  int fill1D( int dimx, double d, double *vector, char*m );
  int ifill( int dimx, int dimy, int i, int **matriz, char**m )/*limpia_int*/;
  int ifill1D( int dimx, int i, int *vector, char*m );
  int cfill( int dimx, int dimy, char val, char **matriz, char**m )/*limpia_char*/;
  int cfill1D( int dimx, char val, char *vector, char*m  ) /*limpia_char_lista*/;

  int copy( int dimx, int dimy, double **fuente, double **destino, char**m )/*asigna*/;
  int copy1D( int dimx, double *fuente, double *destino, char*m) /*asigna_lista*/;
  int ccopy( int dimx, int dimy, char **fuente, char **destino, char**m ) /*asigna_char*/;
  
  int coarseres( int dimx, int dimy, double block, double **source,
		 double **target ) /*coarse_resolution*/;
  int coarseres1D( int dimx, double block, double *source,
		   double *target ) /*coarse_resolution_lista*/;
  int coarseres_mask( int dimx, int dimy, double block, int ivabs,
		      char **source, char **target) /*coarse_resolution_mask*/;
  int ccoarseres( int dimx, int dimy, double block, char **source,
		  char **target ) /*coarse_resolution_char*/;
  int blurres( int dimx, int dimy, int block, double **source, double
	       **target ) /*blur_resolution*/;
  int de4(int dimx, int dimy, double **signal );
  int de4_mirror(int dimx, int dimy, double **signal );
  
  int paddwindow( int dimx, int dimy, int ix0, int iy0, int dx, 
		  int dy, double **data, double **window ) /*recorta_ventana*/;
  int cpaddwindow( int dimx, int dimy, int ix0, int iy0, int dx, 
		   int dy, char **data, char **window)/*recorta_ventana_char*/;
  int unpaddwindow( int dimx, int dimy, int ix0, int iy0, int ix1, 
		    int iy1, int tx, int ty, double **window, 
		    double **total )/*remete_ventana*/; 
  
  int op_diff( int dimx, int dimy, double **fuente, double **destino, char**m ) /*asigna_resta*/;
  int op_diff1D( int dimx, double *fuente, double *destino, char*m ) /*asigna_resta_lista*/;
  int op_add( int dimx, int dimy, double **fuente, double **destino, char**m ) /*asigna_suma*/;
  int op_combine( int dimx, int dimy, double escalar, double **fuente,
		  double **destino, char**m) /*asigna_combina*/;
  int op_fabs(int dimx, int dimy, double **fuente, double **destino, char**m );
  
  int op_shift(int dimx, int dimy, double shift, double **funcion, char**m ) /*asigna_desplaza*/;
  int op_scale(int dimx, int dimy, double scale, double **funcion, char**m ) /*asigna_escala*/;
  int cop_scale( int dimx, int dimy, int scale, char **funcion, char**m ); 
  int op_divide(int dimx, int dimy, double **fuente, double **destino, char**m ) /*asigna_divide*/;
  int op_divide1D( int dimx, double *fuente, double *destino, char*m ) /*asigna_divide_lista*/; 
  int op_multiply(int dimx, int dimy, double **fuente, double **destino, char**m ) 
    /*asigna_multiplica*/;
  int op_multiply1D(int dimx, double *fuente, double *destino, char*m ) 
    /*asigna_multiplica_lista*/;
  int opvec_divide(int dimx, int dimy, double **fx, double **fy, 
		   double **dx, double **dy, char**m ) /*asigna_divide_vec*/;
  int opvec_multiply(int dimx, int dimy, double **fx, double **fy, 
		     double **dx, double **dy, char**m ) /*asigna_multiplica_vec*/;
  int op_modulus( int dimx, int dimy, double **vx, double **vy,
		  double **mod, char**m ) /*asigna_modulo*/;
  int opvec_norma( int dimx, int dimy, double **vx, double **vy, char**m ) /*normaliza_vector*/;
  int opvec_rotate(int dimx, int dimy, double angle, double **vx, double **vy, char**m ) /*gira*/;
  int op_affine( int dimx, int dimy, 
		 double scale, double escalar, double **funcion, char**m ) /*asigna_afin*/;
  int op_addsq( int dimx, int dimy, double **fuentes, double **destino, char**m ) /*asigna_multiplica_cuadra*/;
  int op_sqrtaddsq( int dimx, int dimy, double **fuentes, double **destino, char**m ) 
    /*asigna_multiplica_racina*/;
  int op_square( int dimx, int dimy, double **fuentes, double **destino, char**m ) /*asigna_cuadra*/;
  int op_logaddsq( int dimx, int dimy, double **fuentes, double **destino, char**m ) 
    /*asigna_log_multiplica_cuadra*/;
  int op_log( int dimx, int dimy, double **fuentes, double **destino, char**m ) /*asigna_log*/;
  int opvec_sqdiff( int dimx, int dimy, 
		    double **ux, double **uy, double **vx, double **vy,
		    double **destino, char**m ) /*asigna_cuadradiff_vec*/;
  int opvec_orient( int dimx, int dimy, double **ax, double **ay,
		    double **gx, double **gy, double **orient, char**m ) /*asigna_orientation*/;
  
  
  /* For all the following tests, the test depends on the sign of the
   * variable sign:
   * if sign>0, check:      cont>thres
   * i.e: if cont>thres, mask=TRUE
   *      else           mask=FALSE
   * if sign<0, check:      -cont>-thres <=> cont<thres
   * i.e. the opposite     */
  int op_binary( int dimx, int dimy, double **cont, char** mask,
		 int sign, double thres );
  int op2ptr_threshold( int dimx, int dimy, double **cont, double* c,
			int sign, double thres );
  int op2ptr_binary( int dimx, int dimy, double **cont, char* c,
		     int sign, double thres );
  
  /* if sign>0, check when are both TRUE
   * i.e: if c1=c2=TRUE, then c2=TRUE
   *      else                c2=FALSE
   * if sign<0, check  when are both FALSE
   * i.e: if c1=c2=TRUE, then c2=FALSE
   *      else                c2=TRUE
   */
  int opbin_and( int dimx, int dimy, int sign, 
		 char **fuen, char **dest );
  int op_rdiffratio( double **i1, double**i2, double **i3,
		     int dimx, int dimy,
		     double **o1, double **o2, double **o3 );
  
  int op_find( int dimx, int dimy, double **cont, int nvalues, double *mm);
  int cop_find( int dimx, int dimy, char **m, int nvalues, int *mm);
  
  int op_change( int dimx, int dimy, double **cont, double i0, double i1, char **m );
  int cop_change( int dimx, int dimy, char **cont, int i0, int i1, char **m );

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
