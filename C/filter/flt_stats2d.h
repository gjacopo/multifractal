/*	Version del 14 de Junio, 2002		*/
#ifndef STATS2D_H
#define STATS2D_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

  double cuantil2D( int dimx, int dimy, double **data, double prob0);
  double cuantil2D_sort( int dimx, int dimy, double **data, double prob0);
  double cuantil2D_histo( int dimx, int dimy, double **data, double prob0);
  
  double moda(int dimx, int dimy, double **datos)/*moda_2D*/;
  int mediadisp( int dimx, int dimy, double **cont, 
		 double *media, double *disp, char**m );
  double dispersion( int dimx, int dimy, double **cont, char**m );
  double dispersion_vec( int dimx, int dimy, double **vx, double **vy, char**m );
  double media( int dimx, int dimy, double **cont, char**m ); 
  double media_expon( int dimx, int dimy, double factor, double **cont, 
		      char **m );
  
  int iextrema( int dimx, int dimy, int **cont, int *mm, char **m );
  int cextrema( int dimx, int dimy, char **cont, int *mm, char **m );
  int extrema( int dimx, int dimy, double **cont, double *mm, char**m );
  int extrema_update( int dimx, int dimy, double **in, double *min, double *max, 
		      char **m );
  
  int norma( int dimx, int dimy, double **cont, char**m );
  double anorma( int dimx, int dimy, double **data, char**m );
  double anorma1( int dimx, int dimy, double **cont, char**m );
  double anorma1_line( int dimx, int dimy, double **data, double theta, char**m );
  double anorma2( int dimx, int dimy, double **cont, char**m );
  int denorma( int dimx, int dimy, double norma, double **data, char**m );
  
  double extrema_local( double **cont, int *iix, int *iiy, 
			double mms[2], int winsize );
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
