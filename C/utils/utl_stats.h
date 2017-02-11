#ifndef UTL_STATS_H
#define UTL_STATS_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

  double fit(double *x, double *y, int n, double *a, double *b, double *corr);
  double fit_best(double *x, double *y, int n, int hmin, 
		  double **repx, double **repy, int *N);
  int fit_line( double *yValues, double *xValues, int dimx, double *slope );
  
  int dimensiona( int dim);
  int adimensiona( int size);
  int adimensiona_pos( int size );
  
  double fMax( double a, double b);
  double fMin( double a, double b);
  int Max( int a, int b);
  int Min( int a, int b);
  int Mod( int a, int b);
  int Round( double a);
  
  void C_mult( double a1, double b1, double a2, 
	       double b2, double *a, double *b );
  void C_sqrt( double a0, double b0, double *a, double *b);
  
  double TanH( double x);
  double angulo( double dx, double dy);
  
  int resolution( int winsize, int limit );
  int vecindex( int *vecii, int i, int winsize, int dim, int per );
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
