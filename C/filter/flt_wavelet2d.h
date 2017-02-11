#ifndef   	WAVELET2D_H_
#define   	WAVELET2D_H_


#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

  void wavelet2D_genera( int dimx, int dimy, double sc, double **wave);
  double wavelet2D_escala( int dimx, int dimy, double sc );

  double wavelet2D_escala_line( int dimx, int dimy, double sc );

  void wavelet2D_genera_line( int dimx, int dimy, double sc, double **wave);

  int wt2d_transform( int dimx, int dimy, double a0, double quanto, 
		      double orden, int degree, int npunt,
		      double **muR, double **muI,
		      double **expon, double *minimo, double *maximo)/*wavelet*/;

  void test_wavelet( int dimx, int dimy, double a,  double orden, double **psi);
  int convolverwave_log( int dimx, int dimy, /* int xeff, int yeff, */
			double **muR, double **muI, double scale, int degree, 
			 double orden, double **Tmu);
  void lorentz2D( int dimx, int dimy, double a, int degree, double orden,
		  double **psi );
  void gauss2D( int dimx, int dimy, double a, double **psi );
  void log2D_1( int dimx, int dimy, double a, double **psi );
  void log2D_2( int dimx, int dimy, double a, double **psi );
  void log2D_3( int dimx, int dimy, double a, double **psi );

  void wavelet2D_define_unit( int dimx, int dimy, double sc, int wav, int ord_der, 
			      double D0, double **wavesig )/*wavelet_2D*/;
  void wavelet2D_define( int dimx, int dimy, double sc, int wav, int ord_der, 
		       double **wave );
  double haar2D_coeff( double x, double y );
  double lorentz2D_coeff( double x, double y, double wav, int ord_der );
  double morlet2D_coeff( double x, double y );
  double gauss2D_coeff( double x, double y, int ord_der );


  void wavelet2D_define_line( int dimx, int dimy, double sc, 
			      double theta, int wav, int ord_der,
			      double **wave);
  void wavelet2D_normaliza( int dimx, int dimy, double sc, int ord_der, 
			    double **wave);
  
  double escala2D_lineal( int dimx, int dimy);

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !WAVELET2D_H_ */
