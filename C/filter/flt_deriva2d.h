/*    derivacion.h -  Version del 24 de Agosto, 2004         */

#ifndef DERIVA_2D_H
#define DERIVA_2D_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  /*     Function prototypes   */
  
  int derivaD_bak(int dimx, int dimy, double **gx, double **gy);
  int modderiva_naif( int dimx, int dimy, double **mu);
  int modderiva( int dimx, int dimy, double **mu);
  int modderiva_line( int dimx, int dimy, double **modg, double theta);
  int modgradderiva( int dimx, int dimy, 
		     double **mu, double **gx, double **gy);
 
  void reconstruct2D( int dimx, int dimy, double **gx, double **gy);
  int reconstruct2D_naif( int dimx, int dimy, double **gx, double **gy);
  int reconstruct2D_FFT( int dimx, int dimy, double **gx, double **gy);

  void gradient2D( int dimx, int dimy,  double **gx, double **gy);
  int gradient2D_naif_bak( int dimx, int dimy,  double **gx, double **gy);
  int gradient2D_naif( int dimx, int dimy,  double **gx, double **gy);
  int gradient2D_FFT( int dimx, int dimy,  double **gx, double **gy);
  int gradient2D_naif_inter( int dimx, int dimy,  double **gx, double **gy);

  void gradient_complex( int dimx, int dimy,  double **gx, double **gy);
  void reconstruct_complex( int dimx, int dimy, double **gx, double **gy);
  void deriva_complex( int dimx, int dimy, int mode, double **gx, double **gy);
  void deriva_complex_naif( int dimx, int dimy, int mode, double **gx, 
			    double **gy );
  void deriva_complex_FFT( int dimx, int dimy, int mode, double **gx, 
			   double **gy );

  int filtro2D_bak( int dimx, int dimy, double norma, double expon,
		    double **data);
  int filtro2D( int dimx, int dimy, double expon, double **data); 
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
