/*       FFT.h - Version del 21 de Septiembre 2005 */

#ifndef _FFT2D_H_
#define _FFT2D_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  /*       Function prototypes      */
  
  void Fourier2D(int dimx, int dimy, double **funcionR, double **funcionI, 
		 int signo);
  void convuelve2D( int dimx, int dimy, double **f1, double **f2, 
		    double **salida);
  void convuelto2D( int dimx, int dimy, double **f1, double **f2);
  void convuelto2D_vec( int dimx, int dimy, double **f1x, double **f1y,
			double **f2x, double **f2y);
  void deconvuelve2D( int dimx, int dimy, double **f1, double **f2, 
		       double **salida);
  void deconvuelto2D( int dimx, int dimy, double **f1, double **f2);
  void deconvuelto2D_vec( int dimx, int dimy, double **f1x, double **f1y,
			  double **f2x, double **f2y);
  
  void FFT2D(int dimx, int dimy, double **funcionR, double **funcionI, 
	     int signo);
  void FFThorizontal(int dimx, int dimy, double **funcionhR, double **funcionhI, 
		     int signo);
  void FFTvertical(int dimx, int dimy, double **funcionvR, double **funcionvI, 
		   int signo);
  void FFFT2D(int dimx, int dimy, double **funcionR, double **funcionI, 
	      int signo);
  int FFFThorizontal(int dimx, int dimy, double **funcionhR, double **funcionhI,
		      int signo);
  void PFFThorizontal( int dima, int dimb, int dimc, int iy, 
		       double **uinR, double **uinI, 
		       double **uoutR, double **uoutI, int sign);
  int FFFTvertical(int dimx, int dimy, double **funcionvR, double **funcionvI, 
		    int signo);
  void PFFTvertical( int dima, int dimb, int dimc, int ix, 
		     double **uinR, double **uinI, 
		     double **uoutR, double **uoutI, int sign);
  
  int convuelve_FFT( int dimx, int dimy, double **f1, double **f2, 
		      double **salida);
  int convuelto_FFT( int dimx, int dimy, double **f1, double **f2);
  int convuelto_vec_FFT( int dimx, int dimy, double **f1x, double **f1y,
			  double **f2x, double **f2y);
  int deconvuelve_FFT( int dimx, int dimy, double **f1, double **f2, 
			double **salida);
  int deconvuelto_FFT( int dimx, int dimy, double **f1, double **f2);
  int deconvuelto_vec_FFT( int dimx, int dimy, double **f1x, double **f1y,
			    double **f2x, double **f2y);
  
  int band_pass( int dimx, int dimy, double fmin, double fmax, 
		  double **func);
  
  int convuelve_FFFT( int dimx, int dimy, double **f1, double **f2, 
		       double **salida);
  int convuelto_FFFT( int dimx, int dimy, double **f1, double **f2);
  int convuelto_vec_FFFT( int dimx, int dimy, double **f1x, double **f1y,
			   double **f2x, double **f2y);
  int deconvuelve_FFFT( int dimx, int dimy, double **f1, double **f2, 
			 double **salida);
  int deconvuelto_FFFT( int dimx, int dimy, double **f1, double **f2);
  int deconvuelto_vec_FFFT( int dimx, int dimy, double **f1x, double **f1y,
			     double **f2x, double **f2y);

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
