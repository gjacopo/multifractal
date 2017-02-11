#ifndef   	WAVELET1D_H_
#define   	WAVELET1D_H_

#ifndef WAVEGAUSS
#define WAVEGAUSS 2 
#endif
#ifndef WAVELORENTZ
#define WAVELORENTZ 3
#endif
#ifndef WAVEMORLET
#define WAVEMORLET 1
#endif
#ifndef WAVEHAAR
#define WAVEHAAR 0
#endif

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

  void wavelet1D_genera( int dimx, double sc, double *wave);
  double wavelet1D_escala( int dimx, double sc );

  int wt1D_transform( double *signal, int dimx, int wav, int order,
		      double sc, double *wtrans );
  void wt1D_projection( double *signal,  double *filt, int dimx, double sc,
			double *wtrans );

  
  void wavelet1D_define_unit( int dimx, double sc, int wav, int ord_der, 
			      double D0, double *wavesig )/*wavelet_1D*/;
  void wavelet1D_define( int dimx, double sc, int wav, int ord_der,
			 double *wave );  
  
  double lorentz1D_coeff( double t, double wav, int ord_der );
  double morlet1D_coeff( double t );
  double gauss1D_coeff( double t, int order );
  double haar1D_coeff( double x );
  
  void wavelet1D_normaliza( int dimx, double sc, int ord_der, double *wave);
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !WAVELET1D_H_ */
