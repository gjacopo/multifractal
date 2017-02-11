#ifndef FFT1D_H
#define FFT1D_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

  void Fourier1D(int dimx, double *funcionR, double *funcionI, int signo);
  
  void FFT1D(int dim, double *funcionR, double *funcionI, int signo);
  int FFFT1D( int dim, double *funcionR, double *funcionI, int signo);
  void PFFT1D( int dima, int dimb, int dimc, double *uinR, double *uinI, 
	       double *uoutR, double *uoutI, int sign);
  
  int convuelto1D( int dimx, double *f1, double *f2);
  int convuelve1D( int dimx, double *f1, double *f2, double *salida);

  void FFT(int dim, double *funcionR, double *funcionI, int signo);

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
