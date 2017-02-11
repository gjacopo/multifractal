/*    from derivacion.h -  Version del 24 de Agosto, 2004         */

#ifndef DERIVA_1D_H
#define DERIVA_1D_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  /*     Function prototypes   */
  
  void gradient1D( int dimx, double *gx);
  void reconstruct1D( int dimx, double *gx);
  
  void gradient1D_naif( int dimx, double *gx);
  int reconstruct1D_naif( int dimx, double *gx);
  int gradient1D_FFT( int dimx, double *gx);
  int reconstruct1D_FFT( int dimx, double *gx);
  int gradient1D_naifinter( int dimx, double *gx);
  
  int filtro1D( int dimx,double expon, double *funcion);
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
