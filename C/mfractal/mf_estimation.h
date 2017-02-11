/********************************************************
 * mf_estimation.h               
 ********************************************************/

#ifndef   	_MF_ESTIMATION_H
# define   	_MF_ESTIMATION_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  double calcula_multifractal( int dimx, double **signal, 
			       double **expon, double *mm );
  double calcula1D_multifractal( int dimx, double *signal,
				 double *expon, double *mm_h );
  double calcula2D_multifractal( int dimx, int dimy, double **signal,
				 double **expon, double *mm_h );

  double Dh_registra( int Nr, double sc0, double *h, double *Dh, double *errDh);
  
  void expon_genera( int leff, double **signal, double **expon );
  int expon1D_genera( int leff, double *serie, double *expon );
  int expon2D_genera( int leff, double **image, double **expon );
  
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !MF_ESTIMATION_H */
