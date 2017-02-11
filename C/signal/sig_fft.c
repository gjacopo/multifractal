/* ===================================
** sig_fft.c
** started on Tue Jan 30 11:40:23 2007 
** ===================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <filter.h>
#include <flt_fft2d.h>

#include <signal.h>
#include <sig_fft.h>


/***************************************************************************/
int FFT_signal( Signal* funcR, Signal* funcI, int sign ) {
  /***************************************************************************
   * To calculate Fourier transform by a fast algorithm. This version used the 
   * decomposition in prime factors of the dimensions dimx and dimy of the
   * analyzed data (provided as two matrices, one for the real part and another
   * for the imaginary part). For that reason, no requirement dimx=2^n is
   * imposed. Notice however that if some prime factors are large, the routine
   * will be very slow; in such instances think about padding your signal with
   * constant values up to a nice dimension decomposing in small prime factors.
   * Indeed this could be done inside the program; it is up to you to introduce
   * this improvement. Hints can be requested to the authors.
   ***************************************************************************/
  
  TrackError(compare_dimsignal(funcR,funcI),
	     "Signals dimensions do not match");  
  
  TrackError( FFThorizontal_signal(funcR,funcI,sign),
	      "Error performing horizontal FFT");
  
  if (funcR->dimy > 1)      
    TrackError( FFTvertical_signal(funcR,funcI,sign),
		"Error performing vertical FFT");
  
  return OK;
} // end of FFT_signal


/***************************************************************************/
int FFTvertical_signal( Signal* funcR, Signal* funcI, int sign ) {
  /***************************************************************************/
  /* Auxiliary routine called by FFT */
  Signal *workR,*workI;
  int inu,ix,iy;
  int dimx=funcR->dimx,dimy=funcR->dimy;
  int dima,dimb,dimc;
  
  TrackNullAlloc( workR = create_signal(dimx,dimy,NULL) );
  TrackNullAlloc( workI = create_signal(funcI->dimx,funcI->dimy,NULL) );
  
  for( ix=0; ix<dimx; ix++ ) {
    dima = 1;
    dimb = dimy;
    dimc = 1;
   
    inu = 1; 
    while (dimb > 1) {
      dima = dimc * dima;
      dimc = 2;
      while (dimb%dimc!=0) dimc++;
      dimb = dimb / dimc;
      
      if (inu == 1) 
	PFFTvertical( dima, dimb, dimc, ix, funcR->signal, funcI->signal,
		      workR->signal, workI->signal, sign );
      else        
	PFFTvertical( dima, dimb, dimc, ix, workR->signal, workI->signal,
		      funcR->signal, funcI->signal,sign);
      
      inu = 1 - inu;
    }
    
    if(inu == 0)
      for(iy=0;iy<dimy;iy++)     {
	funcR->signal[iy][ix] = workR->signal[iy][ix];
	funcI->signal[iy][ix] = workI->signal[iy][ix];
      }
    
    for(iy=0;iy<dimy;iy++)    {
      funcR->signal[iy][ix] /= sqrt((double)dimy);
      funcI->signal[iy][ix] /= sqrt((double)dimy);
    }
  }
  
  free_signal(workR);
  free_signal(workI);
  
  return OK;
} // end of FFTvertical_signal


/***************************************************************************/
int FFThorizontal_signal( Signal* funcR, Signal* funcI, int sign ) {
  /***************************************************************************/
  /* Auxiliary routine called by FFT */
  
  Signal *workR,*workI;
  int inu,ix,iy;
  int dimx=funcR->dimx,dimy=funcR->dimy;
  int dima,dimb,dimc;

  TrackNullAlloc( workR = create_signal(dimx,dimy,NULL) );
  TrackNullAlloc( workI = create_signal(funcI->dimx,funcI->dimy,NULL) );

  for (iy=0;iy<dimy;iy++) {
    dima=1;
    dimb=dimx;
    dimc=1;

    inu=1;
    while (dimb>1)    {
      dima=dimc*dima;
      dimc=2;
      while (dimb%dimc!=0) dimc++;
      dimb=dimb/dimc;

      if (inu == 1) 
	PFFThorizontal( dima, dimb, dimc, iy, funcR->signal, funcI->signal,
		       workR->signal, workI->signal, sign);
     else        
       PFFThorizontal( dima, dimb, dimc, iy, workR->signal, workI->signal, 
		       funcR->signal, funcI->signal, sign);

     inu = 1 - inu;
    }
    
    if(inu==0)    
      for (ix=0;ix<dimx;ix++)     {
	funcR->signal[iy][ix] = workR->signal[iy][ix];
	funcI->signal[iy][ix] = workI->signal[iy][ix];
      }
    
    for(ix=0;ix<dimx;ix++)    {
      funcR->signal[iy][ix] /= sqrt((double)dimx);
      funcI->signal[iy][ix] /= sqrt((double)dimx);
    }
  }
  
  free_signal(workR);
  free_signal(workI);
  
  return OK;
} // end of FFThorizontal_signal
