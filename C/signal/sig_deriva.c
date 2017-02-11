/* ===================================
** sig_deriva.c
** started on Tue Jan 30 11:40:23 2007 
** ===================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <filter.h>
#include <flt_deriva2d.h>
#include <flt_fft2d.h>

#include <signal.h>
#include <sig_fft.h>
#include <sig_deriva.h>

#ifdef _PARSE_FILTER_PARAMETERS_
#include <flt_parse.h>
extern ParFILT *p_fil;
#ifndef redimensiona
#define redimensiona(x) (p_fil->flag_memory?(x):dimensiona(x))
#endif

/***************************************************************************/
int modgradient( Signal *data ) {
  /***************************************************************************
   * A simple routine to obtain the modulus of the gradient.
   * High gradients associated to masked points are filtered out.
   ***************************************************************************/
  Matrix gx=_Matrix,gy=_Matrix;
  Matrix filtered=_Matrix;
  double mean_gradient;
  int ix,iy;
  
  if((data->dimx<1) || (data->dimy<1)) Exit(ERROR);
  
  initiate_matrix(data,&gx);
  create_matrix(data.dimx,data.dimy,0,&gy);
  gradient(gx,gy);
  
  for (iy=0;iy<data.dimy;iy++) {
    for (ix=0;ix<data.dimx;ix++) {
      data.signal[iy][ix] = 
	sqrt(gx.signal[iy][ix]*gx.signal[iy][ix]+gy.signal[iy][ix]*gy.signal[iy][ix]);
    }
   }

   destroy_matrix(&gx);
   destroy_matrix(&gy);

   /* Filtering the neighborhood of the mask and the edges */
   create_matrix(data.dimx,data.dimy,1,&filtered);
   assign_matrix(data, filtered);
   
   // Exclude the edges:
   
   for (iy=0,ix=0;iy<data.dimy;iy++)           filtered.Mask[iy][ix] = MASKED;
   for (iy=0,ix=data.dimx-1;iy<data.dimy;iy++) filtered.Mask[iy][ix] = MASKED;
   if (data.dimy>1)     {
     for (iy=0,ix=0;ix<data.dimx;ix++)           filtered.Mask[iy][ix] = MASKED;
     for (iy=data.dimy-1,ix=0;ix<data.dimx;ix++) filtered.Mask[iy][ix] = MASKED;
   }
   
   // Exclude the mask's neighbors:
   if (data.Mask!=NULL)   
     for (iy=1;iy<data.dimy-1;iy++)      
       for (ix=1;ix<data.dimx-1;ix++)	 
         if (data.Mask[iy][ix]==MASKED)         {
	   filtered.Mask[iy-1][ix] = MASKED;
	   filtered.Mask[iy+1][ix] = MASKED;
	   filtered.Mask[iy][ix-1] = MASKED;
	   filtered.Mask[iy][ix+1] = MASKED;
         }

   if (MEANMASK)     {
     // The mean value of unmasked points is assigned to the masked points
     mean_gradient = matrix_mean(filtered);
     for (ix=0; ix<filtered.dimx*filtered.dimy; ix++) 
       if (filtered.pm[ix] == MASKED) filtered.p[ix] = mean_gradient;
   }
   
   assign_matrix(filtered, data);
   
   destroy_matrix(&filtered);
   
   return OK;
} // end of modgradient


/***************************************************************************/
int reconstruct_signal( Signal* gx, Signal* gy){
  /***************************************************************************
   * Common access point to perform gradient reconstruction
   * A routine performing reconstruction from gradient, as in: 
   *        Turiel and del Pozo, IEEE Trans. Im. Proc. (2002). 
   ***************************************************************************/
  int error;
  
#ifdef _PARSE_FILTER_PARAMETERS_
  switch(p_fil->mode_deriva)    
#else
    switch(MODE_DERIVA) 
#endif/* !_PARSE_FILTER_PARAMETERS_ */
      {
      case 0:
	error=reconstruct_signal_FFT(gx,gy);
	break;
      case 1:
	error=reconstruct_signal_naif(gx,gy);
	break;
      default:
	Exit("Unrecognized derivative mode\n");
      }
  
  return error;
} // end of reconstruct_signal


/***************************************************************************/
int reconstruct_signal_FFT( Signal* gx, Signal* gy ) {
  /***************************************************************************/
  Signal *gxI,*gyI;
  double aux,x,y,f;
  int ix,iy;
  int dimx=gx->dimx,dimy=gx->dimy;
  
  TrackError(compare_dimsignal(gx,gy), "Signals dimensions do not match");  
  
  TrackNullAlloc( gxI = create_signal(dimx,dimy,NULL) );
  TrackNullAlloc( gyI = create_signal(dimx,dimy,NULL) );
  
  FFT_signal(gx,gxI,-1);
  FFT_signal(gy,gyI,-1);
  
  for (iy=0;iy<gx.dimy;iy++)    {
    y = ((double)iy) / ((double)gx.dimy);
    if ((iy>0) && (iy>gx.dimy/2)) y -= 1.;
    y = 2.*sin(M_PI*y);

    for(ix=0;ix<gx.dimx;ix++)   {
      x = ((double)ix) / ((double)gx.dimx);
      if ((ix>0) && (ix>gx.dimx/2)) x -= 1.;
      x = 2.*sin(M_PI*x);
      
      f = x*x + y*y;
      if (f>MIN_NUM_RESOLUTION) {
         aux = (x*gxI->signal[iy][ix] + y*gyI->signal[iy][ix]) / f;
         gxI->signal[iy][ix] = -(x*gx->signal[iy][ix] + y*gy->signal[iy][ix]) / f;
         gx->signal[iy][ix] = aux;
      } else {
         gx->signal[iy][ix] = 0.;
         gxI->signal[iy][ix] = 0.;
      }

      gy->signal[iy][ix] = 0.;
      gyI->signal[iy][ix] = 0.;
    }
  }
  
  FFT_signal(gx,gxI,1);
  
  /* Memory release and end  */
  free_signal(gxI);
  free_signal(gyI);
  
  return OK;
} // end of reconstruct_signal_FFT


/***************************************************************************/
int reconstruct_signal_naif( Signal* gx, Signal* gy ) {
  /***************************************************************************/
  Signal *gxI,*gyI;
  double dxR,dxI,dyR,dyI,modx,mody;
  double auxR,auxI,x,y;
  int ix,iy;
  int dimx=gx->dimx,dimy=gx->dimy;

  TrackError(compare_dimsignal(gx,gy), "Signals dimensions do not match");  

  TrackNullAlloc( gxI = create_signal(dimx,dimy,NULL) );
  TrackNullAlloc( gyI = create_signal(dimx,dimy,NULL) );

  FFT_signal(gx,gxI,-1);
  FFT_signal(gy,gyI,-1);
  
  for (iy=0;iy<gx.dimy;iy++) {
    y = ((double)iy) / ((double)gx.dimy);
    dyR = cos(2.*M_PI*y) -1.;
    dyI = -sin(2.*M_PI*y);
    mody = dyR*dyR + dyI*dyI;

    for (ix=0;ix<gx.dimx;ix++)      {
      x = ((double)ix) / ((double)gx.dimx);
      dxR = cos(2.*M_PI*x) -1.;
      dxI = -sin(2.*M_PI*x);
      modx = dxR*dxR + dxI*dxI;
      
      if (modx+mody > MIN_NUM_RESOLUTION)         {
	auxR = (dxR*gx->signal[iy][ix] - dxI*gxI->signal[iy][ix] +
		dyR*gy->signal[iy][ix] - dyI*gyI->signal[iy][ix]) /
	  (modx+mody);
	auxI = (dxR*gxI->signal[iy][ix] + dxI*gx->signal[iy][ix] +
		dyR*gyI->signal[iy][ix] + dyI*gy->signal[iy][ix]) / 
	  (modx+mody);
      } else {
	auxR = auxI = 0.;
      }
      gx->signal[iy][ix] = auxR;
      gxI->signal[iy][ix] = auxI;
      gy->signal[iy][ix] = 0.;
      gyI->signal[iy][ix] = 0.;
    }
  }
  
  FFT_signal(gx,gxI,1);

  /* Memory release and end  */
  free_signal(gxI);
  free_signal(gyI);
  
  return OK;
} // end of reconstruct_signal_naif
