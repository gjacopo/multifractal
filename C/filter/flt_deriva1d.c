#include <stdio.h>
#include <math.h>

/* Personnal libraries */
#include <utils.h>
#include <utl_alloc.h>
#include <utl_stats.h>
#include <utl_operator.h>		

#include <filter.h>
#include <flt_fft1d.h>
#include <flt_stats1d.h>		

#include <flt_deriva1d.h>

#ifdef _PARSE_FILTER_PARAMETERS_
#include <flt_parse.h>
extern ParFILT *p_fil;
#ifndef redimensiona
#define redimensiona(x) (p_fil->flag_memory?(x):dimensiona(x))
#endif

#else
#ifndef redimensiona
#define redimensiona(x) (FLAG_MEMORY?(x):dimensiona(x))
#endif

#endif/* !_PARSE_FILTER_PARAMETERS_ */


/****************************************************************************/
void gradient1D( int dimx,  double *gx) {
/****************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  switch(p_fil->mode_deriva) {
#else
    switch(MODE_DERIVA) {    
#endif/* !_PARSE_FILTER_PARAMETERS_ */
  case 0:
    gradient1D_FFT(dimx,gx);
    break;
  case 1:
    gradient1D_naif(dimx,gx);
    break;
  case 2:
    gradient1D_naifinter(dimx,gx);
    break;
  default:
    Exit("No definition associated to the derivative mode");
  }
}

/****************************************************************************/
void reconstruct1D( int dimx, double *gx) {
/****************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  switch(p_fil->mode_deriva) {
#else
    switch(MODE_DERIVA)  {    
#endif/* !_PARSE_FILTER_PARAMETERS_ */
  case 0:
    reconstruct1D_FFT(dimx,gx);
    break;
  case 1:
    reconstruct1D_naif(dimx,gx);
    break;
  case 2:
    reconstruct1D_naif(dimx,gx);
    break;
  default:
    Exit("No definition associated to the derivative mode");
  }
}

/****************************************************************************/
void gradient1D_naif( int dimx, double *gx) {
/****************************************************************************/

  double cierre;
  int ix;

  cierre=gx[0]-gx[dimx-1];
  for(ix=0;ix<dimx-1;ix++) gx[ix]=gx[ix+1]-gx[ix];

  gx[dimx-1]=cierre;
}

/****************************************************************************/
int reconstruct1D_naif( int dimx, double *gx) {
/****************************************************************************/

  double *aux;
  int ix;

  TrackNullAlloc( aux=(double*)calloc(dimx,sizeof(double)) );
  for(ix=0;ix<dimx;ix++) aux[ix]=gx[ix];
  gx[0]=0.;
  for(ix=1;ix<dimx;ix++)
    gx[ix]=gx[ix-1]+aux[ix-1];

  free(aux);

  return OK;
}


/****************************************************************************/
int gradient1D_FFT( int dimx, double *gx) {
/****************************************************************************/

  double *aux;
  double dR,x;
  int ix;

  TrackNullAlloc( aux=(double*)calloc(dimx,sizeof(double)) );
  
  Fourier1D(dimx,gx,aux,-1);
  
  for(ix=0;ix<dimx;ix++)    {
    
    x=((double)ix)/((double)dimx);
    if(ix>=dimx/2) x-=1.;
    x=sin(M_PI*x);
    dR=-x*aux[ix];
    aux[ix]=x*gx[ix];
    gx[ix]=dR;
  }
  
  Fourier1D(dimx,gx,aux,1);
  
  free(aux);
  return OK;
}


/****************************************************************************/
int reconstruct1D_FFT( int dimx, double *gx) {
/****************************************************************************/
 
  double *aux;
  double dR,x;
  int ix;

  TrackNullAlloc( aux=(double*)calloc(dimx,sizeof(double)) );

  Fourier1D(dimx,gx,aux,-1);

  for(ix=0;ix<dimx;ix++)    {
    
    x=((double)ix)/((double)dimx);
    if(ix>=dimx/2) x-=1.;
    if(x>1e-30) x=1./sin(M_PI*x);
    else x=0.;
    dR=x*aux[ix];
    aux[ix]=-x*gx[ix];
    gx[ix]=dR;
  }
  
  Fourier1D(dimx,gx,aux,1);
  
  free(aux);
  return OK;
}

/****************************************************************************/
int gradient1D_naifinter( int dimx, double *gx) {
/****************************************************************************/

  double *aux;
  double buff;
  int ix,ic;

  TrackNullAlloc( aux=(double*)calloc(dimx,sizeof(double)) );

  for(ix=0;ix<dimx;ix++) aux[ix]=gx[ix];
  gradient1D_naif(dimx,aux);

  for(ix=0;ix<dimx;ix++)    {
      for(ic=1,buff=aux[ix];(ic<10)&&(fabs(buff)<1e-3);ic++)
	buff+=aux[Mod(ix-ic,dimx)]+aux[Mod(ix+ic,dimx)];
      buff/=(double)(2*ic-1);
      gx[ix]=buff;
    }

  free(aux);
  return OK;
}


/****************************************************************************/
int filtro1D( int dimx,double expon, double *funcion) {
/****************************************************************************/
  double *funcR,*funcI;
  double x,f;
  double dx,dy;
  int ix;
  int xeff;

  xeff=redimensiona(dimx);
  TrackNullAlloc( funcR=(double*)calloc(xeff,sizeof(double)) );
  TrackNullAlloc( funcI=(double*)calloc(xeff,sizeof(double)) );
  
  copy1D(dimx,funcion,funcR,NULL);
  Fourier1D(xeff,funcR,funcI,-1);
  
  for(ix=0;ix<xeff;ix++)  {
    x=((double)ix)/((double)xeff);
    if(ix>=xeff/2) x-=1.;
    
    dx=2.*sin(M_PI*x);
    
    f=fabs(dx);
    if(f>1e-30) f=pow(f,expon);
    else f=0.;
    
    funcR[ix]=f*funcR[ix];
    funcI[ix]=f*funcI[ix];
  }

  Fourier1D(xeff,funcR,funcI,1);
  copy1D(dimx,funcR,funcion,NULL);

  free(funcR);

  return OK;
}

