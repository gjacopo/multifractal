#include <stdio.h>
#include <math.h>

/* Personnal libraries */
#include <utils.h>
#include <utl_alloc.h>
#include <utl_stats.h>
#include <utl_operator.h>		

#include <filter.h>
#include <flt_fft1d.h>
#include <flt_fft2d.h>
#include <flt_stats1d.h>		
#include <flt_stats2d.h>		

#include <flt_deriva2d.h>


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


/*     Function declarations       */


/****************************************************************************/
int deriva2D_bak(int dimx, int dimy, double **gx, double **gy) {
  /****************************************************************************/
  /* Computes derivatives of an image by considering the frequency domain.
   * Signal in gxR is derivated. The components of the derivated signal are 
   * stored in gxR and gyR.
   * =========================================================================   
   * The derivative is omputed through the Fourier formula:
   * If FT stands for the Fourier transform and ' for the derivation, we have
   * for any integrable function f:
   *            FT{f'}(y) = iy FT{f}(y)
   * where y is the frequency and i the imaginary unit, so the derivative of f
   * can be expressed wrt its Fourier transform FT{f} by:
   *            f'(x) = IFT{ y->iy FT{f}(y)}
   * where IFT stands for the inverse Fourier transform. */
  /****************************************************************************/

  int ix,iy;
  double x,y,aux,auxI;
  double **gxI, **gyI;
  
  TrackNullAlloc( gxI=matrix2D(dimy,dimx) );
  TrackNullAlloc( gyI=matrix2D(dimy,dimx) );
  /*   fill0(dimx,dimy,gxI, NULL);	
       fill0(dimx,dimy,gyI, NULL);    */
  
  /* Direct Fourier Transform of the signal in gx */
  Fourier2D(dimx,dimy,gx,gxI,-1);
 
  for(iy=0;iy<dimy;iy++)  {
    y=((double)iy)/((double) dimy);
    if(iy>=dimy/2) y-=1.;
    if(iy==dimy/2) y=0.;

    for(ix=0;ix<dimx;ix++)	{
      x=((double)ix)/((double) dimx);
      if(ix>=dimx/2) x-=1.;
      if(ix==dimx/2) x=0.;
      
      /* Vecteur frequence */
#ifdef _PARSE_FILTER_PARAMETERS_
      if(p_fil->mode_freq == FREQSIN)
#else
	if(MODE_FREQ == FREQSIN)
#endif/* !_PARSE_FILTER_PARAMETERS_ */
	x=sin(M_PI*x), y=sin(M_PI*y);
      
      aux = gx[iy][ix];
      auxI = gxI[iy][ix];
      
      /* Multiplication by if where f=(dx,dy) is the frequency vector and i
       * is the imaginary unit */
      gx[iy][ix] = -x*auxI; 
      gxI[iy][ix] = x*aux;
      gy[iy][ix] = -y*auxI;
      gyI[iy][ix] = y*aux;
    }
  }
  /* Inverse Fourier transform of the signals so obtain */
  Fourier2D(dimx,dimy,gx,gxI,1);
  Fourier2D(dimx,dimy,gy,gyI,1);
  
  free_matrix2D(gxI,dimy);
  free_matrix2D(gyI,dimy);

  return OK;
} // end of deriva2D_bak


/***************************************************************************/
int modderiva_naif( int dimx, int dimy, double **mu) {
  /***************************************************************************/
  /* Computes the norm of a naive derivatives */
  /***************************************************************************/

  double *linea_derecha,*linea_abajo;
  double dx,dy,esquina;
  int ix,iy;
  
  TrackNullAlloc( linea_derecha=(double*)calloc(dimy-1,sizeof(double)) );
  TrackNullAlloc( linea_abajo=(double*)calloc(dimx-1,sizeof(double)) );

  for(iy=0;iy<dimy-1;iy++)    {
    dx=mu[iy][0]-mu[iy][dimx-1];
    dy=mu[iy+1][dimx-1]-mu[iy][dimx-1];
    linea_derecha[iy]=sqrt(dx*dx+dy*dy);
  }
  for(ix=0;ix<dimx-1;ix++)  {
    dx=mu[dimy-1][ix+1]-mu[dimy-1][ix];
    dy=mu[0][ix]-mu[dimy-1][ix];
    linea_abajo[ix]=sqrt(dx*dx+dy*dy);
  }
  dx=mu[dimy-1][0]-mu[dimy-1][dimx-1];
  dy=mu[0][dimx-1]-mu[dimy-1][dimx-1];
  esquina=sqrt(dx*dx+dy*dy);
  
  for(iy=0;iy<dimy-1;iy++)    
    for(ix=0;ix<dimx-1;ix++)	{
      dx=mu[iy][ix+1]-mu[iy][ix];
      dy=mu[iy+1][ix]-mu[iy][ix];
      mu[iy][ix]=sqrt(dx*dx+dy*dy);
    }
  
  for(iy=0;iy<dimy-1;iy++) mu[iy][dimx-1]=linea_derecha[iy];
  for(ix=0;ix<dimx-1;ix++) mu[dimy-1][ix]=linea_abajo[ix];
  mu[dimy-1][dimx-1]=esquina;
  
  free(linea_derecha);
  free(linea_abajo);

  return OK;
} // end of modderiva_naif

/***************************************************************************/
int modderiva( int dimx, int dimy, double **modg) {
/***************************************************************************/
  double **aux;
  int ix,iy;

  TrackNullAlloc( aux=matrix2D(dimy,dimx) );
  
  gradient2D(dimx,dimy,modg,aux);
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)
      modg[iy][ix]=sqrt(modg[iy][ix]*modg[iy][ix]+aux[iy][ix]*aux[iy][ix]); 
  
  free_matrix2D(aux,dimy);

  return OK;
} // end of modderiva


/***************************************************************************/
int modderiva_line( int dimx, int dimy, double **modg, double theta ) {
  /***************************************************************************/
  double **aux;
  double ux,uy;
  int ix,iy;

  ux=cos(theta), uy=sin(theta);

  TrackNullAlloc( aux=matrix2D(dimy,dimx) );
  
  gradient2D(dimx,dimy,modg,aux);
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)
      modg[iy][ix]=fabs(ux*modg[iy][ix]+uy*aux[iy][ix]); 
  
  free_matrix2D(aux,dimy);

  return OK;
} // end of modderiva_line


/***************************************************************************/
int modgradderiva( int dimx, int dimy, 
		    double **mu, double **gx, double **gy) {
  /***************************************************************************/
  /* Computes the norm of a naive derivatives */
  /***************************************************************************/

  double *linea_derecha,*linea_abajo;
  double dx,dy,esquina;
  int ix,iy;
  
  TrackNullAlloc( linea_derecha=(double*)calloc(dimy-1,sizeof(double)) );
  TrackNullAlloc( linea_abajo=(double*)calloc(dimx-1,sizeof(double)) );

  for(iy=0;iy<dimy-1;iy++)    {
    gx[iy][dimx-1] = dx=mu[iy][0]-mu[iy][dimx-1];
    gy[iy][dimx-1] = dy=mu[iy+1][dimx-1]-mu[iy][dimx-1];
    linea_derecha[iy]=sqrt(dx*dx+dy*dy);
  }
  for(ix=0;ix<dimx-1;ix++)  {
    gx[dimy-1][ix] = dx=mu[dimy-1][ix+1]-mu[dimy-1][ix];
    gy[dimy-1][ix] = dy=mu[0][ix]-mu[dimy-1][ix];
    linea_abajo[ix]=sqrt(dx*dx+dy*dy);
  }
  gx[dimy-1][dimx-1] = dx=mu[dimy-1][0]-mu[dimy-1][dimx-1];
  gy[dimy-1][dimx-1] = dy=mu[0][dimx-1]-mu[dimy-1][dimx-1];
  esquina=sqrt(dx*dx+dy*dy);
  
  for(iy=0;iy<dimy-1;iy++)    
    for(ix=0;ix<dimx-1;ix++)	{
      gx[iy][ix] = dx=mu[iy][ix+1]-mu[iy][ix];
      gy[iy][ix] = dy=mu[iy+1][ix]-mu[iy][ix];
      mu[iy][ix]=sqrt(dx*dx+dy*dy);
    }
  
  for(iy=0;iy<dimy-1;iy++) mu[iy][dimx-1]=linea_derecha[iy];
  for(ix=0;ix<dimx-1;ix++) mu[dimy-1][ix]=linea_abajo[ix];
  mu[dimy-1][dimx-1]=esquina;
  
  free(linea_derecha);
  free(linea_abajo);

  return OK;
} // end of modgradderiva


/****************************************************************************/
void gradient2D( int dimx, int dimy, double **gx, double **gy) {
/****************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  switch(p_fil->mode_deriva)   
#else
    switch(MODE_DERIVA)  
#endif/* !_PARSE_FILTER_PARAMETERS_ */
      {
      case 0:
	gradient2D_FFT(dimx,dimy,gx,gy);
	break;
      case 1:
	gradient2D_naif(dimx,dimy,gx,gy);
	break;
      case 2:
	gradient2D_naif_inter(dimx,dimy,gx,gy);
	break;
      default:
	Exit("Unrecognized derivative mode");
      }
} // end of gradient2D


/****************************************************************************/
void reconstruct2D( int dimx, int dimy, double **gx, double **gy) {
/****************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  switch(p_fil->mode_deriva)    
#else
    switch(MODE_DERIVA) 
#endif/* !_PARSE_FILTER_PARAMETERS_ */
      {
      case 0:
	reconstruct2D_FFT(dimx,dimy,gx,gy);
	break;
      case 1:
      case 2:
	reconstruct2D_naif(dimx,dimy,gx,gy);
	break;
      default:
	Exit("Unrecognized derivative mode");
      }
} // end of reconstruct2D


/****************************************************************************/
int gradient2D_naif_bak( int dimx, int dimy, double **gx, double **gy) {
/****************************************************************************/

  double *columna_derecha;
  int ix,iy;

  TrackNullAlloc( columna_derecha=(double*)calloc(dimy,sizeof(double)) );

  for(iy=0;iy<dimy;iy++) columna_derecha[iy]=gx[iy][0]-gx[iy][dimx-1];

  for(iy=0;iy<dimy-1;iy++) 
    for(ix=0;ix<dimx;ix++) gy[iy][ix]=gx[iy+1][ix]-gx[iy][ix];
  for(ix=0;ix<dimx;ix++) gy[dimy-1][ix]=gx[0][ix]-gx[dimy-1][ix];
  
  for(ix=0;ix<dimx-1;ix++)
    for(iy=0;iy<dimy;iy++) gx[iy][ix]=gx[iy][ix+1]-gx[iy][ix];
  for(iy=0;iy<dimy;iy++) gx[iy][dimx-1]=columna_derecha[iy];
  
  free(columna_derecha); 

  return OK;
} // end of gradient2D_naif_bak


/****************************************************************************/
int gradient2D_naif( int dimx, int dimy, double **gx, double **gy) {
/****************************************************************************/

  double **gxI,**gyI;
  double dxR,dxI,dyR,dyI;
  double aux,x,y;
  int ix,iy;

  TrackNullAlloc( gxI=matrix2D(dimy,dimx) );
  TrackNullAlloc( gyI=matrix2D(dimy,dimx) );

  /* Direct Fourier Transform of the signal in gx */
  Fourier2D(dimx,dimy,gx,gxI,-1);

  for(iy=0;iy<dimy;iy++)    {
    y = ((double)iy) / ((double)dimy);
    dyR = cos(2.*M_PI*y)-1.;
    dyI = sin(2.*M_PI*y);

    for(ix=0;ix<dimx;ix++)      {
      x = ((double)ix) / ((double)dimx);
      dxR = cos(2.*M_PI*x)-1.;
      dxI = sin(2.*M_PI*x);
      
      gy[iy][ix] = dyR*gx[iy][ix] - dyI*gxI[iy][ix];
      gyI[iy][ix] = dyR*gxI[iy][ix] + dyI*gx[iy][ix];
      aux = dxR*gx[iy][ix] - dxI*gxI[iy][ix];
      gxI[iy][ix] = dxR*gxI[iy][ix] + dxI*gx[iy][ix];
      gx[iy][ix] = aux;
    }  
  }
  
  Fourier2D(dimx,dimy,gx,gxI,1);
  Fourier2D(dimx,dimy,gy,gyI,1);

  free_matrix2D(gxI,dimy);
  free_matrix2D(gyI,dimy);

  return OK;
} // end of gradient2D_naif


/****************************************************************************/
int reconstruct2D_naif( int dimx, int dimy, double **gx, double **gy) {
  /****************************************************************************/

  double **gxI,**gyI;
  double dxR,dxI,dyR,dyI,modx,mody;
  double auxR,auxI,x,y;
  int ix,iy;

  TrackNullAlloc( gxI=matrix2D(dimy,dimx) );
  TrackNullAlloc( gyI=matrix2D(dimy,dimx) );

  Fourier2D(dimx,dimy,gx,gxI,-1);
  Fourier2D(dimx,dimy,gy,gyI,-1);

  for(iy=0;iy<dimy;iy++)    {
    y=((double)iy)/((double)dimy);
    dyR = cos(2.*M_PI*y)-1.;
    dyI = -sin(2.*M_PI*y);
    mody=dyR*dyR+dyI*dyI;

    for(ix=0;ix<dimx;ix++)      {
      x=((double)ix)/((double)dimx);
      dxR = cos(2.*M_PI*x)-1.;
      dxI = -sin(2.*M_PI*x);
      modx=dxR*dxR+dxI*dxI;
      
      if(modx+mody > 1.e-30)	{
	auxR = (dxR*gx[iy][ix]-dxI*gxI[iy][ix]
	      +dyR*gy[iy][ix]-dyI*gyI[iy][ix]) / (modx+mody);
	auxI = (dxR*gxI[iy][ix]+dxI*gx[iy][ix]
	      +dyR*gyI[iy][ix]+dyI*gy[iy][ix]) / (modx+mody);
      }      else	    {
	auxR = 0.;
	auxI = 0.;
      }
      gx[iy][ix] = auxR;
      gxI[iy][ix] = auxI;
    }  
  }

  Fourier2D(dimx,dimy,gx,gxI,1);

  free_matrix2D(gxI,dimy);
  free_matrix2D(gyI,dimy);

  return OK;
} // end of reconstruct2D_naif


/****************************************************************************/
int gradient2D_FFT( int dimx, int dimy, double **gx, double **gy) {
/****************************************************************************/

  double **gxI,**gyI;
  double aux,x,y;
  int ix,iy;

  TrackNullAlloc( gxI=matrix2D(dimy,dimx) );
  TrackNullAlloc( gyI=matrix2D(dimy,dimx) );

  Fourier2D(dimx,dimy,gx,gxI,-1);

  for(iy=0;iy<dimy;iy++)    {
    y=((double)iy)/((double)dimy);
    if((iy>0)&&(iy>=dimy/2)) y-=1.;
    y=sin(M_PI*y);
    for(ix=0;ix<dimx;ix++)      {
      x=((double)ix)/((double)dimx);
      if((ix>0)&&(ix>=dimx/2)) x-=1.;
      x=sin(M_PI*x);
      
      gy[iy][ix]=-y*gxI[iy][ix];
      gyI[iy][ix]=y*gx[iy][ix];
      aux=-x*gxI[iy][ix];
      gxI[iy][ix]=x*gx[iy][ix];
      gx[iy][ix]=aux;
    }
    
  }
  
  Fourier2D(dimx,dimy,gx,gxI,1);
  Fourier2D(dimx,dimy,gy,gyI,1);

  free_matrix2D(gxI,dimy);
  free_matrix2D(gyI,dimy);

  return OK;
} // end of gradient2D_FFT

/****************************************************************************/
int reconstruct2D_FFT( int dimx, int dimy, double **gx, double **gy) {
  /***************************************************************************
   * Reconstruction of a signal c from the essential (vector-valued density) 
   * and the universal propagator g:
   *        c = g ** v
   * where ** stands for the convolution product (being understood as a 
   * scalar products of vectors).
   ***************************************************************************/

  double **gxI,**gyI;
  double aux,x,y,f;
  int ix,iy;

  TrackNullAlloc( gxI=matrix2D(dimy,dimx) );
  TrackNullAlloc( gyI=matrix2D(dimy,dimx) );

  /* ((gx,gxI) (gy,gyI)):  Fourier transform of (gx,gy) 
   * ie., complex bidimensionnal vector field */
  Fourier2D(dimx,dimy,gx,gxI,-1);
  Fourier2D(dimx,dimy,gy,gyI,-1);

  for(iy=0;iy<dimy;iy++)    {
    y=((double)iy)/((double)dimy);
    if((iy>0)&&(iy>=dimy/2)) y-=1.;
    /* if(iy > dimy/2) y -= 1.;
       else if(iy == dimy/2) y = 0.; */

    for(ix=0;ix<dimx;ix++)	{
      x=((double)ix)/((double)dimx);
      if((ix>0)&&(ix>=dimx/2)) x-=1.;
      /* if(ix > dimx/2) x -= 1.;
	 else if(ix == dimx/2) x = 0.; */
      
#ifdef _PARSE_FILTER_PARAMETERS_
      if(p_fil->mode_freq == FREQSIN)
#else
	if(MODE_FREQ == FREQSIN)
#endif/* !_PARSE_FILTER_PARAMETERS_ */
	  x=sin(M_PI*x), y=sin(M_PI*y);

      f = x*x + y*y;
      /* The scalar product with the frequency vector f is computed
       * and, simultaneously, the vector so otained is rotated.
       * This corresponds to a multiplication by the imaginary unit
       * i=sqrt(-1) of the vector viewed as a complex number:
       *     A = a_0 + a_1*i   =>   B = A*i = -a_1 +a_0*i          */
      if(f > 1.e-30)	{
	aux = (x*gxI[iy][ix] + y*gyI[iy][ix]) / f;
	gxI[iy][ix] = -(x*gx[iy][ix] + y*gy[iy][ix]) / f;
	gx[iy][ix] = aux;
      }	  else 
	gx[iy][ix] = gxI[iy][ix] = 0.;
      
      gy[iy][ix] = gyI[iy][ix] = 0.;
    }    
  }

  /* Inverse Fourier transform */
  Fourier2D(dimx,dimy,gx,gxI,1);

  free_matrix2D(gxI,dimy);
  free_matrix2D(gyI,dimy);

  return OK;
} // end of reconstruct2D_FFT


/****************************************************************************/
int gradient2D_naif_inter( int dimx, int dimy, double **gx, double **gy) {
/****************************************************************************/
  double x,y;
  int ix,iy;

  for(iy=0;iy<dimy;iy++)    {
    y=(double)iy;
    if(iy>dimy/2) y-=(double)dimy;
    for(ix=0;ix<dimx;ix++)      {
      x=(double)ix;
      if(ix>dimx/2) x-=(double)dimx;
      gy[iy][ix]=1./(.00001+x*x+y*y);
    }
  }

  convuelto2D(dimx,dimy,gy,gx);
  gradient2D_naif(dimx,dimy,gx,gy);

  return OK;
} // end of gradient2D_naif_inter


/****************************************************************************/
void gradient_complex( int dimx, int dimy, double **gx, double **gy) {
/****************************************************************************/

  fill0(dimx,dimy,gy,NULL);
  deriva_complex(dimx,dimy,0,gx,gy);
} // end of gradient_complex


/****************************************************************************/
void reconstruct_complex( int dimx, int dimy, double **gx, double **gy) {
/****************************************************************************/

  deriva_complex(dimx,dimy,1,gx,gy);
  fill0(dimx,dimy,gy,NULL);
} // end of reconstruct_complex


/****************************************************************************/
void deriva_complex( int dimx, int dimy, int mode, double **gx, double **gy) {
/****************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  switch(p_fil->mode_deriva)   
#else
    switch(MODE_DERIVA) 
#endif/* !_PARSE_FILTER_PARAMETERS_ */
    {
    case 0:
      deriva_complex_FFT(dimx,dimy,mode,gx,gy);
      break;
    case 1:
      deriva_complex_naif(dimx,dimy,mode,gx,gy);
      break;
    default:
      printf("Unrecognized derivative mode\n");
      exit(ERROR);
    }
} // end of deriva_complex
 
 
/****************************************************************************/
 void deriva_complex_naif( int dimx, int dimy, int mode, double **gx, 
			   double **gy) {
/****************************************************************************/

  double dxR,dxI,dyR,dyI;
  double dR,dI,modd;
  double aux,x,y;
  int ix,iy;

  Fourier2D(dimx,dimy,gx,gy,-1);

  for(iy=0;iy<dimy;iy++)    {
    y=((double)iy)/((double)dimy);
    dyR=-(cos(2.*M_PI*y)-1.);
    dyI=-sin(2.*M_PI*y);
    
    for(ix=0;ix<dimx;ix++)	{
      x=((double)ix)/((double)dimx);
      dxR=cos(2.*M_PI*x)-1.;
      dxI=sin(2.*M_PI*x);
      
      dR=dxR+dyI;
      dI=dxI-dyR;
      modd=dR*dR+dI*dI;
      if(mode)	{
	if(modd>1e-30)	    {
	  dR=dR/modd;
	  dI=-dI/modd;
	}	else	  {
	  dR=0.;
	  dI=0;
	}
      }
      
      aux=dR*gx[iy][ix]-dI*gy[iy][ix];
      gy[iy][ix]=dR*gy[iy][ix]+dI*gx[iy][ix];
      gx[iy][ix]=aux;
    }
  }
  
  Fourier2D(dimx,dimy,gx,gy,1);
} // end of deriva_complex_naif


/****************************************************************************/
void deriva_complex_FFT( int dimx, int dimy, int mode, double **gx, 
			 double **gy) {
/****************************************************************************/

  double dxR,dxI,dyR,dyI;
  double dR,dI,modd;
  double aux,x,y;
  int ix,iy;

  Fourier2D(dimx,dimy,gx,gy,-1);
  
  /* Se va a multiplicar por el numero complejo dx+i*conj*dy; de este modo, 
   * si la funcion fuera real se obtendria un numero complejo cuyas componentes 
   * real e imaginaria corresponderian, respectivamente, a las componentes x e 
   * y del gradiente
   * Let us notice that for some kernels, dx and dy are purely imaginary
   */
  for(iy=0;iy<dimy;iy++)    {
    y = ((double)iy)/((double)dimy);
    if(iy > dimy/2) y -= 1.;
    dyR = 0.;
    dyI = -sin(M_PI*y);
    
    for(ix=0;ix<dimx;ix++)	{
      x = ((double)ix)/((double)dimx);
      if(ix > dimx/2) x -= 1.;
      dxR = 0.;
      dxI = sin(M_PI*x);
      
      dR = dxR + dyI;
      dI = dxI - dyR;
      /*     conj = (mode==2) ? 1. : -1.;
       * Definition a medio pixel :
       *     dI = -2. * sin(M_PI*x);
       *     dR = conj * 2. * sin(M_PI*y);
       * Definition a un pixel      
       *    dR = -2.*sin(PI*xs)*sin(PI*xs) - 2.*conj*sin(PI*ys)*cos(PI*ys);
       *    dI = 2.*sin(PI*xs)*cos(PI*xs) - 2.*conj*sin(PI*ys)*sin(PI*ys);
      */
      modd = dR*dR + dI*dI;

      if(mode)	    
	if(modd > 1e-30)
	  dR = dR/modd,  dI = -dI/modd;
       	else	  
	  dR = dI = 0.;
      
      aux = dR*gx[iy][ix] - dI*gy[iy][ix];
      gy[iy][ix] = dR*gy[iy][ix] + dI*gx[iy][ix];
      gx[iy][ix] = aux;
    }
  }
  
  Fourier2D(dimx,dimy,gx,gy,1);
} // end of deriva_complex_FFT


/***************************************************************************/
 int filtro2D_bak( int dimx, int dimy, double norma, double expon,
		   double **data) {
  /***************************************************************************
   * Transformation in the frequency domain: computes the value of a matrice 
   * after having processed a scalar product with the norm of the frequency 
   * vector in the frequency domain.
   * =========================================================================
   * Computes:
   *      IFT( FT(func) * f^expon) )
   * namely:
   *       / FT(data) * f^expon   for non-nul frequencies,
   *       \ FT(data) * norma     for DC-component,
   * where f=(2sin(pi.x), 2sin(pi.y)) stands for the norm of the 
   * frequency vector, and eventually puts the DC-component to 0,
   * and then computes the IFT.
   *
   * Parse:
   *     - data: matrice to be transformed in the frequency 
   *       domain,
   *     - norma: coefficient multiplication of the DC-component 
   *       (0-frequency),
   *     - expon: exponent of the power of f wich will be used
   *       in the multiplication of the FT of data, where f
   *       stands for the frequency vector.
   * Filtrage  dans l'espace des frequences:
   *    - multiplication par norma pour la composante de frequence nulle,
   *    - multiplication par f^expon, ou f designe la norme du vecteur 
   *      frequence, pour les autres composantes.       
  ***************************************************************************/
  int ix,iy;
  double x,y,f,output;
  double **funcI;
  
  TrackNullAlloc( funcI=matrix2D(dimy,dimx) );
  /* fill0(dimx,dimy,funcI, NULL); */
  
  /* Direct Fourier Transform */
  Fourier2D(dimx,dimy,data,funcI,-1);
  
  for(iy=0;iy<dimy;iy++)  {
    y=((double)iy)/((double)dimy);
    if(iy>=dimy/2) y-=1.;
    
    for(ix=0;ix<dimx;ix++)	{
      x=((double)ix)/((double)dimx);
      if(ix>=dimx/2) x-=1.;
      
      /* Vecteur frequence */
#ifdef _PARSE_FILTER_PARAMETERS_
      if(p_fil->mode_freq == FREQSIN)
#else
	if(MODE_FREQ == FREQSIN)
#endif/* !_PARSE_FILTER_PARAMETERS_ */
	x=2.*sin(M_PI*x), y=2.*sin(M_PI*y);
 
      f = sqrt(x*x+y*y);
      if(f>1e-30) f=pow(f,expon);
      else f=0.;
      
      if((ix==0)&&(iy==0))	{
	data[iy][ix] *= norma;
	funcI[iy][ix] *= norma;
      }      else	{
	data[iy][ix] *= f;
	funcI[iy][ix] *= f;
      }
    }
  }
  
  /* Inverse Fourier Transform */
  Fourier2D(dimx,dimy,data,funcI,1);
  
  free_matrix2D(funcI,dimy);

  return OK;
} // end of filtro2D_bak


/****************************************************************************/
int filtro2D( int dimx, int dimy, double expon, double **data) {
  /***************************************************************************
   * Transformation in the frequency domain: computes the value of a matrice 
   * after having processed a scalar product with the norm of the frequency 
   * vector in the frequency domain.
   ***************************************************************************/

  double **auxR,**auxI;
  double x,y,f;
  int xeff,yeff;
  int ix,iy;

  xeff=redimensiona(dimx),yeff=redimensiona(dimy);
  
  TrackNullAlloc( auxR=matrix2D(yeff,xeff) );
  TrackNullAlloc( auxI=matrix2D(yeff,xeff) );

  copy(dimx,dimy,data,auxR,NULL);
  
  Fourier2D(xeff,yeff,auxR,auxI,-1);
  
  for(iy=0;iy<yeff;iy++)    {
    y=((double)iy)/((double)yeff);
    if(iy>=yeff/2) y-=1.;

    for(ix=0;ix<xeff;ix++)      {
      x=((double)ix)/((double)xeff);
      if(ix>=xeff/2) x-=1.;

       /* Vecteur frequence */
#ifdef _PARSE_FILTER_PARAMETERS_
      if(p_fil->mode_freq == FREQSIN)
#else
	if(MODE_FREQ == FREQSIN)
#endif/* !_PARSE_FILTER_PARAMETERS_ */
	  x=sin(M_PI*x), y=sin(M_PI*y);
      
      f=sqrt(x*x+y*y);

      if(f>1e-30) f=pow(f,expon);
      else f=0.;
      auxR[iy][ix] *= f;
      auxI[iy][ix] *= f;
    }
  }
  
  Fourier2D(xeff,yeff,auxR,auxI,1);
  
  copy(dimx,dimy,auxR,data,NULL);

  free_matrix2D(auxR,yeff);
  free_matrix2D(auxI,yeff);
  
  return OK;
} // end of filtro2D

