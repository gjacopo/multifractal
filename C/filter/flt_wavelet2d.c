#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_alloc.h>
#include <utl_stats.h>
#include <utl_operator.h>

#include <filter.h>
#include <flt_wavelet1d.h>
#include <flt_wavelet2d.h>

#ifndef BUEN // criterion on the linear regression coefficient
#define BUEN 0.9
#endif

#ifdef _PARSE_FILTER_PARAMETERS_
#include <flt_parse.h>
extern ParWAV *p_wav;
#endif/* !_PARSE_FILTER_PARAMETERS_ */


/***************************************************************************/
void wavelet2D_genera( int dimx, int dimy, double sc0, double **wavesig) {
  /***************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  wavelet2D_define(dimx,dimy,sc0,p_wav->wav,p_wav->ord_der,wavesig);
  wavelet2D_normaliza(dimx,dimy,sc0,p_wav->ord_der,wavesig);

#else
  wavelet2D_define(dimx,dimy,sc0,WAV, ORDDER,wavesig);
  wavelet2D_normaliza(dimx,dimy,sc0,ORDDER,wavesig);

#endif/* !_PARSE_FILTER_PARAMETERS_ */
} // end of wavelet2D_genera


/***************************************************************************/
void wavelet2D_genera_line( int dimx, int dimy, double sc0, double **wave) {
/***************************************************************************/
 
#ifdef _PARSE_FILTER_PARAMETERS_
  wavelet2D_define_line(dimx,dimy,sc0,
			p_wav->thetau,p_wav->wav,p_wav->ord_der,wave);
  wavelet2D_normaliza(dimx,dimy,sc0,p_wav->ord_der,wave);

#else
  wavelet2D_define_line(dimx,dimy,sc0,THETA0,WAV,ORDDER,wave);
  wavelet2D_normaliza(dimx,dimy,sc0,ORDDER,wave);
#endif/* !_PARSE_FILTER_PARAMETERS_ */
} // end of wavelet2D_genera_line


/***************************************************************************/
double wavelet2D_escala( int dimx, int dimy, double sc0 ) {
  /***************************************************************************/
  double **wavesig;
  double sc;
  double norma;
  double x,y,r;
  int ix,iy;

  TrackNullAlloc( wavesig=matrix2D(dimy,dimx) );

  wavelet2D_define( dimx, dimy, sc0, 
#ifdef _PARSE_FILTER_PARAMETERS_
		    p_wav->wav, p_wav->ord_der,
#else
		   WAV, ORDDER,
#endif/* !_PARSE_FILTER_PARAMETERS_ */
		    wavesig );
  
  norma = media(dimx,dimy,wavesig,NULL);
  sc = sqrt(norma/M_PI);
  
  free_matrix2D(wavesig,dimy);
  
  return sc;
} // end of wavelet2D_escala


/***************************************************************************/
double wavelet2D_escala_line( int dimx, int dimy, double sc0 ) {
/***************************************************************************/
  double **wavesig;
  double x,y,mod,proy;
  double ux,uy;
  double sc=0.;
  int lmax=0;
  int ix,iy;
  
  TrackNullAlloc( wavesig=matrix2D(dimy,dimx) );
  
#ifdef _PARSE_FILTER_PARAMETERS_
  wavelet2D_define_line( dimx, dimy, sc0,
			p_wav->thetau, p_wav->wav, p_wav->ord_der, wavesig );
  ux = cos(p_wav->thetau), uy = sin(p_wav->thetau);
 
#else
  wavelet2D_define_line( dimx, dimy, sc0, THETA0, WAV, ORDDER, wavesig );
  ux = cos(THETA0), uy = sin(THETA0);
#endif/* !_PARSE_FILTER_PARAMETERS_ */

  for(iy=0;iy<dimy;iy++)    {
    y = (double)iy;
    if((iy>0)&&(iy>=dimy/2)) y -= (double)dimy;

    for(ix=0;ix<dimx;ix++)      {
      x = (double)ix;
      if((ix>0)&&(ix>=dimx/2)) x -= (double)dimx;
      mod = sqrt(x*x+y*y);

      if(mod<1e-30)	{
	lmax++;
	sc += 1.;
      } else {
	proy = fabs(uy*x-ux*y) / mod;
	if(proy < 1.e-10)	  {
	  lmax++;
	  sc += fabs(wavesig[iy][ix]);
	}
      }
    }
  }
  sc /= (double)lmax;
  
  free_matrix2D(wavesig,dimy);
  
  return sc;
} // end of wavelet2D_escala_line


/***************************************************************************/
int wt2d_transform( int dimx, int dimy, double a0, double quanto, 
		    double orden, int degree, int npunt,
		    double **muR, double **muI,
		    double **expon, double *minimo, double *maximo ) /*wavelet*/{
  /***************************************************************************/

  double ***Tmu,*x,*y;
  double scale,h,b,corr;
  int ia,ix,iy;
  int N=0;
  
  /* Useful initialisation */
  TrackNullAlloc( Tmu=matrix3D(npunt+1,dimy,dimx) );
  TrackNullAlloc( x=(double*)calloc(npunt+1,sizeof(double)) );
  TrackNullAlloc( y=(double*)calloc(npunt+1,sizeof(double)) );
  
  /*      Desarrollo      */
  scale=a0;
  for(ia=0;ia<=npunt;ia++)    {
    x[ia]=log(scale);
    convolverwave_log( dimx,dimy,muR,muI,scale,degree,orden,Tmu[ia] );
    scale=scale*quanto;
  }
  
  for(iy=0;iy<dimy;iy++)   
    for(ix=0;ix<dimx;ix++)      {
      for(ia=0;ia<=npunt;ia++) y[ia] = Tmu[ia][iy][ix];
      fit(x,y,npunt+1,&h,&b,&corr);
      expon[iy][ix] = h;
      if(fabs(corr) >= BUEN) N++;
      *maximo = fMax(*maximo,h);
      *minimo = fMin(*minimo,h);
    }
  
  /* Freeing memories */
  free_matrix3D(Tmu,npunt+1,dimy);
  free(x);
  free(y);

  return N;
} // end of wt2d_transform


/***************************************************************************/
int convolverwave_log( int dimx, int dimy, /* int xeff, int yeff, */
			double **muR, double **muI, double scale, int degree, 
			double orden, double **Tmu) {
  /***************************************************************************/
  
  double **psiR,**psiI;
  double buffR,buffI;
  int ix,iy;
  
  TrackNullAlloc( psiR=matrix2D(dimy,dimx) );
  TrackNullAlloc( psiI=matrix2D(dimy,dimx) );
  
  /*      Choice of the wavelet           */
  /*--------------------------------------*/
  if(orden < 4.)            lorentz2D(dimx,dimy,scale,degree,orden,psiR);
  else if(orden == 4)       gauss2D(dimx,dimy,scale,psiR);
  else if(orden == 5)       log2D_1(dimx,dimy,scale,psiR);
  else if(orden == 6)       log2D_2(dimx,dimy,scale,psiR);
  else if(orden == 7)       log2D_3(dimx,dimy,scale,psiR);
  else if(orden == 8)       test_wavelet(dimx,dimy,scale,orden,psiR);
  /*--------------------------------------*/

  /* Fourier transform of the wavelet */
  fill0(dimx,dimy,psiI,NULL);
  FFFT2D(dimx,dimy,psiR,psiI,-1);
  
  /* Cross product in the frequency space */
  for(iy=0;iy<dimy;iy++)   {
    for(ix=0;ix<dimx;ix++)	{
      buffR=psiR[iy][ix]*muR[iy][ix]-psiI[iy][ix]*muI[iy][ix];
      buffI=psiR[iy][ix]*muI[iy][ix]+psiI[iy][ix]*muR[iy][ix];
      psiR[iy][ix]=buffR;
      psiI[iy][ix]=buffI;
    }
  }
  
  /* Inverse Fourier transform to obtain the convolution */
  FFFT2D(dimx,dimy,psiR,psiI,1);
  
  for(iy=0;iy<dimy;iy++)    
    for(ix=0;ix<dimx;ix++)
      Tmu[iy][ix]=log(fabs(psiR[iy][ix]));
  
  free_matrix2D(psiR,dimy);
  free_matrix2D(psiI,dimy);

  return OK;
} // end of convolverwave_log


/***************************************************************************/
void test_wavelet( int dimx, int dimy, double a, double orden, double **psi) {
  /***************************************************************************/
  int ix,iy;
  double x,y,r;
  
  for(iy=0;iy<dimy;iy++)    {
    y=(double) iy;
    if(iy>=dimy/2) y-=(double)dimy;
    for(ix=0;ix<dimx;ix++)	{
      x=(double) ix;
      if(ix>=dimx/2) x-=(double)dimx;
      r=sqrt(x*x+y*y)/a;
      
      /*      Here put the Test Wavelet           */
      /* ----------------------------------------*/
      // Log_1
      psi[iy][ix] = log(2.0+r)/pow(1.0+r,4);   
      // Log_2
      // psi[iy][ix] = 1.0/( pow(1.0+r,3)*log(2.0+r) );        
      // Log_3
      // psi[iy][ix] = cos(0.01*r)*pow( log(2.0+r),3 )/pow(1.0+r,5);
      /* ----------------------------------------*/
      
      psi[iy][ix] = psi[iy][ix]/(a*a);
    }
  }
} // end of test_wavelet


/***************************************************************************/
void lorentz2D( int dimx, int dimy, double a, int degree, double orden,
	      double **psi) {
  /***************************************************************************/

  double x,y,r,pref;
  int ix,iy;
  
  for(iy=0;iy<dimy;iy++)    {
    y=(double) iy;
    if(iy>=dimy/2) y-=(double)dimy;
    for(ix=0;ix<dimx;ix++)      {
      x=(double) ix;
      if(ix>=dimx/2) x-=(double)dimx;
      r=sqrt(x*x+y*y)/a;
      switch(degree)	{
      case 0: pref=1.;
	break;
      case 1: pref=-r;
	break;
      case 2: pref=(2.*orden+1.)*r*r-1.;
	break;
      case 3: pref=-r*((2.*orden+1.)*r*r-3.);
	break;
      case 4: pref=
		(2.*orden+1.)*(2.*orden+3.)*r*r*r*r
		-6.*(2.*orden+3.)*r*r+3.;
      break;
      default: pref=1.;
      } 
      psi[iy][ix]=pref*pow(1.+r*r,-(orden+(double)degree))
	/(a*a);
    }
  }
  
} // end of lorentz2D


/***************************************************************************/
void gauss2D( int dimx, int dimy, double a, double **psi) {
  /***************************************************************************/
  
  int ix,iy;
  double x,y,r;
  
  for(iy=0;iy<dimy;iy++)    {
    y=(double) iy;
    if(iy>=dimy/2) y-=(double)dimy;
    for(ix=0;ix<dimx;ix++)      {
      x=(double) ix;
      if(ix>=dimx/2) x-=(double)dimx;

      if(fabs(a) > 1.e-17) r=sqrt(x*x+y*y)/a;
      psi[iy][ix] = exp(-0.5*r*r);
      if(fabs(a) > 1.e-17) psi[iy][ix] /= (a*a);
    }
  }
} // end of gauss2D


/***************************************************************************/
void log2D_1( int dimx, int dimy, double a, double **psi) {
  /***************************************************************************/
  int ix,iy;
  double x,y,r;

  for(iy=0;iy<dimy;iy++)    {
    y=(double) iy;
    if(iy>=dimy/2) y-=(double)dimy;
    for(ix=0;ix<dimx;ix++) {
	x=(double) ix;
	if(ix>=dimx/2) x-=(double)dimx;

	if(fabs(a) > 1.e-17) r=sqrt(x*x+y*y)/a;	
	psi[iy][ix] = log(2.0+r)/pow(1.0+r,4);
	if(fabs(a) > 1.e-17) psi[iy][ix] /= (a*a);
      }
  }
} // end of log2D_1


/***************************************************************************/
void log2D_2( int dimx, int dimy, double a, double **psi) {
  /***************************************************************************/
  int ix,iy;
  double x,y,r;
  
  for(iy=0;iy<dimy;iy++)    {
    y=(double) iy;
    if(iy>=dimy/2) y-=(double)dimy;
    for(ix=0;ix<dimx;ix++)      {
      x=(double) ix;
      if(ix>=dimx/2) x-=(double)dimx;

      if(fabs(a) > 1.e-17) r=sqrt(x*x+y*y)/a;      
      psi[iy][ix] = 1.0/( pow(1.0+r,3)*log(2.0+r) );
      if(fabs(a) > 1.e-17) psi[iy][ix] /= (a*a);
    }
  }
} // end of log2D_2


/***************************************************************************/
void log2D_3( int dimx, int dimy, double a, double **psi) {
  /***************************************************************************/
  int ix,iy;
  double x,y,r;
  
  for(iy=0;iy<dimy;iy++)	{
    y=(double) iy;
    if(iy>=dimy/2) y-=(double)dimy;
    for(ix=0;ix<dimx;ix++)      {
      x=(double) ix;
      if(ix>=dimx/2) x-=(double)dimx;

      if(fabs(a) > 1.e-17) r=sqrt(x*x+y*y)/a;      
      psi[iy][ix] = cos(0.01*r)*pow( log(2.0+r),3 )/pow(1.0+r,5);      
      if(fabs(a) > 1.e-17) psi[iy][ix] /= (a*a);
    }
  }
} // end of log2D_3


/***************************************************************************/
double escala2D_lineal( int dimx, int dimy) {
/***************************************************************************/
  return 1./sqrt((double)dimx*dimy);
} // end of escala2D_lineal


/***************************************************************************/
void wavelet2D_define_unit( int dimx, int dimy, double sc, int wav, int ord_der, 
			    double D0, double **wavesig )/*wavelet_2D*/ {
/***************************************************************************/

  double x,y;
  int ix,iy;

  /* first compute the wavelet at order=0 */
  for(iy=0;iy<dimy;iy++)    {
    y = ((double)iy)/((double)dimy);
    if(y > 0.5) y -= 1.;
    y *= sc*D0;

    for(ix=0;ix<dimx;ix++)	{
      x = ((double)ix)/((double)dimx);
      if(x > 0.5) x -= 1.;
      x *= sc*D0;
    
      if(wav == WAVMORLET) /* morlet */
	wavesig[iy][ix] = morlet2D_coeff(x,y);
      else if(wav == WAVGAUSS) /* gaussian */  
	wavesig[iy][ix] = gauss2D_coeff(x,y,0);
      else if(wav == WAVHAAR) /* haar */ 
	wavesig[iy][ix] = haar2D_coeff(x,y);
      else if(wav >= SWAVLORENTZ)  /* lorentzian */      
	wavesig[iy][ix] = lorentz2D_coeff(x,y,(double)wav/2.+.5,0); 
      /* Lorentzian of order 0.5 : wav=1 => exponent=1.
       * Lorentzian of order 1.  : wav=2 => exponent=1.5
       * Lorentzian of order 1.5 : wav=3 => exponent=2.
     */
    }
  }
  /* then derive: where ord_der intervenes */  
  if(ord_der) filtro2D(dimx,dimy,(double)ord_der,wavesig);

  anorma2(dimx,dimy,wavesig,NULL);
}


/***************************************************************************/
void wavelet2D_define( int dimx, int dimy, double sc, int wav, int ord_der, 
		       double **wavesig ) {
/***************************************************************************/
  double x,y,valor;
  int ix,iy;

  for(iy=0;iy<dimy;iy++)    {
    y=(double)iy;
    if((iy>0)&&(iy>=dimy/2)) y-=(double)dimy;
    y/=sc;

    for(ix=0;ix<dimx;ix++)      {
      x=(double)ix;
      if((ix>0)&&(ix>=dimx/2)) x-=(double)dimx;
      x/=sc;

      if(wav == WAVMORLET) /* morlet */
	wavesig[iy][ix] = morlet2D_coeff(x,y);
      else if(wav == WAVGAUSS) /* gaussian */  
	wavesig[iy][ix] = gauss2D_coeff(x,y,ord_der);
      else if(wav == WAVHAAR) /* haar */ 
	wavesig[iy][ix] = haar2D_coeff(x/(double)dimx,y/(double)dimy);
      else if(wav >= SWAVLORENTZ)  /* lorentzian */      
	wavesig[iy][ix] = lorentz2D_coeff(x,y,(double)wav/2.+.5,ord_der);
      /* Lorentzian of order 0.5 : wav=1 => exponent=1.
       * Lorentzian of order 1.  : wav=2 => exponent=1.5
       * Lorentzian of order 1.5 : wav=3 => exponent=2.
     */
    }
  }
} // end of wavelet2D_define


/***************************************************************************/
double haar2D_coeff( double x, double y ) {
  /***************************************************************************/

  if((fabs(x)>0.5) || (fabs(y)>0.5)) return 0.;
  else if(x*y>0)                     return 1.;
  else                               return -1.;
} // end of haar2D_coeff


/***************************************************************************/
double lorentz2D_coeff( double x, double y, double expon, int ord_der ) {
  /***************************************************************************
   * Called by : wavelet1D_define
   ***************************************************************************/
  double TT=(x*x+y*y),e2=2.*expon;

  switch (ord_der) {
  case 1:    return - e2 * sqrt(TT) * pow( 1. + TT, -expon - 1. );
  case 2:    return e2 * (TT *(e2+1.) - 1.) * pow( 1. + TT, -expon - 2. );
  default:
  case 0:    return pow( 1. + TT, -expon );
  }
} // end of lorentz2D_coeff


/***************************************************************************/
double morlet2D_coeff( double x, double y ) {
/***************************************************************************/
  double TT=(x*x+y*y);
  return exp(-.5*TT) * cos(5.*sqrt(TT));
} // end of morlet2D_coeff


/***************************************************************************/
double gauss2D_coeff( double x, double y, int ord_der ) {
  /***************************************************************************
   * Computes continuous Gaussian wavelet functions (0 to 5th derivative).
   * Called by : wavelet1D_define
   ***************************************************************************/
  double TT=(x*x+y*y), T; 
  
  switch (ord_der) {
  case 1:  return -sqrt(TT) * exp(-.5*TT);
  case 2:  return (TT-1.) * exp(-.5*TT);
  case 3:  return sqrt(TT) * exp(-.5*TT) * (3.-TT);
  case 4:  return exp(-.5*TT) * (pow(TT,2.) - 6.*TT + 3.);
  case 5:  T=sqrt(TT); return -T * exp(-.5*TT)* (pow(T,4.) - 10.*TT + 15.);
  default:
  case 0:  return exp(-.5*TT);
  }
} // end of gauss2D_coeff


/***************************************************************************/
void wavelet2D_define_line( int dimx, int dimy, double sc, double theta,
			    int wav, int ord_der, double **wavesig ) {
  /***************************************************************************/
  double x,y,mod,proy,valor;
  int ix,iy;
  double ux=cos(theta),uy=sin(theta);
  
  for(iy=0;iy<dimy;iy++)   {
    y = (double)iy;
    if((iy>0)&&(iy>=dimy/2)) y -= (double)dimy;
    y /= sc;

    for(ix=0;ix<dimx;ix++)	{
      x = (double)ix;
      if((ix>0)&&(ix>=dimx/2)) x -= (double)dimx;
      x /= sc;
      mod = sqrt(x*x+y*y);

      if(mod < 1e-30) 
	wavesig[iy][ix] = 1.;
      else	{
	proy = fabs(uy*x-ux*y) / mod;
	if(proy < 1e-10)	  {
	  if(wav == WAVMORLET) /* morlet */
	    wavesig[iy][ix] = morlet1D_coeff( mod );
	  else if(wav == WAVGAUSS) /* gaussian */  
	    wavesig[iy][ix] = gauss1D_coeff( mod, ord_der ); 	  
	  else if(wav == WAVHAAR) /* haar */ 
	    wavesig[iy][ix] = haar1D_coeff(mod/(double)dimx);
	  else if(wav >= SWAVLORENTZ)  /* lorentzian */      
	    wavesig[iy][ix] = lorentz1D_coeff( mod, (double)wav/2.+.5, ord_der );
	} else 
	  wavesig[iy][ix]=0.;
      }
    }
  }
  
} // end of wavelet2D_define_line


/***************************************************************************/
void wavelet2D_normaliza( int dimx, int dimy, double sc, int ord_der, 
			  double **wave) {
  /***************************************************************************/
  double norma=0.;
  int ix,iy;
  
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)
      norma += fabs(wave[iy][ix]);
  
  norma *= pow(sc,-ord_der);//  norma *= pow(sc,-ORDDER);
  
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)
      wave[iy][ix] /= norma;
  
} // end of wavelet2D_normaliza


