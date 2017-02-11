#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_alloc.h>
#include <utl_stats.h>
#include <utl_operator.h>

#include <filter.h>
#include <flt_wavelet1d.h>


#ifdef _PARSE_FILTER_PARAMETERS_
#include <flt_parse.h>
extern ParFILT *p_fil;
extern ParWAV *p_wav;
#ifndef redimensiona
#define redimensiona(x) (p_fil->flag_memory?(x):dimensiona(x))
#endif

#else
#ifndef redimensiona
#define redimensiona(x) (FLAG_MEMORY?(x):dimensiona(x))
#endif

#endif/* !_PARSE_FILTER_PARAMETERS_ */



/***************************************************************************/
void wavelet1D_genera( int dimx, double sc0 double *wavesig ) {
  /***************************************************************************/
  
#ifdef _PARSE_FILTER_PARAMETERS_
  wavelet1D_define(dimx, sc0, p_wav->wav,p_wav->ord_der,wavesig);
  wavelet1D_normaliza(dimx, sc0,p_wav->ord_der,wavesig);

#else
  wavelet1D_define(dimx,sc0,WAV, ORDDER,wavesig);
  wavelet1D_normaliza(dimx,sc0,ORDDER,wavesig);

#endif/* !_PARSE_FILTER_PARAMETERS_ */
		      
} // end of wavelet1D_genera


/***************************************************************************/
double wavelet1D_escala( int dimx, double sc0 ) {
 /***************************************************************************/
  double *wavesig;
  double sc;
  double norma;
  int xeff;
  int ix;

  xeff=redimensiona(dimx);
  TrackNullAlloc( wavesig=(double*)calloc(xeff,sizeof(double)) );

  wavelet1D_define( xeff, sc0,
#ifdef _PARSE_FILTER_PARAMETERS_
		   p_wav->wav,p_wav->ord_der,
#else
		   SC0, WAV, ORDDER,
#endif/* !_PARSE_FILTER_PARAMETERS_ */
		   wavesig );

  for(ix=0,norma=0.;ix<xeff;ix++) norma+=wavesig[ix]*wavesig[ix];
  
  // Like that, norma is assumed to be a diameter. It follows
  sc = sqrt(norma/dimx);

  free(wavesig);

  return sc;
} // end of wavelet1D_escala


/***************************************************************************/
int wt1D_transform( double *signal, int dimx, int wav, int order,
		    double sc, double *wtrans ){
  /***************************************************************************
   *
   * Compute the wavelet transform of the signal at scale sc
   *
   * Parse:
   *    - signal : original signal,
   *    - dimx : size of the signal,
   *    - wav : variable used for the choice of the wavelet,
   *    - sc : current scale of analysis,
   *    - maxfilt : maximal size of the wavelet filter.
   * Output:
   *    - wtrans : 1d tabular storing the values of the WT at scale sc,
   *
   * Note : signal is zero padded to remove aliasing. FFT is used
   * Called by : wtmm_pf_compute  
   ***************************************************************************/
 
  int i, j;
  int tempi;
  double *wavesig,*auxsig,*auxwtrans;
 
 /* allocate memory for the wavelet and extended signal and wavele transform */
  TrackNullAlloc( wavesig=(double*)calloc(2*dimx,sizeof(double)) );
  TrackNullAlloc( auxsig=(double*)calloc(2*dimx,sizeof(double)) );
  TrackNullAlloc( auxwtrans=(double*)calloc(2*dimx,sizeof(double)) ); 
  
  copy1D(dimx,signal,auxsig,NULL); // auxsig is zero-padded
  
  /* Create the wavelet filter */
  wavelet1D_define( 2*dimx, sc, wav, order, wavesig ); 
  
  /* Obtaining the wavelet transforms   */
  wt1D_projection( auxsig, wavesig, 2*dimx, sc, auxwtrans );
  
  /*  We just keep the first half */
  copy1D(dimx,auxwtrans,wtrans,NULL);
      
  if(wavesig)  free(wavesig);
  if(auxsig)  free(auxsig);
  if(auxwtrans)  free(auxwtrans);
  
  return OK;
} // end of wt1D_transform


/***************************************************************************/
void wt1D_projection( double *signal,  double *filt, int dimx, double sc, 
		      double *wtrans )  {
  /***************************************************************************
   * Naive convolution + normalization to compute the wavelet coefficient
   *
   * Parse:
   *    - signal : original signal,
   *    - dimx : size of the signal, the filter and the wavelet transform
   *    - sc : current scale of analysis,
   *    - filt : wavelet filter used to convolve the signal at scale sc,
   * Output:
   *    - wtrans : values of the WT at scale sc,
   *
   * Called by : wt1D_transform
   ***************************************************************************/
  int i;
  
  convuelve1D(dimx,signal,filt,wtrans);
  /* absolute value of wavelet coefficient */
  for(i=0;i<dimx;i++) wtrans[i] = fabs(wtrans[i]);  
} // end of wt1D_projection


/***************************************************************************/
void wavelet1D_define_unit( int dimx, double sc, int wav, int ord_der, 
			    double D0, double *wavesig )/*wavelet_1D*/ {
  /***************************************************************************/  
  double x;
  int ix,ix0,ix1;
  
  for(ix=0;ix<dimx;ix++)    {
    x = ((double)ix)/((double)dimx);
    if(x > 0.5) x -= 1.;
    x *= sc*D0;
    
    if(wav == WAVGAUSS) /* gaussian */  
      wavesig[ix] = gauss1D_coeff(x,ord_der); 
    else if(wav == WAVHAAR) /* haar */ 
      wavesig[ix] = haar1D_coeff(x);
    else if(wav >= SWAVLORENTZ)  /* lorentzian */      
      wavesig[ix] = lorentz1D_coeff(x,(double)wav,ord_der); 
  }
  
  if(wav == WAVHAAR)    {
    ix0 = Mod((int)( ((double)dimx)/(2.*sc*D0))-1.,dimx);
    ix1 = Mod(dimx-1-ix0,dimx);
    if(ord_der == 1)      {
      wavesig[0]=1.;
      wavesig[ix0]=wavesig[ix1]=-0.5;
    }    else if(ord_der == 2)      {
      wavesig[0]=-1.;
      wavesig[dimx-1]=1.;
      wavesig[ix0]=wavesig[ix1]=0.5;
      wavesig[Mod(ix0-1,dimx)]=wavesig[Mod(ix1-1,dimx)]=-0.5;
    } 
  }
  
  anorma21D(dimx,wavesig,NULL);
} // end of wavelet1D_define_unit


/***************************************************************************/
void wavelet1D_define( int dimx, double sc, int wav, int ord_der,
		       double *wavesig ) {
  /***************************************************************************
   * Function computing the wavelet function, which can be a lorentzian
   * wavelet (wav>0) or a gaussian wavelet (wav<0)
   *
   * Parse :
   *    - n : dimension of the output wavelet function (this will be
   *      typically TIMESSC*sc),
   *    - sc : scale of analysis, 
   *    - expon : exponent of the lorentzian wavelet when >0, otherwise 
   *      it is a gaussian wavelet,
   *    - order : order of derivative of the gaussian.
   * Returns the wavelet function in wavesig.
   ***************************************************************************/

  int ix;
  double x;
  
  for( ix=0; ix<dimx; ix++ )  {

    x = (double)ix;
    if(ix >= dimx/2) x -= (double)dimx;
    x /= sc;

    if(wav == WAVMORLET) /* morlet */
      wavesig[ix] = morlet1D_coeff( x );
    else if(wav == WAVGAUSS) /* gaussian */  
      wavesig[ix] = gauss1D_coeff( x, ord_der ); 
    else if(wav == WAVHAAR) /* haar */ 
      wavesig[ix] = haar1D_coeff( x/(double)dimx );
    else if(wav >= SWAVLORENTZ)  /* lorentzian */      
      wavesig[ix] = lorentz1D_coeff( x, (double)wav/2., ord_der ); 
    /* Lorentzian of order 0.5 : wav=1 => exponent=0.5.
     * Lorentzian of order 1.  : wav=2 => exponent=1.
     * Lorentzian of order 1.5 : wav=3 => exponent=1.5.
     */
    /*  wave[ix]/=sc; : this is now done through wavelet1D_normaliza */
  }

} // end of wavelet1D_define


/***************************************************************************/
double haar1D_coeff( double x ) {
  /***************************************************************************/

  if((fabs(x)>0.5)) return 0.;
  else if(x>0)      return 1.;
  else              return -1.;
} // end of haar1D_coeff


/***************************************************************************/
double lorentz1D_coeff( double t, double expon, int ord_der ) {
  /***************************************************************************
   * Called by : wavelet1D_define
   ***************************************************************************/
  double TT=t*t,e2=2.*expon;

  switch (ord_der) {
  case 1:    return -e2 * t * pow( 1. + TT, -expon - 1. );
  case 2:    return e2 * (TT *(e2+1.) - 1.) * pow(1. + TT, -expon - 2.);
  default:
  case 0:    return pow(1. + TT, -expon);
  }
} // end of lorentz1D_coefficient


/***************************************************************************/
double morlet1D_coeff( double t ) {
  /***************************************************************************
   * Called by : wavelet1D_define
   ***************************************************************************/
  return exp(-.5*t*t) * cos(5.*t);
} // end of morlet1D_coeff


/***************************************************************************/
double gauss1D_coeff( double t, int ord_der ) {
  /***************************************************************************
   * Computes continuous Gaussian wavelet functions (0 to 5th derivative).
   * Called by : wavelet1D_define
   ***************************************************************************/
  double TT=t*t;
  
  switch (ord_der) {
  case 1:  return -t * exp(-.5*TT);
  case 2:  return (TT-1.) * exp(-.5*TT);
  case 3:  return t * exp(-.5*TT) * (3.-TT);
  case 4:  return exp(-.5*TT) * (pow(TT,2.) - 6.*TT + 3.);
  case 5:  return -t * exp(-.5*TT)* (pow(t,4.) - 10.*TT + 15.);
  default:
  case 0:  return exp(-.5*TT);
  }
} // end of gauss1D_coefficient



/***************************************************************************/
void wavelet1D_normaliza( int dimx, double sc, int ord_der, double *wavesig) {
  /***************************************************************************/
  double norma;
  int ix;
  
  //for(ix=0,norma=0.;ix<dimx;ix++) norma+=fabs(wave[ix]);
  //norma*=pow(sc, 1-ord_der);
  norma = pow( sc, (double)(DIM1D-ord_der) ); // DIM1D=1: dimension of the space
  for(ix=0;ix<dimx;ix++) wavesig[ix] /= norma;
} // end of wavelet1D_normaliza
