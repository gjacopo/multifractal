
/* ===================================
** mf_estimation.c
** started on Mon Jan 29 16:14:11 2007 
** ===================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Personnal libraries */
#include <utils.h>
#include <utl_alloc.h>
#include <utl_operator.h>
#include <utl_stats.h>

#include <filter.h>
#include <flt_stats1d.h>
#include <flt_stats2d.h>
#include <flt_fft1d.h>
#include <flt_fft2d.h>
#include <flt_deriva1d.h>
#include <flt_deriva2d.h>

#include <mfractal.h>
#include <mf_wtmm1d.h>

#include <mf_distribution.h>
#include <mf_estimation.h>

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

#endif /*!_PARSE_FILTER_PARAMETERS_*/

#ifdef _PARSE_FRACTAL_PARAMETERS_
#include <mf_parse.h>
extern ParFRAC *p_frac;
extern ParWTMM *p_wtmm;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */


extern int Verbose;
extern char *typ_name;

#ifndef _DATA_NAME_IN_FORMAT_
#define _DATA_NAME_IN_FORMAT_ "%s-N%05d"
#endif
#ifndef _DH_NAME_IN_FORMAT_
#define _DH_NAME_IN_FORMAT_ "Dh_%s_%s-N%d"
#endif
#ifndef _DH_NAME_OUT_FORMAT_
#define _DH_NAME_OUT_FORMAT_ _DH_NAME_IN_FORMAT_
#endif



/* Parameters for the 1D and 2D wavelet analysis */

extern double D0[3]; 

/* minimum scales (SC0) and scale steps (QS).
 * They depend on the derivative order (first argument) and on
 * the wavelet, defined by WAV (second argument) 
 */

const double EXP_WL_1D[4] = { -1.  , 0.5  , 1.,  1.5   }; 
/* wavelet exponents 1D (negative for gaussian) */

const double SC0_1D[3][4] =  
{ { 1.000, 0.500, 1.000, 0.500}, // Derivative order 0
  { 2.000, 2.000, 2.000, 2.000}, // Derivative order 1
  { 2.000, 1.000, 2.000, 2.000}  // Derivative order 2
}; 

const double QS_1D[3][4] =   
{ { 1.750, 1.125, 1.750, 1.500}, // Derivative order 0
  { 3.000, 3.000, 3.000, 3.000}, // Derivative order 1
  { 3.000, 2.000, 2.000, 3.000}  // Derivative order 2
};
 
const double EXP_WL_2D[4] = { -1.  , 0.5,  1.  , 1.5  }; 
/* wavelet exponents 2D  (negative for gaussian) */

const double SC0_2D[3][4] =  
{ { 1.000, 0.500, 0.500, 1.0000}, // Derivative order 0
  { 2.000, 2.000, 2.000, 2.0000}, // Derivative order 1
  { 2.000, 2.000, 2.000, 2.0000}  // Derivative order 2
}; 

const double QS_2D[3][4] =   
{ { 1.750, 1.125, 1.250, 1.500}, // Derivative order 0
  { 3.000, 3.000, 3.000, 3.000}, // Derivative order 1
  { 3.000, 2.000, 2.000, 3.000}  // Derivative order 2
}; 



/***************************************************************************/
double calcula_multifractal( int dimx, double **signal, double **expon, double *mm) {
/***************************************************************************/

#ifdef _PARSE_FRACTAL_PARAMETERS_
  if(p_frac->dim_space == DIM1D) 
#else
    if(DSPACE == DIM1D) 
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
      calcula1D_multifractal(dimx,signal[0],expon[0],mm);
    else
      calcula2D_multifractal(dimx,dimx,signal,expon,mm); 
}


/***************************************************************************/
double calcula1D_multifractal( int dimx, double *signal, double *expon, double *mm_h) {
/***************************************************************************/
  double **wave;
  double *dmasa;
  double *xx,*yy;
  double sc,sc0;
  double a,b,corr;
  double qs_1D;
  int xeff;
  int Nbuen=0;
  int ix,ip;

  int flag_ana1D, npoints;
  int wav, ord_der;
  double wav_range;
  double thetau;

#ifdef _PARSE_FRACTAL_PARAMETERS_
  flag_ana1D = p_frac->flag_ana1D;
  npoints = p_frac->npoints;
#else
  flag_ana1D = FLAG_ANA1D;
  npoints = NPOINTS;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
#ifdef _PARSE_FILTER_PARAMETERS_
  wav = p_wav->wav, wav_range = p_wav->wav_range; 
  thetau = p_wav->thetau;
  ord_der = p_wav->ord_der; 
#else
  wav = WAV, wav_range = WAVRANGE;
  ord_der = ORDDER;
  thetau = THETA0;
#endif/* !_PARSE_FILTER_PARAMETERS_ */

  xeff=redimensiona(dimx); // xeff=mdimensiona(dimx);
      
  TrackNullAlloc( wave=matrix2D(NPOINTS+1,xeff) );
  TrackNullAlloc( dmasa=(double*)calloc(xeff,sizeof(double)) );
  TrackNullAlloc( xx=(double*)calloc(NPOINTS+1,sizeof(double)) );
  TrackNullAlloc( yy=(double*)calloc(NPOINTS+1,sizeof(double)) );

  copy1D(dimx,signal,dmasa,NULL);
#ifdef _PARSE_FRACTAL_PARAMETERS_
  if(p_frac->flag_holder == FALSE)   
#else
    if(FLAGHOLDER == FALSE)   
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
      {
	gradient1D(xeff,dmasa);
	dmasa[dimx-1] = dmasa[xeff-1] = 0.; // to avoid boundary jumps
	for(ix=0;ix<xeff;ix++) dmasa[ix] = fabs(dmasa[ix]);
      }
  
  /*   Introducing the derivative order  */
  for(ip=0;ip<ord_der;ip++)    {
    IFVERBOSE WarningV("Applying the %d-th derivative\n",ip+1);
    gradient1D(xeff,dmasa);
  }
  
  /*    Calculating the wavelet projections   */
  sc0 = wavelet1D_escala(xeff,SC0_1D[ord_der][wav] ); 
  // EXP_WL_1D[wav] ord_der;
  qs_1D = pow(wavrange,1./((double)npoints));
  for( ip=0,sc=S0*SC0_1D[ord_der][wav];
       ip<=npoints;
       ip++,sc*=qs_1D/*QS_1D[ord_der][wav]*/)    {
    xx[ip]=1./log(sc0*sc);
    wavelet1D_genera(xeff,sc/*SC0_1D[ord_der][wav]*/);
    // EXP_WL_1D[wav],ord_der,wave[ip] 
    convuelto1D(xeff,dmasa,wave[ip]);
  }
  
  for( ix=0; ix<dimx; ix++ ) {
    for( ip=0; ip<=npoints; ip++ )      {
      a = fabs(wave[ip][ix]);
      if(a > 1.e-30) yy[ip] = log(a);
      else           yy[ip] = -30. * log(10.);
      yy[ip] *= xx[ip];
    }
    fit( xx, yy, npoints+1, &a, &b, &corr );
    if(fabs(corr) > NBUEN) Nbuen++;
    expon[ix] = b;
    mm_h[0] = fMin(mm_h[0],b);
    mm_h[1] = fMax(mm_h[1],b);
  }
  
  free(xx), free(yy);

  IFVERBOSE    {
    Warning("Multifractal analysis\n=====================\n");
    switch(wav)	{
    case WAVGAUSS:
      Warning("Gaussian wavelet\n"); break;
    default:
      WarningV("Lorentzian wavelet of exponent %f\n",(double)wav/2./*EXP_WL_1D[wav]*/);
      break;
    }
    if(ord_der) WarningV("in %d-th derivative\n",ord_der);
    WarningVV("Singularities: minimum: %f and maximum: %f\n",mm_h[0],mm_h[1]);
    WarningV("Percentage of good regression points: %f %%\n",
	     100.*((double)Nbuen)/((double)dimx));
  }
  
  free_matrix2D(wave,npoints+1);
  free(dmasa);
	
  return(((double)Nbuen)/((double)dimx));
} // end of calcula1D_multifractal



/***************************************************************************/
double calcula2D_multifractal( int dimx, int dimy, double **signal,
				double **expon, double *mm_h) {
/***************************************************************************/
  double ***wave;
  double **dmasa;
  double *xx,*yy;
  double sc,sc0;
  double a,b,corr;
  double qs_2D;
  
  int xeff,yeff;
  int Nbuen=0;
  int ix,iy,ip;

  int flag_1Danalysis, npoints;
  int wav, ord_der;
  double wav_range;
  double thetau;

#ifdef _PARSE_FRACTAL_PARAMETERS_
  flag_1Danalysis = p_frac->flag_1Danalysis;
  npoints = p_frac->npoints;
#else
  flag_1Danalysis = FLAG_1DANALYSIS;
  npoints = NPOINTS;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
#ifdef _PARSE_FILTER_PARAMETERS_
  wav = p_wav->wav, wav_range = p_wav->wav_range; 
  thetau = p_wav->thetau;
  ord_der = p_wav->ord_der; 
#else
  wav = WAV, wav_range = WAVRANGE;
  ord_der = ORDDER;
  thetau = THETA0;
#endif/* !_PARSE_FILTER_PARAMETERS_ */

  xeff=redimensiona(dimx),yeff=redimensiona(dimy);
  
  TrackNullAlloc( wave=matrix3D(NPOINTS+1,yeff,xeff) );
  TrackNullAlloc( dmasa=matrix2D(yeff,xeff) );
  TrackNullAlloc( xx=(double*)calloc(NPOINTS+1,sizeof(double)) );
  TrackNullAlloc( yy=(double*)calloc(NPOINTS+1,sizeof(double)) );

  copy(dimx,dimy,signal,dmasa,NULL);
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  if(p_frac->flag_holder == FALSE)   
#else
    if(FLAG_HOLDER == FALSE)   
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
      {
	IF(flag_1Danalysis) 
	  modderiva_line(xeff,yeff,dmasa,thetau);
	ELSE 
	  modderiva(xeff,yeff,dmasa);
      }
  if(ord_der) filtro2D(xeff,yeff,(double)ord_der,dmasa);
  
  /*    Calculating the wavelet projections   */  
  IF(flag_1Danalysis) 
    sc0 = wavelet2D_escala_line( dimx,dimy,SC0_2D[ord_der][wav] );
  // thetau,EXP_WL_2D[wav],ord_der;
  ELSE 
    sc0 = wavelet2D_escala( dimx,dimy,SC0_2D[ord_der][wav] );
  //  EXP_WL_2D[wav],ord_der
  
  qs_2D = pow(wavrange,1./((double)npoints));
  for(ip=0,sc=S0*SC0_2D[ordder][wav];
      ip<=npoints;
      ip++,sc*=qs_2D/*QS_2D[ord_der][wav]*/)    {
    xx[ip] = 1. / log(sc*sc0);
    IF(flag_1Danalysis) 
      wavelet2D_genera_line(xeff,yeff,sc, // thetau,EXP_WL_2D[wav],ord_der,
			    wave[ip]);
    ELSE 
      wavelet2D_genera(xeff,yeff,sc, // ord_der,EXP_WL_1D[wav],
		       wave[ip]);
    
    convuelto2D(xeff,yeff,dmasa,wave[ip]);
    anorma1_line(xeff,yeff,wave[ip],thetau,NULL);
  }

  for( iy=0; iy<dimy; iy++ )
    for( ix=0; ix<dimx; ix++ ) {
      for(ip=0;ip<=npoints;ip++)	{
	a = fabs(wave[ip][iy][ix]);
	if(a > 1.e-30) yy[ip] = log(a);
	else           yy[ip] = -30. * log(10.);
	yy[ip] *= xx[ip];
      }
      fit( xx, yy, npoints+1, &a, &b, &corr );
      if(fabs(corr) > NBUEN) Nbuen++;
      expon[iy][ix] = b;
      mm_h[0] = fMin(mm_h[0],b);
      mm_h[1] = fMax(mm_h[1],b);
    }

  free(xx), free(yy);
  
  IFVERBOSE    {
    Warning("Multifractal analysis\n=====================\n");
    switch(wav)      {
    case WAVGAUSS:
      Warning("Gaussian wavelet\n"); break;
    default:
      WarningV("Lorentzian wavelet of exponent %f\n",(double)wav/2.); break;      
    }
    
    if(ord_der) WarningV("in %d-th derivative\n",ord_der);
    WarningVV("Singularities: minimum: %f and maximum: %f\n",mm_h[0],mm_h[1]);
    WarningV("Percentage of good regression points: %f %%\n",
	   100.*((double)Nbuen)/((double)dimx*dimy));
  }
  
  free_matrix3D(wave,npoints+1,yeff);
  free_matrix2D(dmasa,yeff);
  
  return(((double)Nbuen)/((double)dimx*dimy));
} // end of calcula2D_multifractal



/***************************************************************************/
double Dh_registra( int Nr, double sc0, double *h, double *Dh, double *errDh) {
  /***************************************************************************/
  double *errDh_unif, h_unif;
  double Dh_th;
  double hmin,hmax,hrad;
  double mean_err,av_err,quad_err,sigma,norma;
  double weight, y, sc;
  const int Nh=NPOINTS_SPECTRUM; // Number of points to be solved in the spectrum
  int Nerr;
  int ix,ip,ip1,ih,ih1,ih2;
  
  /* Generating close to uniformly sampled D(h) from data  */
  TrackNullAlloc( errDh_unif=(double*)calloc(Nh,sizeof(double)) );
  theoretical_Deltah( &hmin, &hmax );
  hrad = (hmax-hmin) / ((double)Nh); // Step and uncertainty radius

  ip = 0;
  for( ih=0,h_unif=hmin; ih<Nh; ih++,h_unif+=hrad )    {
    Nerr = 0;
    errDh_unif[ih] = 0.;
    if(ip < Nr)      {
      for(; (h[ip]<h_unif-hrad/2)&&(ip<Nr); ip++ ) ;
      for(; (h[ip]<h_unif+hrad/2)&&(ip<Nr); ip++ )	{
	Dh_th = theoretical_Dh( h[ip] );
	if((Dh[ip]<=
#ifdef _PARSE_FRACTAL_PARAMETERS_
	    p_frac->dim_space
#else
	    DSPACE
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
	    ) && (Dh_th>0.)) {  // Reasonable points
	  Nerr++;
	  errDh_unif[ih] += Dh_th-Dh[ip];
	}
      }
    }
    if(Nerr > 0) errDh_unif[ih] /= (double)Nerr;
    else errDh_unif[ih] = theoretical_Dh( h[ip] ); 
  }
  
  /* Calculating errors */
  norma = mean_err = av_err = quad_err = 0.;
  for( ih=0,h_unif=hmin+hrad/2; ih<Nh; ih++,h_unif+=hrad ) {
    weight = 1.;    // exp(-Dh_th*log(p_frac->sc0));
    mean_err += errDh_unif[ih] * weight;
    av_err += fabs(errDh_unif[ih]) * weight;
    quad_err += errDh_unif[ih] * errDh_unif[ih] * weight;
    norma += weight;
  }
  
  mean_err /= norma;
  av_err /= norma;
  sigma = sqrt(quad_err/norma - mean_err*mean_err);
  quad_err = sqrt(quad_err/norma);
  
  IFVERBOSE {
    WarningVV("Mean err: %0.2f; std. deviation: %0.2f",mean_err,sigma);
    WarningVV("Typical deviation: %0.2f; quad. error: %0.2f",av_err,quad_err);
    Warning("");
  }
  
  /* Freeing memory before finishing */ 
  free(errDh_unif);
 
  return av_err/(hmax-hmin);
} // end of registra_Dh



/***************************************************************************/
void expon_genera( int leff, double **signal, double **expon) {
  /***************************************************************************/
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  if(p_frac->dim_space == DIM1D) 
#else
    if(DSPACE == DIM1D) 
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
      expon1D_genera(leff,signal[0],expon[0]);
    else 
      expon2D_genera(leff,signal,expon);
} // end of expon_genera


/***************************************************************************/
int expon1D_genera( int leff, double *serie, double *expon) {
/***************************************************************************/
  double meang;
  int ix;

  copy1D(leff,serie,expon,NULL);
  gradient1D(leff,expon);

  for(ix=0;ix<leff;ix++) expon[ix]=fabs(expon[ix]);
#ifdef _PARSE_FRACTAL_PARAMETERS_
  IF(p_frac->inv_trans) 
#else
    IF(FLAG_INVTRANS) 
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
    meang = media1D(leff,expon,NULL);
  else meang = 1.;
  
  for(ix=0;ix<leff;ix++)    {
    if(expon[ix]/meang > 1.e-30) 
      expon[ix] = -log(expon[ix]/meang)/log((double)leff);
    else expon[ix] = 30.*log(10.)/log((double)leff);
  }
  
  return OK;
} // end of expon1D_genera


/***************************************************************************/
int expon2D_genera( int leff, double **image, double **expon) {
  /***************************************************************************/
  double **gy;
  double meang;
  int ix,iy;

  TrackNullAlloc( gy=matrix2D(leff,leff) );

  copy(leff,leff,image,expon, NULL);
  gradient2D(leff,leff,expon,gy);

  for(iy=0;iy<leff;iy++)
    for(ix=0;ix<leff;ix++)
      expon[iy][ix] = sqrt(expon[iy][ix]*expon[iy][ix]
			 + gy[iy][ix]*gy[iy][ix]);
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  IF(p_frac->inv_trans) 
#else
    IF(FLAG_INVTRANS) 
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
    meang = media(leff,leff,expon, NULL);
  else meang = ((double)leff);
  
  for(iy=0;iy<leff;iy++)
    for(ix=0;ix<leff;ix++)
      if(expon[iy][ix]/meang>1e-30) 
	expon[iy][ix]=-log(expon[iy][ix]/meang)/log((double)leff);
      else expon[iy][ix]=30.*log(10.)/log((double)leff);
  
  free_matrix2D(gy,leff);

  return OK;
} // end of expon2D_genera


