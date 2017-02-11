/********************************************************/
/*                                                      */
/*	FFT.c - Version del 21 de Septiembre, 2005      */
/*                                                      */
/********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Personnal libraries */
#include <utils.h>
#include <utl_alloc.h>
#include <utl_stats.h>
#include <utl_operator.h>

#include <filter.h>
#include <flt_fft2d.h>

#ifndef M_PI
#define M_PI   3.1415926535997932
#endif

#ifdef _PARSE_FILTER_PARAMETERS_
#include <flt_parse.h>
extern ParFILT* p_fil;
#ifndef redimensiona
#define redimensiona(x) (p_fil->flag_memory?(x):dimensiona(x))
#endif

#else
#ifndef redimensiona
#define redimensiona(x) (FLAG_MEMORY?(x):dimensiona(x))
#endif

#endif/* !_PARSE_FILTER_PARAMETERS_ */


/*      Function declarations    */

/***************************************************************************/
void Fourier2D(int dimx, int dimy, double **funcionR, double **funcionI, 
	       int signo) {

  /***************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  IF(p_fil->flag_memory)
#else
    IF(FLAG_MEMORY)
#endif/* !_PARSE_FILTER_PARAMETERS_ */
    FFFT2D(dimx,dimy,funcionR,funcionI,signo);
  else 
    FFT2D(dimx,dimy,funcionR,funcionI,signo);
}


/***************************************************************************/
void convuelve2D( int dimx, int dimy, double **f1, double **f2, 
		   double **salida) {
/***************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  IF(p_fil->flag_memory)
#else
    IF(FLAG_MEMORY)
#endif/* !_PARSE_FILTER_PARAMETERS_ */
    convuelve_FFFT(dimx,dimy,f1,f2,salida);
  else  
    convuelve_FFT(dimx,dimy,f1,f2,salida);
}


/***************************************************************************/
void convuelto2D( int dimx, int dimy, double **f1, double **f2) {
/***************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  IF(p_fil->flag_memory)
#else
    IF(FLAG_MEMORY)
#endif/* !_PARSE_FILTER_PARAMETERS_ */
    convuelto_FFFT(dimx,dimy,f1,f2);
  else  
    convuelto_FFT(dimx,dimy,f1,f2);
}


/***************************************************************************/
void convuelto2D_vec( int dimx, int dimy, double **f1x, double **f1y,
		       double **f2x, double **f2y) {
/***************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  IF(p_fil->flag_memory)
#else
    IF(FLAG_MEMORY)
#endif/* !_PARSE_FILTER_PARAMETERS_ */
    convuelto_vec_FFFT(dimx,dimy,f1x,f1y,f2x,f2y);
  else  
    convuelto_vec_FFT(dimx,dimy,f1x,f1y,f2x,f2y);
}


/***************************************************************************/
void deconvuelve2D( int dimx, int dimy, double **f1, double **f2, 
		     double **salida) {
/***************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  IF(p_fil->flag_memory)
#else
    IF(FLAG_MEMORY)
#endif/* !_PARSE_FILTER_PARAMETERS_ */
    deconvuelve_FFFT(dimx,dimy,f1,f2,salida);
  else  
    deconvuelve_FFT(dimx,dimy,f1,f2,salida);
}


/***************************************************************************/
void deconvuelto2D( int dimx, int dimy, double **f1, double **f2) {
 /***************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  IF(p_fil->flag_memory)
#else
    IF(FLAG_MEMORY)
#endif/* !_PARSE_FILTER_PARAMETERS_ */
    deconvuelto_FFFT(dimx,dimy,f1,f2);
  else  
    deconvuelto_FFT(dimx,dimy,f1,f2);
}


/***************************************************************************/
void deconvuelto2D_vec( int dimx, int dimy, double **f1x, double **f1y,
			 double **f2x, double **f2y) {
/***************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  IF(p_fil->flag_memory)
#else
    IF(FLAG_MEMORY)
#endif/* !_PARSE_FILTER_PARAMETERS_ */
    deconvuelto_vec_FFFT(dimx,dimy,f1x,f1y,f2x,f2y);
  else  
    deconvuelto_vec_FFT(dimx,dimy,f1x,f1y,f2x,f2y);
}


/***************************************************************************/
void FFT2D(int dimx, int dimy, double **funcionR, double **funcionI, 
	   int signo) {
/***************************************************************************/
  FFThorizontal(dimx,dimy,funcionR,funcionI,signo);
  FFTvertical(dimx,dimy,funcionR,funcionI,signo);
}


/***************************************************************************/
void FFThorizontal(int dimx, int dimy, double **funcionhR, double **funcionhI,
		   int signo) {
/***************************************************************************/

  double tempR,tempI,wpasoR,wpasoI,wwR,wwI;
  int iy;
  int ix,je,mm,mmax,istep;

  for(iy=0;iy<dimy;iy++)    {
    
    /*	INICIO DE LINEA HORIZONTAL		*/
    je=1;
    for(ix=0;ix<dimx;ix++)	{
      if(je>ix+1)	    {
	tempR=funcionhR[iy][je-1];
	tempI=funcionhI[iy][je-1];
	funcionhR[iy][je-1]=funcionhR[iy][ix];
	funcionhI[iy][je-1]=funcionhI[iy][ix];
	funcionhR[iy][ix]=tempR;
	funcionhI[iy][ix]=tempI;
      }
      mm=dimx/2;      
      while((mm>1)&&(je>mm))	{
	je=je-mm;
	mm=mm/2;
      }
      je=je+mm;
    }	
    
    mmax=1;
    while(dimx>mmax)	{
      istep=2*mmax;
      wpasoR=cos(M_PI/((double) mmax));
      wpasoI=signo*sin(M_PI/((double) mmax));
      wwR=1.;
      wwI=0.;
      
      for(mm=1;mm<=mmax;mm++)	{
	for(ix=mm-1;ix<dimx;ix+=istep)		{
	  je=ix+mmax;
	  C_mult(wwR,wwI,funcionhR[iy][je],funcionhI[iy][je],
		 &tempR,&tempI);
	  funcionhR[iy][je]=funcionhR[iy][ix]-tempR;
	  funcionhI[iy][je]=funcionhI[iy][ix]-tempI;
	  funcionhR[iy][ix]=funcionhR[iy][ix]+tempR;
	  funcionhI[iy][ix]=funcionhI[iy][ix]+tempI;
	}
	C_mult(wwR,wwI,wpasoR,wpasoI,&wwR,&wwI);
      }
      mmax=istep;
      
    }
    
    for(ix=0;ix<dimx;ix++)     {
      funcionhR[iy][ix]=funcionhR[iy][ix]/sqrt((double)dimx);
      funcionhI[iy][ix]=funcionhI[iy][ix]/sqrt((double)dimx);
    }
    
    /*	FIN DE LINEA HORIZONTAL			*/  
  }
}


/***************************************************************************/
void FFTvertical(int dimx, int dimy, double **funcionvR, double **funcionvI, 
		 int signo) {
/***************************************************************************/

  double tempR,tempI,wpasoR,wpasoI,wwR,wwI;
  int ix;
  int iy,je,mm,mmax,istep;

  for(ix=0;ix<dimx;ix++)    {
    
    /*	INICIO DE LINEA VERTICAL		*/
    je=1;
    for(iy=0;iy<dimy;iy++)      {
      if(je>iy+1)	    {
	tempR=funcionvR[je-1][ix];
	tempI=funcionvI[je-1][ix];
	funcionvR[je-1][ix]=funcionvR[iy][ix];
	funcionvI[je-1][ix]=funcionvI[iy][ix];
	funcionvR[iy][ix]=tempR;
	funcionvI[iy][ix]=tempI;
      }
      mm=dimy/2;  
      while((mm>1)&&(je>mm))	{
	je=je-mm;
	mm=mm/2;
      }
      je=je+mm;      
    }	
    
    mmax=1;
    while(dimy>mmax)      {
      istep=2*mmax;
      wpasoR=cos(M_PI/((double) mmax));
      wpasoI=signo*sin(M_PI/((double) mmax));
      wwR=1.;
      wwI=0.;
      
      for(mm=1;mm<=mmax;mm++)	    {
	for(iy=mm-1;iy<dimy;iy+=istep)	  {
	  je=iy+mmax;
	  C_mult(wwR,wwI,funcionvR[je][ix],funcionvI[je][ix],
		 &tempR,&tempI);
	  funcionvR[je][ix]=funcionvR[iy][ix]-tempR;
	  funcionvI[je][ix]=funcionvI[iy][ix]-tempI;
	  funcionvR[iy][ix]=funcionvR[iy][ix]+tempR;
	  funcionvI[iy][ix]=funcionvI[iy][ix]+tempI;
	}
	C_mult(wwR,wwI,wpasoR,wpasoI,&wwR,&wwI);
      }
      mmax=istep;
      
    }
    
    for(iy=0;iy<dimy;iy++)	{
      funcionvR[iy][ix]=funcionvR[iy][ix]/sqrt((double)dimy);
      funcionvI[iy][ix]=funcionvI[iy][ix]/sqrt((double)dimy);
    }    
    
    /*	FIN DE LINEA VERTICAL			*/
  }
}


/***************************************************************************/
void FFFT2D( int dimx, int dimy, double **funcionR, double **funcionI, 
	     int signo) {
  /***************************************************************************/
  FFFThorizontal(dimx,dimy,funcionR,funcionI,signo);
  FFFTvertical(dimx,dimy,funcionR,funcionI,signo);
}


/***************************************************************************/
int FFFThorizontal( int dimx, int dimy,  double **funcionR, double **funcionI, 
		     int signo) {
/***************************************************************************/
  double **workR,**workI;

  int inu;
  int dima,dimb,dimc;
  int ix,iy,i1,i2,i3;

  TrackNullAlloc( workR=matrix2D(dimy,dimx) );
  TrackNullAlloc( workI=matrix2D(dimy,dimx) );

  for(iy=0;iy<dimy;iy++)    {
    
    dima=1;
    dimb=dimx;
    dimc=1;    
    inu=1;
    
    while(dimb>1)	{
      dima=dimc*dima;
      dimc=2;
      while(dimb%dimc!=0) dimc++;
      
      dimb=dimb/dimc;
      
      if(inu==1) PFFThorizontal(dima,dimb,dimc,iy,funcionR,funcionI,
				workR,workI,signo);
      else PFFThorizontal(dima,dimb,dimc,iy,workR,workI,funcionR,funcionI,
			  signo);
      inu=1-inu;
    }
    
    if(inu==0)	
      for(ix=0;ix<dimx;ix++)	{
	funcionR[iy][ix]=workR[iy][ix];
	funcionI[iy][ix]=workI[iy][ix];
      }
    
    for(ix=0;ix<dimx;ix++)      {
      funcionR[iy][ix]=funcionR[iy][ix]/sqrt((double)dimx);
      funcionI[iy][ix]=funcionI[iy][ix]/sqrt((double)dimx);
    } 
  }
  
  free_matrix2D(workR,dimy);
  free_matrix2D(workI,dimy);

  return OK;
}


/***************************************************************************/
void PFFThorizontal( int dima, int dimb, int dimc, int iy, 
		     double **uinR, double **uinI, 
		     double **uoutR, double **uoutI, int sign) {
  /***************************************************************************
   * Auxiliary routine called by FFThorizontal
   ***************************************************************************/
  double angle;
  double deltaR,deltaI;
  double omegaR=1.,omegaI=0.;
  double sumR,sumI;
  int ia,ib,ic,ix,jcr,jc;

  angle=2.*M_PI/((double)(dima*dimc));

  deltaR = cos(angle);
  deltaI = ((double)sign)*sin(angle);

  for(ic=0;ic<dimc;ic++)    
    for(ia=0;ia<dima;ia++)      {
      for(ib=0;ib<dimb;ib++)	{
	sumR = sumI = 0.;
	for(jcr=0;jcr<dimc;jcr++)	  {
	  jc = dimc - 1 - jcr;
	  ix = ib + dimb * (jc+ia*dimc);
	  C_mult(omegaR,omegaI,sumR,sumI,&sumR,&sumI);
	  sumR += uinR[iy][ix];
	  sumI += uinI[iy][ix];
	}
	ix = ib + dimb * (ia+ic*dima);
	uoutR[iy][ix] = sumR;
	uoutI[iy][ix] = sumI;
      }
      C_mult(omegaR,omegaI,deltaR,deltaI,&omegaR,&omegaI);
    }
  
}


/***************************************************************************/
int FFFTvertical( int dimx, int dimy,  double **funcionR, double **funcionI, 
		   int signo) {
/***************************************************************************/
  double **workR,**workI;

  int inu;
  int dima,dimb,dimc;
  int ix,iy,i1,i2,i3;

  TrackNullAlloc( workR=matrix2D(dimy,dimx) );
  TrackNullAlloc( workI=matrix2D(dimy,dimx) );
  
  for(ix=0;ix<dimx;ix++)    {
    dima=1;
    dimb=dimy;
    dimc=1;
    inu=1;
    
    while(dimb>1)	{
      dima=dimc*dima;
      dimc=2;
      while(dimb%dimc!=0) dimc++;
      
      dimb=dimb/dimc;
      
      if(inu==1) PFFTvertical(dima,dimb,dimc,ix,funcionR,funcionI,
			      workR,workI,signo);
      else PFFTvertical(dima,dimb,dimc,ix,workR,workI,funcionR,funcionI,
			signo);
      inu=1-inu;
    }
    
    if(inu==0)      
      for(iy=0;iy<dimy;iy++)	  {
	funcionR[iy][ix]=workR[iy][ix];
	funcionI[iy][ix]=workI[iy][ix];
      }
    
    for(iy=0;iy<dimy;iy++)	{
      funcionR[iy][ix]=funcionR[iy][ix]/sqrt((double)dimy);
      funcionI[iy][ix]=funcionI[iy][ix]/sqrt((double)dimy);
    }
    
  }
  
  free_matrix2D(workR,dimy);
  free_matrix2D(workI,dimy);

  return OK;
}


/***************************************************************************/
void PFFTvertical( int dima, int dimb, int dimc, int ix,
		   double **uinR, double **uinI, 
		   double **uoutR, double **uoutI, int sign) {
  /***************************************************************************
   * Auxiliary routine called by FFTvertical
   ***************************************************************************/
  double angle;
  double deltaR,deltaI;
  double omegaR=1.,omegaI=0.;
  double sumR,sumI;
  int ia,ib,ic,iy,jcr,jc;

  angle = 2.*M_PI / ((double)(dima*dimc));

  deltaR = cos(angle);
  deltaI = ((double)sign)*sin(angle);

  for(ic=0;ic<dimc;ic++)   
    for(ia=0;ia<dima;ia++)   {
      for(ib=0;ib<dimb;ib++)	{
	sumR = sumI = 0.;
	for(jcr=0;jcr<dimc;jcr++)	  {
	  jc = dimc - 1 - jcr;
	  iy = ib + dimb * (jc+ia*dimc);
	  C_mult(omegaR,omegaI,sumR,sumI,&sumR,&sumI);
	  sumR += uinR[iy][ix];
	  sumI += uinI[iy][ix];
	}
	iy = ib + dimb * (ia+ic*dima);
	uoutR[iy][ix] = sumR;
	uoutI[iy][ix] = sumI;
      }
      C_mult(omegaR,omegaI,deltaR,deltaI,&omegaR,&omegaI);
    }
  
}


/***************************************************************************/
int convuelve_FFT( int dimx, int dimy, double **f1, double **f2, 
		    double **salida) {
  /***************************************************************************/
  double **If1,**If2,**Isalida;
  double norma;
  int ix,iy;
  
  norma=sqrt((double)(dimx*dimy));
  TrackNullAlloc( If1=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2=matrix2D(dimy,dimx) );
  TrackNullAlloc( Isalida=matrix2D(dimy,dimx) );

  fill0(dimx,dimy,If1,NULL);
  fill0(dimx,dimy,If2,NULL);
  fill0(dimx,dimy,Isalida,NULL);
  
  FFT2D(dimx,dimy,f1,If1,-1);
  FFT2D(dimx,dimy,f2,If2,-1);

  for(iy=0;iy<dimy;iy++)   
    for(ix=0;ix<dimx;ix++)      {
      salida[iy][ix]=norma*(f1[iy][ix]*f2[iy][ix]
			    -If1[iy][ix]*If2[iy][ix]);
      Isalida[iy][ix]=norma*(f1[iy][ix]*If2[iy][ix]
			     +If1[iy][ix]*f2[iy][ix]);
    }
  
  FFT2D(dimx,dimy,salida,Isalida,1);
  FFT2D(dimx,dimy,f1,If1,1);
  FFT2D(dimx,dimy,f2,If2,1);
  
  free_matrix2D(If1,dimy);
  free_matrix2D(If2,dimy);
  free_matrix2D(Isalida,dimy);
 
 return OK;
}


/***************************************************************************/
int convuelto_FFT( int dimx, int dimy, double **f1, double **f2) {
/***************************************************************************/
  int ix,iy;
  double **If1,**If2;
  double buffR,buffI;
  double norma;

  norma=sqrt((double)(dimx*dimy));
  TrackNullAlloc( If1=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2=matrix2D(dimy,dimx) );

  fill0(dimx,dimy,If1,NULL);
  fill0(dimx,dimy,If2,NULL);

  FFT2D(dimx,dimy,f1,If1,-1);
  FFT2D(dimx,dimy,f2,If2,-1);

  for(iy=0;iy<dimy;iy++)   
    for(ix=0;ix<dimx;ix++)      {
      buffR=f1[iy][ix]*f2[iy][ix]
	-If1[iy][ix]*If2[iy][ix];
      buffI=f1[iy][ix]*If2[iy][ix]
	+If1[iy][ix]*f2[iy][ix];
      f2[iy][ix]=norma*buffR;
      If2[iy][ix]=norma*buffI;
    }

  FFT2D(dimx,dimy,f1,If1,1);
  FFT2D(dimx,dimy,f2,If2,1);

  free_matrix2D(If1,dimy);
  free_matrix2D(If2,dimy);

  return OK;
}


/***************************************************************************/
int convuelto_vec_FFT( int dimx, int dimy, double **f1x, double **f1y,
			double **f2x, double **f2y) {
/***************************************************************************/
  int ix,iy;
  double **If1x,**If1y,**If2x,**If2y;
  double buffRx,buffIx;
  double buffRy,buffIy;
  double norma;

  norma=sqrt((double)(dimx*dimy));

  TrackNullAlloc( If1x=matrix2D(dimy,dimx) );
  TrackNullAlloc( If1y=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2x=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2y=matrix2D(dimy,dimx) );

  fill0(dimx,dimy,If1x,NULL);
  fill0(dimx,dimy,If1y,NULL);
  fill0(dimx,dimy,If2x,NULL);
  fill0(dimx,dimy,If2y,NULL);

  FFT2D(dimx,dimy,f1x,If1x,-1);
  FFT2D(dimx,dimy,f1y,If1y,-1);
  FFT2D(dimx,dimy,f2x,If2y,-1);
  FFT2D(dimx,dimy,f2y,If2y,-1);

  for(iy=0;iy<dimy;iy++)    
    for(ix=0;ix<dimx;ix++)      {
      buffRx=f1x[iy][ix]*f2x[iy][ix]-If1x[iy][ix]*If2x[iy][ix]
	-(f1y[iy][ix]*f2y[iy][ix]-If1y[iy][ix]*If2y[iy][ix]);
      buffIx=f1x[iy][ix]*If2x[iy][ix]+If1x[iy][ix]*f2x[iy][ix]
	-(f1y[iy][ix]*If2y[iy][ix]+If1y[iy][ix]*f2y[iy][ix]);
      
      buffRy=f1x[iy][ix]*f2y[iy][ix]-If1x[iy][ix]*If2y[iy][ix]
	+f1y[iy][ix]*f2x[iy][ix]-If1y[iy][ix]*If2x[iy][ix];
      buffIy=f1x[iy][ix]*If2y[iy][ix]+If1x[iy][ix]*f2y[iy][ix]
	+f1y[iy][ix]*If2x[iy][ix]+If1y[iy][ix]*f2x[iy][ix];
            
      f2x[iy][ix]=norma*buffRx;
      f2y[iy][ix]=norma*buffRy;
      If2x[iy][ix]=norma*buffIx;
      If2y[iy][ix]=norma*buffIy;
    }

  FFT2D(dimx,dimy,f1x,If1y,1);
  FFT2D(dimx,dimy,f1y,If1y,1);
  FFT2D(dimx,dimy,f2x,If2x,1);
  FFT2D(dimx,dimy,f2y,If2y,1);

  free_matrix2D(If1x,dimy);
  free_matrix2D(If1y,dimy);
  free_matrix2D(If2x,dimy);
  free_matrix2D(If2y,dimy);
 
 return OK;
}


/***************************************************************************/
int deconvuelve_FFT( int dimx, int dimy, double **f1, double **f2, 
		      double **salida) {
/***************************************************************************/
  int ix,iy;
  double **If1,**If2,**Isalida;
  double mod;

  TrackNullAlloc( If1=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2=matrix2D(dimy,dimx) );
  TrackNullAlloc( Isalida=matrix2D(dimy,dimx) );
  fill0(dimx,dimy,If1,NULL);
  fill0(dimx,dimy,If2,NULL);

  FFT2D(dimx,dimy,f1,If1,-1);
  FFT2D(dimx,dimy,f2,If2,-1);

  for(iy=0;iy<dimy;iy++)    
    for(ix=0;ix<dimx;ix++)      {
      mod=f2[iy][ix]*f2[iy][ix]+If2[iy][ix]*If2[iy][ix];
      if(mod>1e-30)	    {
	salida[iy][ix]=(f1[iy][ix]*f2[iy][ix]
			+If1[iy][ix]*If2[iy][ix])/mod;
	Isalida[iy][ix]=(-f1[iy][ix]*If2[iy][ix]
			 +If1[iy][ix]*f2[iy][ix])/mod;
      }  else {
	salida[iy][ix]=0.;
	Isalida[iy][ix]=0.;
      }
    }

  FFT2D(dimx,dimy,salida,Isalida,1);
  FFT2D(dimx,dimy,f1,If1,1);
  FFT2D(dimx,dimy,f2,If2,1);

  free_matrix2D(If1,dimy);
  free_matrix2D(If2,dimy);
  free_matrix2D(Isalida,dimy);

  return OK;
}


/***************************************************************************/
int deconvuelto_FFT( int dimx, int dimy, double **f1, double **f2) {
  /***************************************************************************/
  int ix,iy;
  double **If1,**If2;
  double buffR,buffI,mod;

  TrackNullAlloc( If1=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2=matrix2D(dimy,dimx) );
  fill0(dimx,dimy,If1,NULL);
  fill0(dimx,dimy,If2,NULL);
  
  FFT2D(dimx,dimy,f1,If1,-1);
  FFT2D(dimx,dimy,f2,If2,-1);
  
  for(iy=0;iy<dimy;iy++)  
    for(ix=0;ix<dimx;ix++)     {
      mod=f1[iy][ix]*f1[iy][ix]+If1[iy][ix]*If1[iy][ix];
      if(mod>1e-30)	{
	buffR=(f1[iy][ix]*f2[iy][ix]
	       +If1[iy][ix]*If2[iy][ix])/mod;
	buffI=(f1[iy][ix]*If2[iy][ix]
	       -If1[iy][ix]*f2[iy][ix])/mod;
	f2[iy][ix]=buffR;
	If2[iy][ix]=buffI;
      }	  else	    {
	f2[iy][ix]=0.;
	If2[iy][ix]=0.;
      }
    }

  FFT2D(dimx,dimy,f1,If1,1);
  FFT2D(dimx,dimy,f2,If2,1);

  free_matrix2D(If1,dimy);
  free_matrix2D(If2,dimy);

  return OK;
}


/***************************************************************************/
int deconvuelto_vec_FFT( int dimx, int dimy, double **f1x, double **f1y,
			  double **f2x, double **f2y) {
 /***************************************************************************/
  int ix,iy;
  double **If1x,**If1y,**If2x,**If2y;
  double buffRx,buffIx;
  double buffRy,buffIy;
  double modR,modI,modd;
  double norma;

  norma=sqrt((double)(dimx*dimy));

  TrackNullAlloc( If1x=matrix2D(dimy,dimx) );
  TrackNullAlloc( If1y=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2x=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2y=matrix2D(dimy,dimx) );
  
  fill0(dimx,dimy,If1x,NULL);
  fill0(dimx,dimy,If1y,NULL);
  fill0(dimx,dimy,If2x,NULL);
  fill0(dimx,dimy,If2y,NULL);

  FFT2D(dimx,dimy,f1x,If1x,-1);
  FFT2D(dimx,dimy,f1y,If1y,-1);
  FFT2D(dimx,dimy,f2x,If2y,-1);
  FFT2D(dimx,dimy,f2y,If2y,-1);

  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)      {
      
      modR=(f1x[iy][ix]*f1x[iy][ix]
	    -If1x[iy][ix]*If1x[iy][ix]
	    +f1y[iy][ix]*f1y[iy][ix]
	    -If1y[iy][ix]*If1y[iy][ix])*norma*norma;
      modI=(2.*f1x[iy][ix]*If1x[iy][ix]
	    +2.*f1y[iy][ix]*If1y[iy][ix])*norma*norma;
      
      modd=modR*modR+modI*modI;
      
      buffRx=(f1x[iy][ix]*f2x[iy][ix]
	      -If1x[iy][ix]*If2x[iy][ix]
	      +f1y[iy][ix]*f2y[iy][ix]
	      -If1y[iy][ix]*If2y[iy][ix])*norma;
      buffIx=(f1x[iy][ix]*If2x[iy][ix]
	      +If1x[iy][ix]*f2x[iy][ix]
	      +f1y[iy][ix]*If2y[iy][ix]
	      +If1y[iy][ix]*f2y[iy][ix])*norma;
      
      buffRy=(f1x[iy][ix]*f2y[iy][ix]
	      -If1x[iy][ix]*If2y[iy][ix]
	      -f1y[iy][ix]*f2x[iy][ix]
	      +If1y[iy][ix]*If2x[iy][ix])*norma;
      buffIy=(f1x[iy][ix]*If2y[iy][ix]
	      +If1x[iy][ix]*f2y[iy][ix]
	      -f1y[iy][ix]*If2x[iy][ix]
	      -If1y[iy][ix]*f2x[iy][ix])*norma;
      
      if(modd>1e-30)	    {
	f2x[iy][ix]=(buffRx*modR+buffIx*modI)/modd;
	If2x[iy][ix]=(buffIx*modR-buffRx*modI)/modd;
	
	f2y[iy][ix]=(buffRy*modR+buffIy*modI)/modd;
	If2y[iy][ix]=(buffIy*modR-buffRy*modI)/modd;	
      }	  else	    {
	f2x[iy][ix]=0.;
	If2x[iy][ix]=0.;
	f2y[iy][ix]=0.;
	If2y[iy][ix]=0.;
      }
      
    }
  
  FFT2D(dimx,dimy,f1x,If1y,1);
  FFT2D(dimx,dimy,f1y,If1y,1);
  FFT2D(dimx,dimy,f2x,If2x,1);
  FFT2D(dimx,dimy,f2y,If2y,1);

  free_matrix2D(If1x,dimy);
  free_matrix2D(If1y,dimy);
  free_matrix2D(If2x,dimy);
  free_matrix2D(If2y,dimy);

  return OK;
}


/***************************************************************************/
int band_pass( int dimx, int dimy, double fmin, double fmax, 
		double **func) {
  /***************************************************************************/
  double **fR,**fI;
  double x,y,f;
  int xeff,yeff;
  int ix,iy;

  xeff=dimensiona(dimx);
  yeff=dimensiona(dimy);

  TrackNullAlloc( fR=matrix2D(yeff,xeff) );
  TrackNullAlloc( fI=matrix2D(yeff,xeff) );
  
  fill0(dimx,dimy,fR,NULL);
  fill0(dimx,dimy,fI,NULL);
  
  copy(dimx,dimy,func,fR,NULL);
  FFT2D(xeff,yeff,fR,fI,-1);
  
  for(iy=0;iy<yeff;iy++)    {
    y=((double)iy)/((double)yeff);
    if(iy>=yeff/2) y-=1.;
    for(ix=0;ix<xeff;ix++)      {
      x=((double)ix)/((double)xeff);
      if(ix>=xeff/2) x-=1.;
      f=sqrt(x*x+y*y);
      if((f<fmin)||(f>fmax))	{
	fR[iy][ix]=0.;
	fI[iy][ix]=0.;
      }
    }
  }
  
  FFT2D(xeff,yeff,fR,fI,1);
  copy(dimx,dimy,fR,func,NULL);

  free_matrix2D(fR,yeff);
  free_matrix2D(fI,yeff);

  return OK;
}


/***************************************************************************/
int convuelve_FFFT( int dimx, int dimy, double **f1, double **f2, 
		     double **salida) {
/***************************************************************************/
  double **If1,**If2,**Isalida;
  double norma;
  int ix,iy;

  norma=sqrt((double)(dimx*dimy));
  TrackNullAlloc( If1=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2=matrix2D(dimy,dimx) );
  TrackNullAlloc( Isalida=matrix2D(dimy,dimx) );

  fill0(dimx,dimy,If1,NULL);
  fill0(dimx,dimy,If2,NULL);
  fill0(dimx,dimy,Isalida,NULL);

  FFFT2D(dimx,dimy,f1,If1,-1);
  FFFT2D(dimx,dimy,f2,If2,-1);

  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)      {
      salida[iy][ix]=norma*(f1[iy][ix]*f2[iy][ix]
			    -If1[iy][ix]*If2[iy][ix]);
      Isalida[iy][ix]=norma*(f1[iy][ix]*If2[iy][ix]
			     +If1[iy][ix]*f2[iy][ix]);
    }

  FFFT2D(dimx,dimy,salida,Isalida,1);
  FFFT2D(dimx,dimy,f1,If1,1);
  FFFT2D(dimx,dimy,f2,If2,1);

  free_matrix2D(If1,dimy);
  free_matrix2D(If2,dimy);
  free_matrix2D(Isalida,dimy);

  return OK;
}


/***************************************************************************/
int convuelto_FFFT( int dimx, int dimy, double **f1, double **f2) {
/***************************************************************************/
  int ix,iy;
  double **If1,**If2;
  double buffR,buffI;
  double norma;

  norma=sqrt((double)(dimx*dimy));
  TrackNullAlloc( If1=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2=matrix2D(dimy,dimx) );

  fill0(dimx,dimy,If1,NULL);
  fill0(dimx,dimy,If2,NULL);

  FFFT2D(dimx,dimy,f1,If1,-1);
  FFFT2D(dimx,dimy,f2,If2,-1);

  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)      {
      buffR=f1[iy][ix]*f2[iy][ix]
	-If1[iy][ix]*If2[iy][ix];
      buffI=f1[iy][ix]*If2[iy][ix]
	+If1[iy][ix]*f2[iy][ix];
      f2[iy][ix]=norma*buffR;
      If2[iy][ix]=norma*buffI;
    }
  
  FFFT2D(dimx,dimy,f1,If1,1);
  FFFT2D(dimx,dimy,f2,If2,1);

  free_matrix2D(If1,dimy);
  free_matrix2D(If2,dimy);

 return OK;
}


/***************************************************************************/
int convuelto_vec_FFFT( int dimx, int dimy, double **f1x, double **f1y,
			 double **f2x, double **f2y) {
/***************************************************************************/
  int ix,iy;
  double **If1x,**If1y,**If2x,**If2y;
  double buffRx,buffIx;
  double buffRy,buffIy;
  double norma;

  norma=sqrt((double)(dimx*dimy));

  TrackNullAlloc( If1x=matrix2D(dimy,dimx) );
  TrackNullAlloc( If1y=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2x=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2y=matrix2D(dimy,dimx) );
  
  fill0(dimx,dimy,If1x,NULL);
  fill0(dimx,dimy,If1y,NULL);
  fill0(dimx,dimy,If2x,NULL);
  fill0(dimx,dimy,If2y,NULL);

  FFFT2D(dimx,dimy,f1x,If1x,-1);
  FFFT2D(dimx,dimy,f1y,If1y,-1);
  FFFT2D(dimx,dimy,f2x,If2y,-1);
  FFFT2D(dimx,dimy,f2y,If2y,-1);

  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)	{
      buffRx=f1x[iy][ix]*f2x[iy][ix]-If1x[iy][ix]*If2x[iy][ix]
	-(f1y[iy][ix]*f2y[iy][ix]-If1y[iy][ix]*If2y[iy][ix]);
      buffIx=f1x[iy][ix]*If2x[iy][ix]+If1x[iy][ix]*f2x[iy][ix]
	-(f1y[iy][ix]*If2y[iy][ix]+If1y[iy][ix]*f2y[iy][ix]);
      
      buffRy=f1x[iy][ix]*f2y[iy][ix]-If1x[iy][ix]*If2y[iy][ix]
	+f1y[iy][ix]*f2x[iy][ix]-If1y[iy][ix]*If2x[iy][ix];
      buffIy=f1x[iy][ix]*If2y[iy][ix]+If1x[iy][ix]*f2y[iy][ix]
	+f1y[iy][ix]*If2x[iy][ix]+If1y[iy][ix]*f2x[iy][ix];
      
      f2x[iy][ix]=norma*buffRx;
      f2y[iy][ix]=norma*buffRy;
      If2x[iy][ix]=norma*buffIx;
      If2y[iy][ix]=norma*buffIy;
    }

  FFFT2D(dimx,dimy,f1x,If1y,1);
  FFFT2D(dimx,dimy,f1y,If1y,1);
  FFFT2D(dimx,dimy,f2x,If2x,1);
  FFFT2D(dimx,dimy,f2y,If2y,1);

  free_matrix2D(If1x,dimy);
  free_matrix2D(If1y,dimy);
  free_matrix2D(If2x,dimy);
  free_matrix2D(If2y,dimy);

 return OK;
}


/***************************************************************************/
int deconvuelve_FFFT( int dimx, int dimy, double **f1, double **f2, 
		       double **salida) {
/***************************************************************************/
  int ix,iy;
  double **If1,**If2,**Isalida;
  double mod;

  TrackNullAlloc( If1=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2=matrix2D(dimy,dimx) );
  TrackNullAlloc( Isalida=matrix2D(dimy,dimx) );

  fill0(dimx,dimy,If1,NULL);
  fill0(dimx,dimy,If2,NULL);

  FFFT2D(dimx,dimy,f1,If1,-1);
  FFFT2D(dimx,dimy,f2,If2,-1);

  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)      {
      mod=f2[iy][ix]*f2[iy][ix]+If2[iy][ix]*If2[iy][ix];
      if(mod>1e-30)	{
	salida[iy][ix]=(f1[iy][ix]*f2[iy][ix]
			+If1[iy][ix]*If2[iy][ix])/mod;
	Isalida[iy][ix]=(-f1[iy][ix]*If2[iy][ix]
			 +If1[iy][ix]*f2[iy][ix])/mod;
      }	  else	    {
	salida[iy][ix]=0.;
	Isalida[iy][ix]=0.;
      }
    }

  FFFT2D(dimx,dimy,salida,Isalida,1);
  FFFT2D(dimx,dimy,f1,If1,1);
  FFFT2D(dimx,dimy,f2,If2,1);

  free_matrix2D(If1,dimy);
  free_matrix2D(If2,dimy);
  free_matrix2D(Isalida,dimy);

 return OK;
}


/***************************************************************************/
int deconvuelto_FFFT( int dimx, int dimy, double **f1, double **f2) {
/***************************************************************************/
  int ix,iy;
  double **If1,**If2;
  double buffR,buffI,mod;

  TrackNullAlloc( If1=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2=matrix2D(dimy,dimx) );

  fill0(dimx,dimy,If1,NULL);
  fill0(dimx,dimy,If2,NULL);

  FFFT2D(dimx,dimy,f1,If1,-1);
  FFFT2D(dimx,dimy,f2,If2,-1);

  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)      {
      mod=f1[iy][ix]*f1[iy][ix]+If1[iy][ix]*If1[iy][ix];
      if(mod>1e-30)	{
	buffR=(f1[iy][ix]*f2[iy][ix]
	       +If1[iy][ix]*If2[iy][ix])/mod;
	buffI=(f1[iy][ix]*If2[iy][ix]
	       -If1[iy][ix]*f2[iy][ix])/mod;
	f2[iy][ix]=buffR;
	If2[iy][ix]=buffI;
      }	  else	    {
	f2[iy][ix]=0.;
	If2[iy][ix]=0.;
      }
    }

  FFFT2D(dimx,dimy,f1,If1,1);
  FFFT2D(dimx,dimy,f2,If2,1);

  free_matrix2D(If1,dimy);
  free_matrix2D(If2,dimy);

 return OK;
}


/***************************************************************************/
int deconvuelto_vec_FFFT( int dimx, int dimy, double **f1x, double **f1y,
			   double **f2x, double **f2y) {
/***************************************************************************/
  int ix,iy;
  double **If1x,**If1y,**If2x,**If2y;
  double buffRx,buffIx;
  double buffRy,buffIy;
  double modR,modI,modd;
  double norma;

  norma=sqrt((double)(dimx*dimy));

  TrackNullAlloc( If1x=matrix2D(dimy,dimx) );
  TrackNullAlloc( If1y=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2x=matrix2D(dimy,dimx) );
  TrackNullAlloc( If2y=matrix2D(dimy,dimx) );

  fill0(dimx,dimy,If1x,NULL);
  fill0(dimx,dimy,If1y,NULL);
  fill0(dimx,dimy,If2x,NULL);
  fill0(dimx,dimy,If2y,NULL);

  FFFT2D(dimx,dimy,f1x,If1x,-1);
  FFFT2D(dimx,dimy,f1y,If1y,-1);
  FFFT2D(dimx,dimy,f2x,If2y,-1);
  FFFT2D(dimx,dimy,f2y,If2y,-1);

  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)	{
      
      modR=(f1x[iy][ix]*f1x[iy][ix]
	    -If1x[iy][ix]*If1x[iy][ix]
	    +f1y[iy][ix]*f1y[iy][ix]
	    -If1y[iy][ix]*If1y[iy][ix])*norma*norma;
      modI=(2.*f1x[iy][ix]*If1x[iy][ix]
	    +2.*f1y[iy][ix]*If1y[iy][ix])*norma*norma;
      
      modd=modR*modR+modI*modI;
      
      buffRx=(f1x[iy][ix]*f2x[iy][ix]
	      -If1x[iy][ix]*If2x[iy][ix]
	      +f1y[iy][ix]*f2y[iy][ix]
	      -If1y[iy][ix]*If2y[iy][ix])*norma;
      buffIx=(f1x[iy][ix]*If2x[iy][ix]
	      +If1x[iy][ix]*f2x[iy][ix]
	      +f1y[iy][ix]*If2y[iy][ix]
	      +If1y[iy][ix]*f2y[iy][ix])*norma;
      
      buffRy=(f1x[iy][ix]*f2y[iy][ix]
	      -If1x[iy][ix]*If2y[iy][ix]
	      -f1y[iy][ix]*f2x[iy][ix]
	      +If1y[iy][ix]*If2x[iy][ix])*norma;
      buffIy=(f1x[iy][ix]*If2y[iy][ix]
	      +If1x[iy][ix]*f2y[iy][ix]
	      -f1y[iy][ix]*If2x[iy][ix]
	      -If1y[iy][ix]*f2x[iy][ix])*norma;
      
      if(modd>1e-30)	    {
	f2x[iy][ix]=(buffRx*modR+buffIx*modI)/modd;
	If2x[iy][ix]=(buffIx*modR-buffRx*modI)/modd;
	
	f2y[iy][ix]=(buffRy*modR+buffIy*modI)/modd;
	If2y[iy][ix]=(buffIy*modR-buffRy*modI)/modd;	
      }  else	{
	f2x[iy][ix]=0.;
	If2x[iy][ix]=0.;
	f2y[iy][ix]=0.;
	If2y[iy][ix]=0.;
      }
      
    }

  FFFT2D(dimx,dimy,f1x,If1y,1);
  FFFT2D(dimx,dimy,f1y,If1y,1);
  FFFT2D(dimx,dimy,f2x,If2x,1);
  FFFT2D(dimx,dimy,f2y,If2y,1);

  free_matrix2D(If1x,dimy);
  free_matrix2D(If1y,dimy);
  free_matrix2D(If2x,dimy);
  free_matrix2D(If2y,dimy);

  return OK;
}



