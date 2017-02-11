/********************************************************/
/*                                                      */
/*	fft1d.c - Version del 21 de Septiembre, 2005      */
/*                                                      */
/********************************************************/


#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_alloc.h>
#include <utl_stats.h>
#include <utl_operator.h>

#include <filter.h>
#include <flt_fft1d.h>

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



/***************************************************************************/
void Fourier1D(int dimx, double *funcionR, double *funcionI, int signo) {
  /***************************************************************************/

#ifdef _PARSE_FILTER_PARAMETERS_
  IF(p_fil->flag_memory)
#else
    IF(FLAG_MEMORY)
#endif/* !_PARSE_FILTER_PARAMETERS_ */
    FFFT1D(dimx,funcionR,funcionI,signo);
  else 
    FFT1D(dimx,funcionR,funcionI,signo);
}


/***************************************************************************/
void FFT1D(int dim, double *funcionR, double *funcionI, int signo) {
/***************************************************************************/
  
  double tempR,tempI,wpasoR,wpasoI,wwR,wwI;
  int ix,je,mm,mmax,istep;
  
  je=1;
  for(ix=0;ix<dim;ix++)    {
    if(je>ix+1)	{
      tempR=funcionR[je-1];
      tempI=funcionI[je-1];
      funcionR[je-1]=funcionR[ix];
      funcionI[je-1]=funcionI[ix];
      funcionR[ix]=tempR;
      funcionI[ix]=tempI;
    }
    mm=dim/2;
    
    while((mm>1)&&(je>mm))	{
      je=je-mm;
      mm=mm/2;
    }
    je=je+mm;  
  }	

  mmax=1;
  while(dim>mmax)    {
    istep=2*mmax;
    wpasoR=cos(M_PI/((double) mmax));
    wpasoI=signo*sin(M_PI/((double) mmax));
    wwR=1.;
    wwI=0.;
    
    for(mm=1;mm<=mmax;mm++)	{
      for(ix=mm-1;ix<dim;ix+=istep)	{
	je=ix+mmax;
	C_mult(wwR,wwI,funcionR[je],funcionI[je],
	       &tempR,&tempI);
	funcionR[je]=funcionR[ix]-tempR;
	funcionI[je]=funcionI[ix]-tempI;
	funcionR[ix]=funcionR[ix]+tempR;
	funcionI[ix]=funcionI[ix]+tempI;
      }
      C_mult(wwR,wwI,wpasoR,wpasoI,&wwR,&wwI);
    }
    mmax=istep;
    
  }

  for(ix=0;ix<dim;ix++)  {
    funcionR[ix]=funcionR[ix]/sqrt((double)dim);
    funcionI[ix]=funcionI[ix]/sqrt((double)dim);
  }
 
}


/***************************************************************************/
int FFFT1D( int dim, double *funcionR, double *funcionI, int signo) {
  /***************************************************************************/
  double *workR,*workI;
  
  int inu;
  int dima,dimb,dimc;
  int i,i1,i2,i3;

  TrackNullAlloc( workR=(double*)calloc(dim,sizeof(double)) );
  TrackNullAlloc( workI=(double*)calloc(dim,sizeof(double)) );

  dima=1;
  dimb=dim;
  dimc=1;
  inu=1;

  while(dimb>1)    {
    dima=dimc*dima;
    dimc=2;
    while(dimb%dimc!=0) dimc++;
    
    dimb=dimb/dimc;
    
    if(inu==1) PFFT1D(dima,dimb,dimc,funcionR,funcionI,
		      workR,workI,signo);
    else PFFT1D(dima,dimb,dimc,workR,workI,funcionR,funcionI,
		signo);
    inu=1-inu;
  }

  if(inu==0)    
    for(i=0;i<dim;i++)	{
      funcionR[i]=workR[i];
      funcionI[i]=workI[i];
    }

  for(i=0;i<dim;i++)    {
    funcionR[i]=funcionR[i]/sqrt((double)dim);
    funcionI[i]=funcionI[i]/sqrt((double)dim);
  }
  
  free(workR);
  free(workI);

  return OK;
}


/***************************************************************************/
void PFFT1D( int dima, int dimb, int dimc, double *uinR, double *uinI, 
	     double *uoutR, double *uoutI, int sign) {
/***************************************************************************/
  double angle;
  double deltaR,deltaI;
  double omegaR,omegaI;
  double sumR,sumI;

  int ia,ib,ic,i,jcr,jc;

  angle=2.*M_PI/((double)(dima*dimc));
  omegaR=1.;
  omegaI=0.;

  deltaR=cos(angle);
  deltaI=((double)sign)*sin(angle);

  for(ic=0;ic<dimc;ic++)
    for(ia=0;ia<dima;ia++)  {
      for(ib=0;ib<dimb;ib++)    {
	i=ib+dimb*(dimc-1+ia*dimc);
	sumR=uinR[i];
	sumI=uinI[i];
	for(jcr=1;jcr<dimc;jcr++)		{
	  jc=dimc-1-jcr;
	  i=ib+dimb*(jc+ia*dimc);
	  C_mult(omegaR,omegaI,sumR,sumI,&sumR,&sumI);
	  sumR+=uinR[i];
	  sumI+=uinI[i];
	}
	i=ib+dimb*(ia+ic*dima);
	uoutR[i]=sumR;
	uoutI[i]=sumI;
      }
      C_mult(omegaR,omegaI,deltaR,deltaI,&omegaR,&omegaI);
    }

}


/* void convuelto1D( int dimx, double *f1, double *f2)
   {
   double *If1,*If2;
   double buffR,buffI;
   int xeff;
   int ix;
   
   xeff=redimensiona(dimx);
   If1=(double *)calloc(xeff,sizeof(double));
   If2=(double *)calloc(xeff,sizeof(double));
   
   Fourier1D(xeff,f2,If2,-1);
   Fourier1D(xeff,f1,If1,-1);
   for(ix=0;ix<xeff;ix++)
   {
   buffR=f1[ix]*f2[ix]-If1[ix]*If2[ix];
   buffI=f1[ix]*If2[ix]+If1[ix]*f2[ix];
   f2[ix]=buffR*sqrt(dimx);
   If2[ix]=buffI*sqrt(dimx);
   }
   Fourier1D(xeff,f1,If1,1);
   Fourier1D(xeff,f2,If2,1);
   
   free(If1);
   free(If2);
   }
*/


/***************************************************************************/
int convuelto1D( int dimx, double *f1, double *f2) {
/***************************************************************************/
    double *Rf1,*Rf2,*If1,*If2;
    double buffR,buffI;
    int xeff;
    int ix;
    
    xeff=redimensiona(dimx);
    TrackNullAlloc( Rf1=(double*)calloc(xeff,sizeof(double)) );
    TrackNullAlloc( Rf2=(double*)calloc(xeff,sizeof(double)) );
    TrackNullAlloc( If1=(double*)calloc(xeff,sizeof(double)) );
    TrackNullAlloc( If2=(double*)calloc(xeff,sizeof(double)) );
    
    copy1D(dimx,f1,Rf1,NULL);
    copy1D(dimx,f2,Rf2,NULL);
    
    Fourier1D(xeff,Rf2,If2,-1);
    Fourier1D(xeff,Rf1,If1,-1);

    for(ix=0;ix<xeff;ix++)    {
      buffR=Rf1[ix]*Rf2[ix]-If1[ix]*If2[ix];
      buffI=Rf1[ix]*If2[ix]+If1[ix]*Rf2[ix];
      Rf2[ix]=buffR*sqrt(dimx);
      If2[ix]=buffI*sqrt(dimx);
    }

    Fourier1D(xeff,Rf1,If1,1);
    Fourier1D(xeff,Rf2,If2,1);

    copy1D(dimx,Rf2,f2,NULL);

    free(Rf1);
    free(Rf2);
    free(If1);
    free(If2);

    return OK;
}


/***************************************************************************/
int convuelve1D( int dimx, double *f1, double *f2, double *salida) {
/***************************************************************************/
    double *Rf1,*Rf2,*If1,*If2;
    double buffR,buffI;
    int xeff;
    int ix;
    
    xeff=redimensiona(dimx);
    TrackNullAlloc( Rf1=(double *)calloc(xeff,sizeof(double)) );
    TrackNullAlloc( Rf2=(double *)calloc(xeff,sizeof(double)) );
    TrackNullAlloc( If1=(double *)calloc(xeff,sizeof(double)) );
    TrackNullAlloc( If2=(double *)calloc(xeff,sizeof(double)) );
    
    copy1D(dimx,f1,Rf1,NULL);
    copy1D(dimx,f2,Rf2,NULL);

    Fourier1D(xeff,Rf2,If2,-1);
    Fourier1D(xeff,Rf1,If1,-1);

    for(ix=0;ix<xeff;ix++)    {
      buffR=Rf1[ix]*Rf2[ix]-If1[ix]*If2[ix];
      buffI=Rf1[ix]*If2[ix]+If1[ix]*Rf2[ix];
      Rf2[ix]=buffR*sqrt(dimx);
      If2[ix]=buffI*sqrt(dimx);
    }

    Fourier1D(xeff,Rf1,If1,1);
    Fourier1D(xeff,Rf2,If2,1);

    copy1D(dimx,Rf2,salida,NULL);

    free(Rf1);
    free(Rf2);
    free(If1);
    free(If2);

  return OK;
}


/***************************************************************************/
void FFT(int dim, double *funcionR, double *funcionI, int signo) {
/***************************************************************************/

  double tempR,tempI,wpasoR,wpasoI,wwR,wwI;
  int ix,je,mm,mmax,istep;
  
  je=1;
  for(ix=0;ix<dim;ix++)    {
    if(je>ix+1)	{
      tempR=funcionR[je-1];
      tempI=funcionI[je-1];
      funcionR[je-1]=funcionR[ix];
      funcionI[je-1]=funcionI[ix];
      funcionR[ix]=tempR;
      funcionI[ix]=tempI;
    }
    mm=dim/2;
    while((mm>1)&&(je>mm))      {
      je=je-mm;
      mm=mm/2;
    }
    je=je+mm;
  }	
  
  mmax=1;
  while(dim>mmax)    {
    istep=2*mmax;
    wpasoR=cos(M_PI/((double) mmax));
    wpasoI=signo*sin(M_PI/((double) mmax));
    wwR=1.;
    wwI=0.;
    
    for(mm=1;mm<=mmax;mm++)	{
      for(ix=mm-1;ix<dim;ix+=istep)	{
	je=ix+mmax;
	C_mult(wwR,wwI,funcionR[je],funcionI[je],
	       &tempR,&tempI);
	funcionR[je]=funcionR[ix]-tempR;
	funcionI[je]=funcionI[ix]-tempI;
	funcionR[ix]=funcionR[ix]+tempR;
	funcionI[ix]=funcionI[ix]+tempI;
      }
      C_mult(wwR,wwI,wpasoR,wpasoI,&wwR,&wwI);
    }
    mmax=istep;
  }
  
  for(ix=0;ix<dim;ix++)    {
    funcionR[ix]=funcionR[ix]/sqrt((double)dim);
    funcionI[ix]=funcionI[ix]/sqrt((double)dim);
  }

}

