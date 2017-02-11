/********************************************************/
/*                                                      */
/*              mf_simulation.c                */
/*            Version: 18 de Octubre, 2004              */
/*                                                      */
/********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Personnal libraries */

#include <utils.h>
#include <utl_alloc.h>
#include <utl_operator.h>
#include <utl_stats.h>

#include <mfractal.h>
#include <mf_simulation.h>

extern int flag_verbose;

extern double *prob_levi; 
extern double m_th[2]; 
extern double acprob;  

/* */

#ifdef _PARSE_FRACTAL_PARAMETERS_
#include <mf_parse.h>
extern ParFRAC *p_frac;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */



#ifdef _PARSE_FRACTAL_PARAMETERS_
/* See file mfractal.h for default definitions and parameters adjustment
 * when the flag _PARSE_FRACTAL_PARAMETERS_ is not present */
/***************************************************************************/
int prepara_multifractal(  ) {
  /***************************************************************************/
  
  double val,beta;
  /* these parameters may be modified within this function */
  double mu=p_frac->mu, sigma=p_frac->sigma;
  double alpha=p_frac->alpha;
  double hinf=p_frac->hinf, h1=p_frac->h1;
  
  switch(p_frac->mf_type)	{

  case TYPLOGPOISSON: // 0: Log-Poisson
    if(p_frac->inv_trans)  // Constraint of translational invariance
      IFVERBOSE Warning("Log-Poisson MF are scale invariant, flag 'inv_trans' ignored");
    if(p_frac->codinf > 1.e-30) beta = 1.+ hinf / p_frac->codinf;
    else beta=1.;
    if(p_frac->dim_space == DIM1D) val=exp(1.);
    else val=3.59;
    m_th[0] = hinf;
    m_th[1] = hinf - val * log(beta) * p_frac->codinf;
    break;

  case TYPLOGNORMAL: // 1: Log-Normal
    if(p_frac->inv_trans)   {
      p_frac->mu = mu = sigma * sigma / 4.;
      IFVERBOSE	WarningV("Mean singularity changed to %f",mu);
    }
    m_th[0] = mu - p_frac->tch*sigma;
    m_th[1] = mu + p_frac->tch*sigma; 

    break; 

  case TYPLOGLEVI: // 2: Log-Levi
    if(alpha < 0.01)	{
      /* Change the values of the analyzing parameters */
      p_frac->alpha = alpha = 0.01;
      IFVERBOSE Warning("Too small alpha value. Quantizing...");
    }
    if(p_frac->inv_trans)   {
      if(fabs(alpha-1.) > 1.e-30)	    {
	p_frac->mu = mu = pow(sigma/alpha, 1.+1./(alpha-1.)) * (alpha-1.);
	IFVERBOSE WarningV("Mean singularity changed to %f",mu);
      }
      ELSE 
	IFVERBOSE Warning("Scale invariance impossible to implement for this MF");
    }
    m_th[0] = mu - p_frac->tch*sigma;
    m_th[1] = mu + p_frac->tch*sigma;
    break;

  case TYPBINOMIAL: // 3: Binomial
    if(h1 <= hinf)      {
      p_frac->h1 = h1 = hinf + p_frac->dh;
      IFVERBOSE WarningV("Unacceptable greater exponent; changing it to %f", h1);
    }
    if(p_frac->inv_trans)   
      if(hinf > 1.)	{
	p_frac->h1 = h1 = -log(1-pow(0.5,hinf));
	IFVERBOSE WarningV("Maximum singularity changed to %0.2f",h1);
      } else	{
	if((h1>0.) && (h1<1.))	  {
	  p_frac->h1 = h1 = 0.75;
	  IFVERBOSE WarningV("Maximum singularity changed to %0.2f",h1);
	}
	p_frac->hinf = hinf = -log(1-pow(0.5,h1));
	IFVERBOSE WarningV("Minimum singularity changed to %0.2f",hinf);
      }  
    m_th[0] = hinf;
    m_th[1] = h1;
    break;
    
  case TYPMONOFRACTAL:  // 4: Monofractal
    if(p_frac->inv_trans)   {
      p_frac->hinf = hinf = -p_frac->codinf;
      IFVERBOSE WarningV("Least singularity changed to %f",hinf);
    }
    m_th[0] = hinf - p_frac->dh; // A conventional error around the value is accepted 
    m_th[1] = hinf + p_frac->dh;
    break;
  default:
    m_th[0] = -1.;
    m_th[1] = 2.;
    break;
  }
 
  return OK;
}
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/


/***************************************************************************/
void genera_multifractal( int leff, double sc0, int wav, int ord_der, double *prob_exp, 
			   double **signal) {
   /***************************************************************************/

#ifdef _PARSE_FRACTAL_PARAMETERS_
  if(p_frac->dim_space == DIM1D) 
#else
    if(DSPACE == DIM1D) 
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
      genera1D_multifractal(leff,sc0, wav,ord_der,prob_exp,signal[0]);
  else 
    genera2D_multifractal(leff,sc0, wav,ord_der,prob_exp,signal);
} // end of genera_multifractal


/***************************************************************************/
int genera1D_multifractal( int leff, double sc0, int wav, int ord_der, 
			   double *prob_exp, double *serie) {
  /***************************************************************************/
  double *wavesig,*aux,*serie2;
  double *lalpha0;
  int dimalpha,linfalpha,Nwav;
  int dima,bl,db;
  int ix,ip,iwav;

  double block=
#ifdef _PARSE_FRACTAL_PARAMETERS_
    pow(2.,(double)p_frac->outres);
#else
  pow(2.,(double)OUTRES);
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/

  Nwav = adimensiona(leff) - 1;
  dimalpha = leff - 1;

  TrackNullAlloc( lalpha0=(double*)calloc(dimalpha,sizeof(double)) );
  TrackNullAlloc( wavesig=(double*)calloc(leff,sizeof(double)) );
  TrackNullAlloc( aux=(double*)calloc(leff,sizeof(double)) );
  TrackNullAlloc( serie2=(double*)calloc(leff,sizeof(double)) );
  
  genera1D_alpha( dimalpha, Nwav, prob_exp, lalpha0 );
  
  for( iwav=0,dima=1,linfalpha=0; iwav<=Nwav; iwav++,dima*=2 )   {
    wavelet1D_define_unit(leff,(double)dima,wav,ord_der,sc0,wavesig);
    fill1D0(leff,aux,NULL);
    bl = leff / dima;
    db = bl / 2;
    for(ix=0;ix<dima;ix++) aux[bl*ix+db] = lalpha0[linfalpha+ix];
    convuelto1D(leff,wavesig,aux);
    for(ix=0;ix<leff;ix++) serie2[ix] += aux[ix];
    linfalpha += dima;
  }
  
  coarseres1D( leff, block, serie2, serie );
  
  free(lalpha0);
  free(wavesig);
  free(serie2);
  free(aux);

  return OK;
} // end of genera1D_multifractal


/***************************************************************************/
int genera2D_multifractal( int leff, double sc0, int wav, int ord_der, 
			   double *prob_exp, double **image) {
  /***************************************************************************/
  double **wavesig,**aux,**image2;
  double *lalpha0;
  int dimalpha,linfalpha,Nwav;
  int dima,bl,db;
  int ix,iy,iwav;

  double block=
#ifdef _PARSE_FRACTAL_PARAMETERS_
    pow(2.,(double)p_frac->outres);
#else
  pow(2.,(double)OUTRES);
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
  
  Nwav=adimensiona(leff)-1;
  dimalpha=(leff*leff-1)/3;

  TrackNullAlloc( lalpha0=(double*)calloc(dimalpha,sizeof(double)) );  
  TrackNullAlloc( wavesig=matrix2D(leff,leff) );
  TrackNullAlloc( aux=matrix2D(leff,leff) );
  TrackNullAlloc( image2=matrix2D(leff,leff) );

  genera2D_alpha( dimalpha, Nwav, prob_exp, lalpha0 );
  
  for( iwav=0,dima=1,linfalpha=0; iwav<=Nwav; iwav++,dima*=2 )  {
    wavelet2D_define_unit( leff, leff, (double)dima, wav, ord_der, sc0, wavesig );
    fill0( leff, leff, aux, NULL );
    bl = leff / dima;
    db = bl / 2;
    for( iy=0; iy<dima; iy++ )     
      for( ix=0; ix<dima; ix++ )   
	aux[bl*iy+db][bl*ix+db] = lalpha0[linfalpha+ix+iy*dima];
    convuelto2D( leff, leff, wavesig, aux );
    op_add( leff, leff, aux, image2, NULL );
    
    linfalpha += dima * dima;
  }
  coarseres( leff, leff, block, image2, image );
  
  free(lalpha0);
  free_matrix2D(wavesig,leff);
  free_matrix2D(image2,leff);
  free_matrix2D(aux,leff);

  return OK;
} // end of genera2D_multifractal


/***************************************************************************/
void genera_alpha( int dimeta, int Nwav, double *prob_exp, double *alpha) {
  /***************************************************************************/

#ifdef _PARSE_FRACTAL_PARAMETERS_
  if(p_frac->dim_space == DIM1D) 
#else
    if(DSPACE == DIM1D) 
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
    genera1D_alpha(dimeta,Nwav,prob_exp,alpha);
  else 
    genera2D_alpha(dimeta,Nwav,prob_exp,alpha);
} // end of genera_alpha


/***************************************************************************/
int genera1D_alpha( int dimeta, int Nwav, double *prob_exp, double *lalpha0 ) {
  /***************************************************************************/
  double *lalpha,*llalpha;
  double llalphaexp;
  double heff;
  double signo;
  int linfalpha,linfalpha0;
  int dima,dima0,Nit;
  int ix,ix0,iwav;
  int ip,it;
  double *prob_levi=NULL;
  
  int mftype, nbox;
  double alpha,dens,tch;
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
    nbox = p_frac->nbox, tch = p_frac->tch;
    mftype = p_frac->type_mfsim;
    alpha = p_frac->alpha, dens = p_frac->dens;
#else
  nbox = NBOX, tch = TCH;
  mftype = TYPE_MFSIM;
  alpha=ALPHA, dens = DENS;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/

  TrackNullAlloc( alpha=(double*)calloc(dimeta,sizeof(double)) );
  TrackNullAlloc( llalpha=(double*)calloc(dimeta,sizeof(double)) );

  if(mftype == TYPLOGLEVI) { // prepare the log-levi distribution
    TrackNullAlloc( prob_levi=(double*)calloc(nbox,sizeof(double)) );
    genera_levi( prob_levi, alpha, nbox, tch );
  }

  if(mftype == TYPMONOFRACTAL) {  // Monofractal to be calculated
    Nit=(int)(dens*pow(4.,(double)Nwav+1.));
    if(Nit<1) Nit=1;		  
  } else Nit=1;
  
  for( it=0; it<Nit; it++ )    {    
    lalpha[0] = 1.;
    llalpha[0] = 0.;
    linfalpha = 0;
    dima = 1;
    
    for( iwav=1; iwav<=Nwav; iwav++ )	{      
      /*         Update of parameters              */ 
      linfalpha0 = linfalpha;
      linfalpha += dima;
      dima0 = dima;
      dima *= 2;
      
      /*       Generation                 */
      for( ix=0; ix<dima; ix++ )	{
	ix0 = ix / 2;
	
	/*        Producing the accompanying multiplicative factor     */	
	if((signo=((double)random())/RAND_PER) < 0.5) signo = -1.;
	else                                          signo = 1.;

	lalpha[linfalpha+ix] = signo * sqrt(0.5) * lalpha[linfalpha0+ix0];
	
	/*       Producing the logarithm of alphas            */	
	llalpha[linfalpha+ix] = genera_h(prob_levi) + llalpha[linfalpha0+ix0];	

	/*      Accumulating alphas      */	
	lalpha0[linfalpha+ix] +=
	  lalpha[linfalpha+ix] * exp(-log(2.) * llalpha[linfalpha+ix]);
	
	/*       Accounting for experimental distribution of h     */	
	if((prob_exp!=NULL) && (iwav==Nwav))		{
	  if(mftype != TYPMONOFRACTAL)  
	    llalphaexp = -log(exp(-log(2.)*llalpha[linfalpha+ix])) / log(2.);
	  else 
	    llalphaexp = llalpha[linfalpha+ix];
	  heff = llalphaexp / ((double)Nwav);
	  if((mftype==TYPMONOFRACTAL) && (heff>1e5)) 
	    ip = nbox - 1;
	  else 
	    ip = (int)(((double)nbox) * (heff-m_th[0]) / (m_th[1]-m_th[0]));
	  
	  if(ip < 0) ip=0;
	  if(ip > nbox-1) ip=nbox-1;
	  prob_exp[ip]+=1.;
	}	
      }
      
    } //     End of scale loop
  } //     End of iteration loop 
  
  free(llalpha);
  free(lalpha);
  if(mf_type == TYPLOGLEVI) free(prob_levi);

  return OK;
} // end of genera1D_alpha


/***************************************************************************/
int genera2D_alpha( int dimeta, int Nwav, double *prob_exp, double *lalpha0 ) {
/***************************************************************************/
  double *lalpha,*llalpha;
  double llalphaexp;
  double heff;
  double signo;
  int linfalpha,linfalpha0;
  int dima,dima0;
  int ix,iy,ix0,iy0,iwav;
  int ip,it,Nit;
  double *prob_levi=NULL;

  int mftype, nbox;
  double alpha,dens,tch;
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
    nbox = p_frac->nbox;
    mftype = p_frac->type_mfsim;
    alpha = p_frac->alpha, dens = p_frac->dens;
    tch = p_frac->tch;
#else
  nbox = NBOX;
  mftype = TYPE_MFSIM;
  alpha=ALPHA, dens = DENS;
  tch = TCH;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/

  TrackNullAlloc( lalpha=(double*)calloc(dimeta,sizeof(double)) );
  TrackNullAlloc( llalpha=(double*)calloc(dimeta,sizeof(double)) );

  if(mftype == TYPLOGLEVI) { // prepare the log-levi distribution
    TrackNullAlloc( prob_levi=(double*)calloc(nbox,sizeof(double)) );
    genera_levi( prob_levi, alpha, nbox, tch );
  }
  
  if(mftype == TYPMONOFRACTAL) { // Monofractal to be calculated
    Nit=(int)(dens*pow(4.,(double)Nwav+1.));
    if(Nit < 1) Nit = 1;		  
  } else Nit = 1;
  
  for( it=0; it<Nit; it++ )    {    
    lalpha[0] = 1.;
    llalpha[0] = 0.;
    linfalpha = 0;
    dima = 1;
    
    for( iwav=1; iwav<=Nwav; iwav++ )	{      
      /*         Update of parameters              */       
      linfalpha0 = linfalpha;
      linfalpha += dima*dima;
      dima0 = dima;
      dima *= 2;
      /*       Generation                 */      
      for( iy=0; iy<dima; iy++ )	    {
	iy0 = iy / 2;

	for( ix=0; ix<dima;ix++ )	  {
	  ix0 = ix / 2;
	  
	  /*        Producing the accompanying multiplicative factor     */	  
	  if((signo=((double)random())/RAND_PER) < 0.5) signo = -1.;
	  else                                          signo = 1.;
	  
	  lalpha[linfalpha+ix+iy*dima] =
	    signo * lalpha[linfalpha0+ix0+iy0*dima0];
	  
	  /*       Producing the logarithm of alphas            */	  
	  llalpha[linfalpha+ix+iy*dima] = genera_h(prob_levi) +
	    llalpha[linfalpha0+ix0+iy0*dima0];
	  
	  /*      Accumulating alphas      */
	  lalpha0[linfalpha+ix+iy*dima] += lalpha[linfalpha+ix+iy*dima] *
	    exp(-log(2.) * llalpha[linfalpha+ix+iy*dima]);
	  
	  /*       Accounting for experimental distribution of h     */
	  if((prob_exp!=NULL) && (iwav==Nwav)) {
	    if(mftype != TYPMONOFRACTAL) 
	      llalphaexp = - log( exp(-log(2.)*llalpha[linfalpha+ix+iy*dima])) /
		log(2.); 
	    else 
	      llalphaexp = llalpha[linfalpha+ix+iy*dima];
	    heff = llalphaexp / ((double)Nwav);
	    if((mftype==TYPMONOFRACTAL) && (heff>1.e5)) 
	      ip = nbox-1;
	    else 
	      ip = (int)(((double)nbox) *(heff-m_th[0]) / (m_th[1]-m_th[0]));
	    if(ip < 0) ip = 0;
	    if(ip > nbox-1) ip = nbox-1;
	    prob_exp[ip] += 1.;
	  }  	  
	}	
      }	 
    } //     End of scale loop
    
  } //     End of iteration loop 
  
  free(lalpha);
  free(llalpha);
  if(mf_type == TYPLOGLEVI) free(prob_levi);
  
  return OK;
} // end of genera2D_alpha


/***************************************************************************/
void genera_base( int leff, char *base ) {
/***************************************************************************/
  char aux[MAXCHARLENGTH];

  int dim_space, ndata;
  double mu, sigma, alpha;
  double hinf, codinf, h1;

#ifdef _PARSE_FRACTAL_PARAMETERS_
  dim_space = p_frac->dim_space,  ndata = p_frac->ndata;
  mu = p_frac->mu, sigma = p_frac->sigma, alpha = p_frac->alpha;
  hinf = p_frac->hinf, codinf = p_frac->codinf, h1 = p_frac->h1;
#else
  dim_space =  DSPACE, ndata = NDATA;
  mu = MU, sigma = SIGMA, alpha = ALPHA;
  hinf = HINF, codinf = CODINF, h1 = H1;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
  
  IFVERBOSE WarningVV("Analyzing for %d input %dD signals ",ndata,dim_space);
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  switch(p_frac->type_mfsim)
#else
	 switch(TYPE_MFSIM)
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
      { 
      case TYPLOGPOISSON: // Log-Poisson
	sprintf(base,"Log-Poisson.h%0.2f-coD%0.2f",hinf,codinf);
	IFVERBOSE 
	  WarningVV("of the type Log-Poisson; hinf= %0.2f, Dinf= %0.2f\n",
		    hinf,(double)dim_space-codinf);
	break;
      case TYPLOGNORMAL: // Log-Normal
	sprintf(base,"Log-Normal.mean%0.2f-sigma%0.2f",mu,sigma);
	IFVERBOSE 
	  WarningVV("of the type Log-normal; mean= %0.2f, sigma= %0.2f\n",
		    mu,sigma);
	break;
      case TYPLOGLEVI: // Log-Levi
	sprintf(base,"Log-Levi_mean%0.2f-sigma%0.2f-alpha%0.2f",mu,sigma,alpha);
	IFVERBOSE 
	  WarningVV("of the type Log-Levi; mean= %0.2f, sigma= %0.2f, alpha= %0.2f\n",
		    mu,sigma,alpha);
	break;
      case TYPBINOMIAL:  // Binomial
	sprintf(base,"Binomial_h0%0.2f-h1%0.2f",hinf,h1);
	if(verbose)
	  IFVERBOSE WarningVV("of the type binomial; h0= %0.2f, h1= %0.2f\n",
			      hinf,h1);
	break;
      case TYPMONOFRACTAL: // Monofractal
	sprintf(base,"Monomf_h0%0.2f-cod%0.2f",hinf,codinf);
	IFVERBOSE WarningVV("of the type monofractal; h0= %0.2f, D0= %0.2f\n",
			    hinf,(double)dim_space-codinf);
	break;
      default:
	sprintf(base,"bugged");
	break;
      }
  
  if(leff>0) {
    sprintf(aux,"-size%d",leff);
    strcat(base,aux);
    IFVERBOSE WarningV("and size of %d points\n",leff);
  }
  
  if(dim_space == DIM1D) strcat(base,"_1D");
  else                   strcat(base,"_2D");
  
  return OK;
} // end of genera_base

