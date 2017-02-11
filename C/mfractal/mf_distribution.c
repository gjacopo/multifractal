#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mfractal.h>
#include <mf_distribution.h>

#ifdef _PARSE_FRACTAL_PARAMETERS_
#include <mf_parse.h>
extern ParFRAC *p_frac;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/


/***************************************************************************/
double genera_h( double* prob_levi ) {
  /***************************************************************************/
  double h;
  double beta,s;
  double prob;
  int ip;
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  switch(p_frac->mf_type) 
#else
    switch(TYPE_MFSIM) 
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
      {
      case TYPLOGPOISSON:   /* Log-Poisson */
#ifdef _PARSE_FRACTAL_PARAMETERS_
	beta = 1+ p_frac->hinf/ p_frac->codinf;
	s = p_frac->codinf * log(2.);
	h = p_frac->hinf - poisson(s) * log(beta) / log(2.);
#else
	beta = 1 + HINF / CODINF;
	s = CODINF * log(2.);
	h = HINF - poisson(s) * log(beta) / log(2.);
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	break;
	
      case TYPLOGNORMAL: /* Log-Normal */
#ifdef _PARSE_FRACTAL_PARAMETERS_
	h= p_frac->mu + p_frac->sigma * normal_standard(p_frac->nbox,p_frac->tch) / sqrt(log(2.));
#else
	h = MU + SIGMA * normal_standard(NBOX,TCH) / sqrt(log(2.));
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	break;
	
      case TYPLOGLEVI: /* Log-Levi */
#ifdef _PARSE_FRACTAL_PARAMETERS_
	h = p_frac->mu + 
	  p_frac->sigma * levi_standard(prob_levi,p_frac->nbox,p_frac->tch) * 
	  pow(log(2.), -1./p_frac->alpha);
#else
	h = MU + 
	  SIGMA * levi_standard(prob_levi,NBOX,TCH) * 
	  pow(log(2.), -1./ALPHA);
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	break;
	
      case TYPBINOMIAL: /* Binomial */
#ifdef _PARSE_FRACTAL_PARAMETERS_
	h = binomial(p_frac->hinf, p_frac->h1);
#else
	h = binomial(HINF, H1);
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	break;
	
      case TYPMONOFRACTAL: /* Monofractal */
#ifdef _PARSE_FRACTAL_PARAMETERS_
	prob = pow(2., -p_frac->codinf);
	h = monofractal(p_frac->hinf, prob);
#else
	prob = pow(2., -CODINF);
	h = monofractal(HINF, prob);
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	break;
	
      default:
	h=0.;
	break;
      }
  
  return h;
} // end of genera_h


/***************************************************************************/
double theoretical_Dh( double h ) {
/***************************************************************************/
  double Dh_th;
  double beta,omega,y;
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  switch(p_frac->mf_type) 
#else
    switch(TYPE_MF) 
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
      {
      case TYPLOGPOISSON:   /* Log-Poisson */
#ifdef _PARSE_FRACTAL_PARAMETERS_
	beta = 1+ p_frac->hinf/ p_frac->codinf;
	if( (omega=-(h-p_frac->hinf)/(p_frac->codinf*log(beta))) > 1.e-30) 
	  Dh_th = (double)p_frac->dim_space + p_frac->codinf*(omega-omega*log(omega)-1.);
#else
	beta = 1 + HINF/CODINF;
	if( (omega=-(h-HINF)/(CODINF*log(beta))) > 1.e-30) 
	  Dh_th = (double)DSPACE + CODINF*(omega-omega*log(omega)-1.);
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	else Dh_th = 0.;
	break;
      case TYPLOGNORMAL: /* Log-Normal */
#ifdef _PARSE_FRACTAL_PARAMETERS_
	Dh_th = (double)p_frac->dim_space - 
	  (h-p_frac->mu)*(h-p_frac->mu)/(2.*p_frac->sigma*p_frac->sigma);
#else
	Dh_th = (double)DSPACE - (h-MU)*(h-MU)/(2.*SIGMA*SIGMA);
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	break;
      case TYPLOGLEVI: /* Log-Levi */
#ifdef _PARSE_FRACTAL_PARAMETERS_
	Dh_th = (double)p_frac->dim_space - 
	  pow(fabs(h-p_frac->mu)/p_frac->sigma,p_frac->alpha)/p_frac->alpha;
#else
	Dh_th = (double)DSPACE - pow(fabs(h-MU)/SIGMA,ALPHA)/ALPHA;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	break;
      case TYPBINOMIAL: /* Binomial */
#ifdef _PARSE_FRACTAL_PARAMETERS_
	y = (p_frac->h-p_frac->hinf) / (p_frac->h1-p_frac->hinf);
	if((y<0.01) || (y>0.99))  Dh_th = 0.;
	else Dh_th = (double)p_frac->dim_space -1. -(y*log(y)+(1-y)*log(1-y))/log(2.);
#else
	y = (h-HINF) / (H1-HINF);
	if((y<0.01) || (y>0.99))  Dh_th = 0.;
	else Dh_th = (double)DSPACE -1. -(y*log(y)+(1-y)*log(1-y))/log(2.);
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	break;
      case TYPMONOFRACTAL: /* Monofractal */
	break;
      default:
	Dh_th = 0.; // (double)DSPACE
	break;
      }
  
  return Dh_th;
} // end of theoretical_Dh


/***************************************************************************/
int theoretical_Deltah( double *hmin, double *hmax ) {
/***************************************************************************/
  const double steph=0.01;
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  switch(p_frac->mf_type) 
#else
    switch(TYPE_MF) 
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
      {
      case TYPLOGPOISSON:   /* Log-Poisson */
#ifdef _PARSE_FRACTAL_PARAMETERS_
	*hmin = p_frac->hinf;
#else
	*hmin = HINF;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	/* As the equation is trancendent, we look for the second zero-crossing 
	 * in an iterative way 
	 */
	for( *hmax=0.; theoretical_Dh(*hmax)>0.; *hmax+=steph );
	*hmax -= steph;
	break;

      case TYPLOGNORMAL: /* Log-Normal */
#ifdef _PARSE_FRACTAL_PARAMETERS_
	*hmin = p_frac->sigma * sqrt(2.*(double)p_frac->dim_space);
	*hmax = p_frac->mu + (*hmin);
#else
	*hmin = SIGMA * sqrt(2.*(double)DSPACE);
	*hmax = MU + (*hmin);
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	*hmin = *hmax - (*hmin) * 2.;
	    break;

      case TYPLOGLEVI: /* Log-Levi */
#ifdef _PARSE_FRACTAL_PARAMETERS_
	*hmin = p_frac->sigma * pow(p_frac->alpha*(double)p_frac->dim_space,1./p_frac->alpha);
	*hmax = p_frac->mu + (*hmin);
#else
	*hmin = SIGMA * pow(ALPHA*(double)DSPACE,1./ALPHA);
	*hmax = MU + (*hmin);
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	*hmin = *hmax - (*hmin) * 2.;
	break;
	
      case TYPBINOMIAL: /* Binomial */
#ifdef _PARSE_FRACTAL_PARAMETERS_
	*hmin=p_frac->hinf,     *hmax=p_frac->h1;
#else
	*hmin=HINF,     *hmax=H1;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	break;
      
      case TYPMONOFRACTAL: /* Monofractal */
	break;
      
      default:
	*hmin=-1.,	*hmax=1.;
	break;
      }
  
  return OK;
} // end of theoretical_Deltah


/***************************************************************************/
double poisson( double lambda ) {
/***************************************************************************/
  double tirada,prob,weight;
  int ip;

  tirada = ((double)random()) / RAND_PER;
  weight = exp(-lambda);
  ip = 0;

  for( prob=weight; ((prob<tirada)&&(ip<100)); ip++,prob+=weight )
    weight *= lambda / ((double)(ip+1));
  
  return((double)ip);
} // end of poisson


/***************************************************************************/
double normal_standard( int nbox, double tch ) {
/***************************************************************************/
  double norma,dx;
  double tirada,prob;
  double x;
  int ib; 
 
  dx = 2. * tch / ((double)nbox);
  norma = 0.;
  x = -tch;
  for( ib=0; ib<nbox; ib++,x+=dx ) norma += exp(-0.5*x*x);
  norma = 1. / norma;

  tirada = ((double)random()) / RAND_PER; 
  x = -tch;
  prob = 0.;
  for( ib=0; (ib<nbox)&&(prob<=tirada); ib++,x+=dx )
    prob += norma * exp(-0.5*x*x);
 
  return x;
} // end of normal_standard


/***************************************************************************/
int genera_levi( double*prob_levi, double alpha, int nbox, double tch ) {
  /***************************************************************************
   * To keep a copy of a Levi distribution of parameter alpha
   ***************************************************************************/
  double *prob_leviI;
  double x,dx,norma;
  int ix;

  TrackNullAlloc( prob_leviI=(double*)calloc(nbox,sizeof(double)) );

  for(ix=0;ix<nbox;ix++)    {
    x = (double)ix;
    if(ix>nbox/2) x -= (double)nbox;
    dx = fabs(x*M_PI/(2.*tch));
    if(dx>1e-30) prob_levi[ix] = exp(-alpha*pow(dx,alpha));
    prob_leviI[ix] = prob_levi[ix]*sin(M_PI*x);
    prob_levi[ix] = prob_levi[ix]*cos(M_PI*x);
  }

  FFFT1D(nbox,prob_levi,prob_leviI,-1);

  for(ix=0,norma=0.;ix<nbox;ix++) norma = fMin(norma,prob_levi[ix]);
  for(ix=0;ix<nbox;ix++) prob_levi[ix] -= norma;

  for(ix=0,norma=0.;ix<nbox;ix++) norma += prob_levi[ix];
  for(ix=0;ix<nbox;ix++) prob_levi[ix] /= norma;

  free(prob_leviI);

  return OK;
} // end of genera_levi


/***************************************************************************/
double levi_standard( double *prob_levi, int nbox, double tch ) {
 /***************************************************************************/
  double norma;
  double tirada,prob,x,dx;
  int ib; 

  tirada=((double)random())/RAND_PER;

  x = -tch;
  dx = 2.*tch/((double)nbox);
  prob = 0.;
  for(ib=0;(ib<nbox)&&(prob<tirada);ib++,x+=dx)
    prob += prob_levi[ib];

  return x;
} // end of levi_standard


/***************************************************************************/
double monofractal( double h0, double prob ) {
/***************************************************************************/
  double tirada;
  double h;

  acprob += prob;
  tirada = ((double)random())/RAND_PER;
  if(tirada <= acprob)    {
    h = h0;
    acprob -= 1.;
  }  else h=1e30;

  return h;
} // end of monofractal


/***************************************************************************/
double binomial( double h0, double h1 ) {
  /***************************************************************************/
  double tirada;
  double h;
  
  tirada = ((double)random())/RAND_PER;
  if(tirada < 0.5) h = h0;
  else             h = h1;

  return h;
} // end of binomial
 
