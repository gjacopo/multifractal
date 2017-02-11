
/* ===================================
** mf_analysis.c
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
#include <mf_ghgmwp.h>
#include <mf_moments.h>
#include <mf_wtmm1d.h>

#include <mf_simulation.h>
#include <mf_estimation.h>
#include <mf_inout.h>

#include <mf_analysis.h>

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



/***************************************************************************/
int create_multifractal( int nseries, char *base ) {
/***************************************************************************/

  /*	QUANTITIES TO BE COMPUTED	*/
  double **signal;
  
  /*	AUXILIAR VARIABLES	*/
  int Nh_qth;
  double *h_qth, *Dh_qth, *errDh_qth;
  int block,leff,dimy;
  int in,ip,ir;
  char name[90];
  double *histo_r, h_r;
  int *width;
  double mu, sigma, alpha, dens;
  double hinf, codinf;
  double h1;
  double tch, dh;
  int outres;
  int dim_space, nbox;
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  dim_space = p_frac->dim_space;
  mu = p_frac->mu, sigma = p_frac->sigma;
  alpha = p_frac->alpha, dens = p_frac->dens;
  hinf = p_frac->hinf, codinf = p_frac->codinf, h1 = p_frac->h1;
  lmax = p_frac->lmax, outres = p_frac->outres;
  nbox = p_frac->nbox;
#else
  dim_space =  DSPACE;
  mu = MU, sigma = SIGMA;
  alpha = ALPHA, dens = DENS;
  hinf = HINF, codinf = CODINF, h1 = H1;
  lmax = LMAX, outres = OUTRES;
  nbox = NBOX;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/

  block = (int)pow(2.,(double)outres);
  leff = dimensiona(lmax)/block;
  dimy = (dim_space==DIM1D) ? 1 : leff;
  
  /*      Initializating memory           */
  TrackNullAlloc( prob_exp=(double*)calloc(nbox,sizeof(double)) );
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  /* Preparing some multifractal parameters */
  prepara_multifractal();
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
  
  /*		Data reading		*/  
  IFVERBOSE WarningVV("Processing for %d output %dD signals",
		      nseries,dim_space);
  
  if(leff<2)       Error("Too reduced resolution!!!\n");

  TrackNullAlloc( signal=matrix2D(dimy,leff) );

  /* Composing output file names */
  genera_base( base, leff);
  
  /* Generating the multifractals */
  for(in=0;in<nseries;in++)    {
    
    genera_multifractal(block*leff,WAVBASE,DERWAVBASE,prob_exp,signal);
    
    sprintf(name,"%s-N%05d",base,in);
    /* #ifdef _PARSE_FRACTAL_PARAMETERS_
       IF(p_frac->flag_savefloat)
       #else
       IF(FLAG_SAVEFLOAT)
       #endif 
      write_data_infloat(leff,dimy,name,signal);
      else */
    write_data(leff,dimy,name,signal);
    IFVERBOSE
      if(dim_space == DIM1D) 
	write_serie(leff,dimy,name,signal[0]);
      else 
	write_image_foto(leff,leff,name,signal);
  }
  
  free_matrix2D(signal,dimy);

  /*	Freeing memory before finishing	*/
 free(prob_exp);

  
  return leff;
} // end of create_multifractal




/***************************************************************************/
int analize_multifractal( /*inputs*/char *name_in, int dimx, 
			  double *Moms, int nmoms, int *Dist, int ndists ) {
  /***************************************************************************/
  double *h, *Dh, *errDh;
  double shift;
  int N_h;
  double output=-1.;
  
  int dim_space, ndata, type;
  int dimy, nbox;
  int flag_fromdh, flag_holder;  
 
#ifdef _PARSE_FRACTAL_PARAMETERS_
  dim_space = p_frac->dim_space, ndata = p_frac->ndata;
  type = p_frac->type_mfana;
  dimy = ((DSCALE==DIM1D) ? 1 : dimx);
  nbox = p_frac->nbox;
  flag_fromdh = p_frac->flag_fromdh, flag_holder = p_frac->flag_holder;
#else
  dim_space = DSPACE, ndata = NDATA;
  type = TYPE_MFANA;
  dimy = ((p_frac->dim_scale==DIM1D) ? 1 : dimx);
  nbox = NBOX;
  flag_fromdh = FLAG_FROMDH, flag_holder = FLAG_HOLDER;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
  
  IF(flag_fromdh) {	
    sprintf( name, _DH_NAME_IN_FORMAT_,typ_name[type], gen_name, ndata );
    /* first pass to get the size */
    TrackError( N_h=Dh_read(name,NULL,NULL,NULL) );
    
  } else {
    switch(type) {
    case TYPWTMM:
      N_h = nmoms; break;
    case TYPMOM:
      N_h = nmoms - 1; break;
    default:      
    case TYPGRAD:
      N_h = nbox;
      break;
    case TYPGMWP:
      N_h = nbox; sc = 1./((double)dimx);
      if(dim_space == DIM1D) sc = S0 * wavelet1D_escala(dimx);
      else                   sc = S0 * wavelet2D_escala(dimx,dimx);
    }
  }
    
  /* Allocation for histograms and error bars */
  TrackNullAlloc( h=(double*)calloc(N_h,sizeof(double)) );
  TrackNullAlloc( Dh=(double*)calloc(N_h,sizeof(double)) );
  TrackNullAlloc( errDh=(double*)calloc(N_h,sizeof(double)) );  
  
  IF(flag_fromdh) {
    sprintf( name, _DH_NAME_IN_FORMAT_,typ_name[type], gen_name, ndata );
    Dh_read( name, h, Dh, errDh );
    
  } else {
    switch(type)    {
    case TYPWTMM:
      Dh_estima_wtmm( gen_name, dimx, Moms, Nmom, h, Dh, errDh );
      break;
    case TYPMOM:
      Dh_estima_moments( gen_name, ndata, dimx, dimy, 
			 Moms, nmoms, Dist, ndists, 
			 h, Dh, errDh );
      break;
    case TYPGRAD:
      Dh_estima_gh( gen_name, dimx, h, Dh, errDh );
      break;
    default:      
    case TYPGMWP:
      Dh_estima_gmwp( name, dimx, &shift_w, &shift_g, h_w, Dh_w, errDh_w );
      IF(!flag_holder) {
	shift = shiftg - shiftw;
	IFVERBOSE
	  WarningV("Correction exponent due to the lack of translational invariance: %f\n",
		   shift);
	for(ir=0;ir<Nh_w;ir++) h_w[ir] += shift;
      }
    }	
  }
    
  /*     Summarizing results        */
  IFVERBOSE
    WarningV(" Quality estimation for %s method"
	     "\n =================================",typ_name[type]);

  /* Save the Dh */   
  if(flag_fromdh) {
    sprintf(name,_DH_NAME_OUT_FORMAT_,typ_name[type],gen_name,ndata);
    Dh_write( name, N_h, h, Dh, errDh );
  }
  
  /* Compute the quality */
  output = Dh_registra( N_h, sc, h, Dh, errDh );

  /* Freeing memory before finishing */
  Free(h);
  Free(Dh); 
  Free(errDh);
    
  return output;
} // end of analiza_data



/****************************************************************************/
int old_analiza_series( int leff, int dimy, char *base) {
  /****************************************************************************/
  
  char nombre[90];
  double *h_g,*h_w,*Dh_g,*Dh_w,*errDh_g,*errDh_w;
  double sc;

  int Nh_g,Nh_w;

  int nbox;

  TrackNullAlloc( h_g=(double*)calloc(NBOX,sizeof(double)) ); 
  TrackNullAlloc( h_w=(double*)calloc(NBOX,sizeof(double)) );
  TrackNullAlloc( Dh_g=(double*)calloc(NBOX,sizeof(double)) );
  TrackNullAlloc( Dh_w=(double*)calloc(NBOX,sizeof(double)) );
  TrackNullAlloc( errDh_g=(double*)calloc(NBOX,sizeof(double)) );
  TrackNullAlloc( errDh_w=(double*)calloc(NBOX,sizeof(double)) );
  
  estima_Dh_file(leff,dimy,base,&Nh_g,h_g,Dh_g,errDh_g,&Nh_w,h_w,Dh_w,errDh_w);
  
  Warning("Estimation for derivative method\n================================\n");
  sprintf(nombre,"Dh_der_%s",base);
  registra_Dh_file(Nh_g,nombre,1./((double)leff),h_g,Dh_g,errDh_g);

  if(D_space==1) 
    sc=S0*wavelet1D_escala(leff,SC0_1D[ORDDER][WAV],WAV,ORDDER);
  else 
    sc=S0*wavelet2D_escala(leff,leff,SC0_2D[ORDDER][WAV],ORDDER,EXP_WL_2D[WAV]);

  Warning("Estimation for wavelet method\n=============================\n");
  sprintf(nombre,"Dh_sing_%s",base);
  registra_Dh_file(Nh_w,nombre,sc,h_w,Dh_w,errDh_w);

  /*	FREEING MEMORY BEFORE FINISHING */

  free(h_w);
  free(h_g);
  free(Dh_g);
  free(Dh_w); 
  free(errDh_g);
  free(errDh_w); 

  return OK;
}



/****************************************************************************/
int old_Dh_estima( int nseries, int leff, int dimy, char *base, 
	       int *Nh_g, double *h_g, double *Dh_g, double *errDh_g, 
	       int *Nh_w, double *h_w, double *Dh_w, double *errDh_w ) {
  /****************************************************************************/
  char name[MAXNAMELENGTH];
  double **signal;
  double **expon_g,**expon_w;
  double *histo_g,*histo_w;

  double m_g[2],m_w[2];
  double m_g0[2],m_w0[2];
  double shift,sc;

  int in,ir;

  TrackNullAlloc( signal=matrix2D(dimy,leff) );
  TrackNullAlloc( expon_g=matrix2D(dimy,leff) );
  TrackNullAlloc( expon_w=matrix2D(dimy,leff) );
  
  TrackNullAlloc( histo_g=(double*)calloc(NBOX,sizeof(double)) );
  TrackNullAlloc( histo_w=(double*)calloc(NBOX,sizeof(double)) );
  
  /*     First passage is to compute the maxima      */
  m_g[0] = m_w[0] = 1.e30;
  m_g[1] = m_w[1] = -1.e30;

  for(in=0;in<nseries;in++)    {
    sprintf(name,"%s-N%03d",base,in);
    read_data(leff,(D_space==1)?1:leff,name,signal);
    genera_expon(leff,signal,expon_g);
#ifdef _PARSE_FRACTAL_PARAMETERS_
    if(p_frac->dim_space == DIM1D) 
#else
      if(DSPACE == DIM1D) 
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
	extrema1D(leff,expon_g[0],&m_g0[0],NULL);
      else extrema(leff,leff,expon_g,&m_g0[0],NULL);
    m_g[0] = fMin(m_g[0],m_g0[0]);
    m_g[1] = fMax(m_g[1],m_g0[1]);
    
    m_w0[0] = 1.e30,    m_w0[1] = -1.e30;
      
    calcula_multifractal(leff,signal,expon_w,m_w0);
    
    m_w[0] = fMin(m_w[0],m_w0[0]);
    m_w[1] = fMax(m_w[1],m_w0[1]);    
  }
  
  /*      Second passage: the histograms are calculated    */  
  shift=0.;
  for(in=0;in<nseries;in++)    {
    sprintf(name,"%s-N%03d",base,in);
    read_data(leff,(D_space==1)?1:leff,name,signal);
    genera_expon(leff,signal,expon_g);
    
    m_w0[0]=1e30;
    m_w0[1]=-1e30;
    calcula_multifractal(leff,signal,expon_w,m_w0);
    IFVERBOSE Warning("\n");
    
    /* G-histograms are well referenced, so they are already accumulated 
     * H-histograms need to be re-referrenced if they are computed from   
     *      the measure                                                        */
    accumula_histogram(leff,m_g,expon_g,histo_g);
    accumula_histogram(leff,m_w,expon_w,histo_w);
    if(!(HOLDER||FLAG_INVTRANS))
      if(D_space==1) 
	shift += moda1D(leff,expon_g[0]) - moda1D(leff,expon_w[0]);    
      else 
	shift += moda(leff,leff,expon_g) - moda(leff,leff,expon_w);    
  }
  
  *Nh_g=calcula_Dh_histo(1./((double)leff),&m_g[0],histo_g,
			 h_g,Dh_g,errDh_g);
  
  if(D_space==1) 
    sc=S0*escala1D_wavelet(leff,SC0_1D[ORDDER][WAV],ORDDER,EXP_WL_1D[WAV]);
  else 
    sc=S0*escala2D_wavelet(leff,leff,SC0_2D[ORDDER][WAV],ORDDER,EXP_WL_2D[WAV]);

   *Nh_w=calcula_Dh_histo(sc,&m_w[0],histo_w,h_w,Dh_w,errDh_w);
   
   /*     If necessary, correcting by shift measure wavelet projetcions    */   
   if(!(HOLDER||FLAG_INVTRANS))     {
     shift/=(double)NSERIES;
     WarningV("Correction exponent due to lack of translational invariance: %f\n",shift);     
     for(ir=0;ir<*Nh_w;ir++) h_w[ir]+=shift;
   }
   
   /*	FREEING MEMORY BEFORE FINISHING */
  free(histo_g);
  free(histo_w);

  free_matrix2D(expon_g,dimy);
  free_matrix2D(expon_w,dimy);
  free_matrix2D(signal,dimy);

  return OK;
} // end of estima_Dh_file

