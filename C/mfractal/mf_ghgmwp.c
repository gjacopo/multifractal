/* ===================================
** mf_ghgmwp.c
** started on Wed Jan 31 18:30:31 2007 
** ===================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mf_ghgmwp.h>


#ifdef _PARSE_FRACTAL_PARAMETERS_
#include <mf_parse.h>
extern ParFRAC *p_frac;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */




/***************************************************************************/
int Dh_estima_gh( /*inputs*/char *name_in, int dimx, 
		 /*outputs*/double *shift_g, 
		 double *h_g, double *Dh_g, double *errDh_g ) {
  /***************************************************************************/
  char name[MAXNAMELENGTH];
  double **signal;
  double **expon_g;
  double *histo_g;
  double m_g[2],m_g0[2];
  double sc;
  int in,ir;
  double *histo_g_r, *h_g_r;
  int *width_g;
  
  int dimy, ndata, nbox;
  
#ifdef _PARSE_FILTER_PARAMETERS_
  ndata = p_frac->ndata,  dim_space = p_frac->dim_space;
  nbox = p_frac->nbox;
#else
  ndata = NDATA,  dim_space = DSPACE;
  nbox = NBOX;
#endif /*!_PARSE_FILTER_PARAMETERS_*/
  dimy = ((dim_space==DIM1D) ? 1 : dimx);
  
  TrackNullAlloc( signal=matrix2D(dimy,dimx) );
  TrackNullAlloc( expon_g=matrix2D(dimy,dimx) );
  TrackNullAlloc( histo_g=(double*)calloc(nbox,sizeof(double)) );
  
  /* First passage is to compute the maxima      */
  m_g[0] = 1.e30,    m_g[1] = -1.e30;
  for( in=0; in<ndata; in++ ) {
    if(ndata > 1) sprintf(name,"%s-N%05d",name_in,in);
    else          strcpy(name,name_in);
    read_data( dimx, dimy, name, signal );
    
    expon_genera( dimx, signal, expon_g );
    
    if(dim_space == DIM1D) extrema1D( dimx, expon_g[0], m_g0, NULL );
    else                   extrema( dimx, dimy, expon_g, m_g0, NULL );
    
    m_g[0] = fMin(m_g[0],m_g0[0]),  m_g[1] = fMax(m_g[1],m_g0[1]);
  }
  
  if(dim_space == DIM1D)  sc = 1. / (double)dimx;
  else                    sc = escala2D_lineal(dimx, dimy);

  /* Second passage: the histograms are calculated    */
  for( in=0; in<ndata; in++ ) {
    if(ndata > 1) sprintf(name,_DATA_NAME_IN_FORMAT_,name_in,in);
    else          strcpy(name,name_in);
    read_data( dimx, dimy, name, signal );
    
    expon_genera( dimx, signal, expon_g );
    
    /* If this method is explicitly invoked, we take profit of that to
     * save time in the estimation of the shift for the GMWP  */
    histogram_accumula( dimx, dimy, m_g, expon_g, histo_g );
  }
  
  TrackNullAlloc( histo_g_r=(double*)calloc(nbox,sizeof(double)) );
  TrackNullAlloc( h_g_r=(double*)calloc(nbox,sizeof(double)) );
  TrackNullAlloc( width_g=(int*)calloc(nbox,sizeof(int)) );
  ifill1D(nbox,1,width_g,NULL); 
 
  Nh_g = Dh_filter( histo_g, m_g, histo_g_r, h_g_r );
  Dh_errorbar_weight( 1./((double)dimx), Nh_g, histo_g_r, h_g_r, width_g, 
		      h_g, Dh_g, errDh_g );
  
  (*shift_g) = moda1D_histo( nbox, m_g, histo_g );
  
  /* Freeing memory before finishing */
  free(histo_g_r);
  free(h_g_r);
  free(width_g); 
  free(histo_g);
  free_matrix2D(expon_g,dimy);
  free_matrix2D(signal,dimy);
  
  return OK;
} // end of Dh_estima_g


/***************************************************************************/
int Dh_estima_gmwp( /*inputs*/char *name_in, int dimx, 
		    /*outputs*/double *shift_w,
		    double **h_w, double **Dh_w, double **errDh_w ){
  /***************************************************************************/
  char name[MAXNAMELENGTH];
  double **signal;
  double **expon_w;
  double *histo_w;
  double m_w[2],m_w0[2];
  double sc;
  int Nh_w;
  int dimy;
  int in,ir;
  
  int dimy, ndata, nbox;
  
#ifdef _PARSE_FILTER_PARAMETERS_
  ndata = p_frac->ndata;
  dimy = ((p_frac->dim_scale==DIM1D) ? 1 : dimx);
  nbox = p_frac->nbox;
  sc0 = p_frac->sc0;
#else
  ndata = NDATA;
  dimy = ((DSCALE==DIM1D) ? 1 : dimx);
  nbox = NBOX;
  sc0 = MINSCALE;
#endif /*!_PARSE_FILTER_PARAMETERS_*/

  TrackNullAlloc( signal=matrix2D(dimy,dimx) );
  TrackNullAlloc( expon_w=matrix2D(dimy,dimx) );
  TrackNullAlloc( histo_w=(double*)calloc(nbox,sizeof(double)) );
  
  // Implicit call to GH to obtain shift_g
  /* estima_Dh_g(name_in,&Nh_g,&h_g,&Dh_g,&errDh_g);
     free(h_g);
     free(Dh_g);
     free(errDh_g);
  */
  
  /*     First passage is to compute the maxima      */
  m_w[0]=1e30, m_w[1]=-1e30;
  for( in=0; in<ndata; in++ ) {
    if(ndata > 1) sprintf(name,_DATA_NAME_IN_FORMAT_,name_in,in);
    else          strcpy(name,name_in);
    read_data( dimx, dimy, name, signal );
    
    m_w0[0]=1.e30,      m_w0[1]=-1.e30;
    calcula_multifractal( dimx, signal, expon_w, m_w0 );
    
    m_w[0] = fMin(m_w[0], m_w0[0]),   m_w[1] = fMax(m_w[1], m_w0[1]);
  }
  
  /* Second passage: the histograms are calculated and the shift, corrected */
  if(dim_space==DIM1D) sc = wavelet1D_escala(dimx,sc0);
  else                 sc = wavelet2D_escala(dimx,dimy,sc0);
  
  for( in=0; in<ndata; in++ ) {
    if(ndata > 1) sprintf(name,"%s-N%05d",name_in,in);
    else          strcpy(name,name_in);
    read_data( dimx, dimy, name, signal );
    
    m_w0[0]=1.e30,      m_w0[1]=-1.e30;
    calcula_multifractal( dimx, signal, expon_w, m_w0 );
   
    histogram_accumula( dimx, dimy, m_w, expon_w, histo_w );
  }
  
  TrackNullAlloc( histo_w_r=(double*)calloc(nbox,sizeof(double)) );
  TrackNullAlloc( h_w_r=(double*)calloc(nbox,sizeof(double)) );
  TrackNullAlloc( width_w=(int*)calloc(nbox,sizeof(int)) );
  ifill1D(nbox,1,width_w,NULL); 
 
  Nh_w = Dh_filter( histo_w, m_w, histo_w_r, h_w_r );
  Dh_errorbar_weight( sc, Nh_w, histo_w_r, h_w_r, width_w, 
		      h_w, Dh_w, errDh_w );

  (*shift_w) = moda1D_histo( nbox, m_w, histo_w );

  /* Freeing memory before finishing */
  free(histo_w_r);
  free(h_w_r);
  free(width_w); 
  free(histo_w);
  free_matrix2D(expon_w,dimy);
  free_matrix2D(signal,dimy);
  
  return OK;
} // end of Dh_estima_gmwp


/***************************************************************************/
int Dh_filter_weight( /*input*/ double *histo, double *mh, 
		      /*outputs*/ double *histo_r, double *h_r, int *width ) {
  /***************************************************************************
  .* Called by Dh_estima_gh
  ***************************************************************************/
  int ip,ir=0;
  double cump=0.;
  
  int nbox=
#ifdef _PARSE_FRACTAL_PARAMETERS_
    p_frac->nbox;
#else
  NBOX;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
  
  /* Filtering histogram to avoid low-probability distorsions    */
  for( ip=0; ip<nbox; ip++ )    {
    cump += histo[ip];
    h_r[ir] += (mh[0]+(0.5+(double)ip)/((double)nbox)*(mh[1]-mh[0]))
      * histo[ip];
    if(cump > MINEV)      {
      width[ir] = ip - ip0;
      histo_r[ir] = cump / ((double)width[ir]);
      h_r[ir] /= cump;
      cump = 0.;
      ip0 = ip;
      ir++;
    }
  }
  
  if(cump > 0.)    {
    width[ir] = nbox - 1 - ip0;
    histo_r[ir] = cump / ((double)width[ir]);
    h_r[ir] /= cump;
    ir++;	  
  }
  
  /* Note: necessary ir<=nbox */
  return ir; 
} // end of Dh_filter_weight


/***************************************************************************/
int Dh_filter( /*input*/ double *histo, double *mh, 
	       /*outputs*/ double *histo_r, double *h_r ) {
  /***************************************************************************
  .* Called by Dh_estima_gh
   *           Dh_estima_gmwp
  ***************************************************************************/
  int ip,ir=0;
  double cump=0.;
  
  int nbox=
#ifdef _PARSE_FRACTAL_PARAMETERS_
    p_frac->nbox;
#else
  NBOX;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
  
  /* Filtering histogram to avoid low-probability distorsions    */
  for(ip=0;ip<nbox;ip++)    {
    cump += histo[ip];
    h_r[ir] += (mh[0]+(0.5+(double)ip)/((double)nbox)*(mh[1]-mh[0]))
      * histo[ip];
    if(cump > MINEV)      {
      histo_r[ir] /= cump;
      h_r[ir] /= cump;
      cump = 0.;
      ir++;
    }
  }
  
  if(cump > 0.)    {
    histo_r[ir] /= cump;
    h_r[ir] /= cump;
    ir++;	  
  }
  
  return ir;
} // end of Dh_filter


/***************************************************************************/
int Dh_errorbar_weight( double sc0, int Nh,
			/*inputs*/ double *histo_r, double *h_r, int *width, 
			/*outputs*/ double *h, double *Dh, double *errDh ) {
  /***************************************************************************
  .* Called by Dh_estima_gh
   *           Dh_estima_gmwp
  ***************************************************************************/
  double dim_space;
  double prob, maxprob, dp, Nev;
  int ip, ir;

#ifdef _PARSE_FRACTAL_PARAMETERS_
  dim_space = (double) p_frac->dim_space;
  //sc0 =  p_frac->minscale;
#else
  dim_space = (double)DSPACE;
  // sc0 = MINSCALE;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */

  /*                Producing the error bars                         */ 

  //  First, we compute the total number of events
  Nev = 0.;
  for(ip=0;ip<Nh;ip++) Nev += histo_r[ip] * (double)width[ip];
  
  /* The confidence range is taken as +- SIGCONFIDENCE sigmas in the  
   * distribution of probability boxes, which is a renormalized binomial. 
   * In this routine we directly propagate the confidence range to the D(h)
   * calculation               */
  for(ip=0;ip<Nh;ip++)   {
    h[ip] = h_r[ip]; // associated singularity value
    prob = histo_r[ip] * (double)width[ip] / Nev; // probability of that interval
    if(prob > 1.e-30) dp = SIGCONFIDENCE * sqrt((1.-prob) / (prob*Nev)); 
    // SIGCONFIDENCE sigmas in the prob. distribution
    else dp = 0.;
    
    /* The error bar is constructed by propagation; as the propagated interval
     * is asymetric, we take the bar as the maximum of the two distances to
     * the  central value. This is given by the lower bound in the logarithm
     */
    if(dp >= 1.) errDh[ip] = dim_space;
    else errDh[ip] = log(1.-dp) / log(sc0);
  } 
  
  /* Finding and normalizing by the mode  */
  maxprob = histo_r[0];
  for( ip=1; ip<Nh; ip++ ) maxprob = fMax(maxprob,histo_r[ip]);
  if(maxprob > 1.e-30) for( ip=0; ip<Nh; ip++ ) histo_r[ip] /= maxprob;
  
  /* Evaluating experimental D(h) */
  /*old: ir = 0; */
  for( ip=0; ip<Nh; ip++ )    {
    /*old: h[ip] = mh[0] + (0.5*(double)width[ip]+(double)ir)
       / ((double)nbox) * (mh[1]-mh[0]); */
    if(histo_r[ip] > 1.e-30) 
      Dh[ip] = dim_space - log(histo_r[ip]) / log(sc0);
    else 
      Dh[ip] = 0.; // -dim_space;
    /*old: ir += width[ip]; */
  }
  
  return OK;
} // end of Dh_errorbar_weight



/***************************************************************************/
int Dh_histo_weight( double sc0, double *mh, double *histo ) {
  /***************************************************************************
  .* Called by Dh_estima_gh
   *           Dh_estima_gmwp
  ***************************************************************************/
  double *h,  *Dh, *errDh;
  int Nh;

  int nbox=
#ifdef _PARSE_FRACTAL_PARAMETERS_
    p_frac->nbox;
#else
  NBOX;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */

  /*          Initialization                    */
  TrackNullAlloc( histo_r=(double*)calloc(nbox,sizeof(double)) );
  TrackNullAlloc( h_r=(double*)calloc(nbox,sizeof(double)) );
  TrackNullAlloc( width=(int*)calloc(nbox,sizeof(int)) );

  Nh = Dh_filter_weight( histo, histo_r, h_r, width );
  
  /*  The number of exponents is now known; we allocate memory accordingly */
  TrackNullAlloc( h=(double*)calloc(Nh,sizeof(double)) );
  TrackNullAlloc( Dh=(double*)calloc(Nh,sizeof(double)) );
  TrackNullAlloc( errDh=(double*)calloc(Nh,sizeof(double)) );

  Dh_errorbar_weight(sc0, Nh, histo_r, h_r, width, h, Dh, errDh );
  
  free(histo_r);
  free(h_r);
  free(width);
  
  return Nh;
} // end of Dh_histo_weight


/***************************************************************************/
int Dh_histo( double sc0, double *mh, double *histo ) {
  /***************************************************************************
  .* Called by Dh_estima_gh
   *           Dh_estima_gmwp
  ***************************************************************************/
  double *h,  *Dh, *errDh;
  int Nh;

  int nbox=
#ifdef _PARSE_FRACTAL_PARAMETERS_
    p_frac->nbox;
#else
  NBOX;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */

  /*          Initialization                    */
  TrackNullAlloc( histo_r=(double*)calloc(nbox,sizeof(double)) );
  TrackNullAlloc( h_r=(double*)calloc(nbox,sizeof(double)) );

  Nh = Dh_filter( histo, histo_r, h_r );
  
  /*  The number of exponents is now known; we allocate memory accordingly */
  TrackNullAlloc( h=(double*)calloc(Nh,sizeof(double)) );
  TrackNullAlloc( Dh=(double*)calloc(Nh,sizeof(double)) );
  TrackNullAlloc( errDh=(double*)calloc(Nh,sizeof(double)) );

  /* create a mute variable */
  TrackNullAlloc( width=(int*)calloc(nbox,sizeof(int)) );
  ifill1D(nbox,1,width,NULL); 
  
  Dh_errorbar_weight(sc0, Nh, histo_r, h_r, h, Dh, errDh );
  
  free(histo_r);
  free(h_r);
  free(width);
  
  return Nh;
} // end of Dh_histo



/***************************************************************************/
void histogram_accumula( int dimx, int dimy, double *mh, double **expon, 
			 double *histo) {
  /***************************************************************************
   * Called by Dh_estima_gh 
   *           Dh_estima_gmwp
   ***************************************************************************/
  int ix,iy,ip;
  
  int nbox=
#ifdef _PARSE_FRACTAL_PARAMETERS_
    p_frac->nbox;
#else
  NBOX;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
  
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)      {
      ip = (int)((expon[iy][ix]-mh[0]) / (mh[1]-mh[0])*nbox);
      if(ip < 0) ip = 0;
      if(ip >= nbox) ip = nbox-1;
      histo[ip] += 1.;		
    }

} // end of histogram_accumula

