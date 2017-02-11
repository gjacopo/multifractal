/* ===================================
** mf_moments.c
** started on Wed Jan 31 18:05:57 2007 
** ===================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <utils.h>
#include <utl_alloc.h>
#include <utl_stats.h>
#include <utl_operator.h>

#include <filter.h>
#include <flt_deriva1d.h>

#include <mf_moments.h>


/***************************************************************************/
int Dh_estima_moments( /*inputs*/char *name_in, int ndata,
		       int dimx, int dimy,
		       double *Moms, int nmoms, int *Dist, int ndists, 
		       /*outputs*/double *h_m, double *Dh_m, double *errDh_m ) {
  /***************************************************************************/
  char name[MAXNAMELENGTH];
  double **signal;
  double **moments, *taup;
  double sc;
  int in;
  
  TrackNullAlloc( signal=matrix2D(dimy,dimx) );
  TrackNullAlloc( moments=matrix2D(nmoms,ndists) );
  TrackNullAlloc( taup=(double*)calloc(nmoms,sizeof(double)) );
  
  IFVERBOSE {
    Warning("Compute Moments");
    if(ndata > 1) 
      WarningV("The statistics of the %d signals are accumulated.",ndata);
  }

  /* First step is to compute the moments      */  
  for( in=0; in<ndata; in++ ) {
    if(ndata > 1) sprintf(name,_DATA_NAME_IN_FORMAT_,name_in,in);
    else          strcpy(name,name_in);
    read_data( dimx, dimy, name, signal );
    
    moments_accumula( dimx, dimy, signal, Mom, nmoms, Dist, ndists,
		      moments );
  }
  
  /* It is necessary to normalize by the number of series */
  op_scale( ndists, nmoms, 1./(double)ndata, moments, NULL );
  
  /* Second step: taup is evaluated    */
  calcula_taup( moments, nmoms, Dist, ndists, taup );
  
  /* Third step: from taup, Dh is evaluated    */
  legendre_transform( Moms, nmoms-1, taup, h_m, Dh_m, errDh_m );
    
  /* Freeing memory before finishing */
  free(taup);
  free_matrix2D(moments,nmoms);
  free_matrix2D(signal,dimy);
  
    return OK;
} // end of Dh_estima_moments


/***************************************************************************/
int moments_accumula( /*inputs*/int dimx, int dimy, double **signal,
		      double *Mom, int nmoms, int *Dist, int ndists,
		      /*output*/double **moments ) {
  /***************************************************************************
  .* Called by Dh_estima_moments
  ***************************************************************************/
  double *partial_mom;
  double *grad;
  double eps;
  int xeff, dmax;
  int ix,j,iy,ip;
  int it;
  
  /*     Initialization    */
  TrackNullAlloc( grad=(double*)calloc(dimx,sizeof(double)) );
  TrackNullAlloc( partial_mom=(double*)calloc(nmoms,sizeof(double)) );
  
  for(it=0;it<ndists;it++)    {
    
    /* Defining the number of intervals of size t that can be included    */
    dmax = dimx - Dist[it] - 1;
    
    /* Cleaning the partial moments in which the results at this distance
     * are accumulated. We proceed in this way to avoid losing contributions
     * which are relevant for the single series under process but what could 
     * be lost when compared to all the accumulated course. */
    for(ip=0;ip<nmoms;ip++) partial_mom[ip] = 0.;
    
    /* We will now run all the (horizontal) ranges of point of lenght t, and for 
     * them we evaluate their contribution to the different moments at this size. */
    for(iy=0;iy<dimy;iy++)      {
      
      /*    Initializating the weight list   */
      copy1D(dimx,signal[iy],grad,NULL);
      gradient1D(dimx,grad);
      for(ix=0;ix<dimx;ix++) grad[ix] = fabs(grad[ix]);
      
      /*    Computing the contributions of all the ranges of size dist[it]  */
      /* To increase computational speed, to compute the contribution of a new 
       * interval we update it with respect to the previous one; so, we just need 
       * to sum the contribution of the new point entering and to discount that
       * of the first point in the previous interval */
      
      /*      Contribution for ix=0    */
      eps = 0.;
      for( j=0; j<Dist[it]; j++ ) eps += grad[j];
      eps /= (double)Dist[it];
      for(ip=0;ip<nmoms;ip++)	
	if(eps > 1.e-10) partial_mom[ip] += pow(eps,Moms[ip]);
      
      for (ix=0;ix<dmax-1;ix++)	{
	eps += (grad[ix+Dist[it]]-grad[ix]) / ((double)Dist[it]);
	for(ip=0;ip<nmoms;ip++)	  
	  if(eps > 1.e-10) partial_mom[ip] += pow(eps,Moms[ip]);
      }
    }
    
    /*   Transferring the partial moments to the accumulated moments  */
    for(ip=0;ip<nmoms;ip++) 
      moments[ip][it] += partial_mom[ip] / ((double)dimy*dmax);
  }
  
  /*         Releasing memory before ending   */
  free(grad);
  free(partial_mom);
  
  return OK;
} // end of moments_accumula


/***************************************************************************/
void moments_taup( /*inputs*/double **moments, int nmoms, int *Dist, int ndists,
		   /*output*/double *taup ) {
  /***************************************************************************
  .* Called by Dh_estima_moments
  ***************************************************************************/
  double *lgd,*rows;
  double a,b,corr;
  int ip,it;
  
  TrackNullAlloc( lgd=(double*)calloc(ndists,sizeof(double)) );
  TrackNullAlloc( rows=(double*)calloc(ndists,sizeof(double)) );

  for( ip=0; ip<nmoms; ip++ ) {
    for ( it=0; it<ndists; it++ )	{
      lgd[it] = log(Dist[it]);
      rows[it] = log(moments[ip][it]);
    }
    fit(lgd,rows,ndists,&a,&b,&corr);
    taup[ip] = a;
  }

  free(lgd);
  free(rows);

  return OK;
} // end of moments_taup


/***************************************************************************/
int legendre_transform( double *moms, int N_m, 
			double *taup, double **h_m, double **Dh_m,
			double **errDh_m ) {
  /***************************************************************************
  .* Called by Dh_estima_moments
  ***************************************************************************/
  int ip;
  
  /* Obtaining the spectrum and its singularities   */
  h_m[0][0] = (taup[1]-taup[0]) / (moms[1]-moms[0]);
  Dh_m[0][0] = moms[0]*h_m[0][0] + 1 - taup[0];
  errDh_m[0][0] = 0.2; // conventionally
  
  for( ip=1; ip<N_m; ip++ )    {
    h_m[0][N_m-1-ip] = (taup[ip+1]-taup[ip-1]) / (moms[ip+1]-moms[ip-1]);
    Dh_m[0][N_m-1-ip] = moms[ip]*h_m[0][N_m-1-ip]+1-taup[ip];
    errDh_m[0][N_m-1-ip] = 0.2; // conventionally
  }
  
  return OK;
} // end of legendre_transform
