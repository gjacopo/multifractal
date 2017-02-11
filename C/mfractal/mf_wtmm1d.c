#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_alloc.h>
#include <utl_stats.h>
#include <utl_operator.h>

#include <filter.h>
#include <flt_deriva1d.h>
#include <flt_wavelet1d.h>

#include <mfractal.h>
#include <mf_wtmm1d.h>

extern int Verbose;

#ifdef _PARSE_FILTER_PARAMETERS_
extern ParWAV *p_wav;
#endif /*!_PARSE_FILTER_PARAMETERS_*/

#ifdef _PARSE_FRACTAL_PARAMETERS_
#include <mf_parse.h>
extern ParWTMM *p_wtmm;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/


/***************************************************************************/
int Dh_estima_wtmm( /*inputs*/
		   char *name_in, int dimx, double *Moms, int nmoms, 
		    /*outputs*/
		    double *h_wtmm, double *Dh_wtmm, double *errDh_wtmm ) {
  /***************************************************************************/
  char name[MAXNAMELENGTH];
  double **ExtWTlis, **Z;
  double *scArray;
  double *taup;
  int in,ix,iy,ip;
  double sc;
  int *n_ext;
  int ind_q, ind_sc, nsc=0;
  double **sTq=NULL, **sTqLogT=NULL;
  double **signal;

  int dimy, ndata;
  double minscale, maxscale; 
  double scstep;
  
#ifdef _PARSE_FILTER_PARAMETERS_
  minscale = p_wav->minscale,  scstep = p_wav->scstep;
  maxscale = p_wav->scratio*((double)dimx);
#else
  minscale = MINSCALE,  scstep = SCSTEP;
  maxscale = SCRATIO*((double)dimx);
  dimy = ((DSCALE==DIM1D) ? 1 : dimx);
#endif /*!_PARSE_FILTER_PARAMETERS_*/
#ifdef _PARSE_FRACTAL_PARAMETERS_
  ndata = p_frac->ndata;
  dimy = ((p_frac->dim_scale==DIM1D) ? 1 : dimx);
#else
  ndata = NDATA;
  dimy = ((p_frac->dim_scale==DIM1D) ? 1 : dimx);
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
  
  /*      Initialization of the dimensions    */
  for( sc=minscale; sc<=maxscale; sc*=scstep ) nsc++;
  
  /*     Memory initialization     */
  TrackNullAlloc( signal=matrix2D(dimy,dimx+1) );
  TrackNullAlloc( ExtWTlis=matrix2D(nsc,dimx+1) );
  TrackNullAlloc( Z=matrix2D(nsc,nmoms+1) );
  TrackNullAlloc( scArray=(double*)calloc(nsc,sizeof(double)) );
  TrackNullAlloc( n_ext=(int*)calloc(nsc+1,sizeof(int)) );
  TrackNullAlloc( taup=(double*)calloc(nmoms,sizeof(double)) );
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  if(p_wtmm->flag_methodwtmm == CANONICAL)   
#else
    if(FLAG_METHODWTMM == CANONICAL)   
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
      {
	TrackNullAlloc( sTq=matrix2D(nsc,nmoms) );
	TrackNullAlloc( sTqLogT=matrix2D(nsc,nmoms) );
      }
  
  /*     Value initialization    */
  for (ind_sc=1, scArray[0]=minscale; ind_sc<nsc; ind_sc++) 
    scArray[ind_sc] = scArray[ind_sc-1] * scstep;
  
  IFVERBOSE {
    Warning("Compute Wavelet Transform, Wavelet Extrema and Partition Function"
	    " over all scales");
    if(ndata > 1) 
      WarningV("The statistics of the %d signals are accumulated.",ndata);
  }
  
  for( in=0; in<ndata; in++ ) {
 if(ndata > 1) sprintf(name,_DATA_NAME_IN_FORMAT_,name_in,in);
   else          strcpy(name,name_in);
   read_data( dimx, dimy, name, signal ); 
    
    /* Generate and track maxima lines for extrema extraction */
    wtmm1d_compute( dimx, dimy, signal, nsc, Moms, nmoms, ExtWTlis, Z, n_ext, 
		    scArray, sTq, sTqLogT );
  }
  
  /* free temporay used memory */
  free(n_ext);
  free_matrix2D(ExtWTlis,nsc);
  
  IFVERBOSE Warning("Approximate Legendre Transform for exponents estimation");
  
  /* Compute the exponents tau(q) and the spectrum (h,D(h)) from the values of
   * the WT over the detected extrema **/
  wtmm1d_estima( dimx, nsc, Moms, nmoms, Z, n_ext, scArray,
		 h_wtmm, Dh_wtmm, errDh_wtmm, taup, sTq, sTqLogT );
  
  /*	Freeing memory before finishing */
  free(scArray);
  free(taup);
#ifdef _PARSE_FRACTAL_PARAMETERS_
  if(p_wtmm->flag_methodwtmm == CANONICAL)
#else
    if(FLAG_METHODWTMM == CANONICAL)
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
      {
	free_matrix2D(sTqLogT,nsc);
	free_matrix2D(sTq,nsc);
      }  
  free_matrix2D(Z,nsc);
  free_matrix2D(signal,dimy);
  
  return OK;
} // end of Dh_estima_wtmm


/***************************************************************************/
int wtmm1d_compute( int dimx, int dimy, double **signal, 
		    int nsc, double *Moms, int nmoms, 
		    double **ExtWTlis, double **Z, int *n_ext, 
		    double *scArray, double **sTq, double**sTqLogT ) {
  /***************************************************************************
  .* Called by Dh_estima_wtmm
  ***************************************************************************/
  int iy;
  
  /* Generate and track maxima lines for extrema extraction */
  for( iy=0; iy<dimy; iy++ )      {
    /* clean here the WTMM first : ExtWTlis = 0 everywhere, even if we
       rewrite on it  */
    fill0( dimx+1, nsc, ExtWTlis, NULL );
    
    /* Note: even the input signal is 2D, the analysis performed is still
     * a 1D analysis */
    wtmm1D_pf( signal[iy], dimx, 
#ifdef _PARSE_FRACTAL_PARAMETERS_
	       (double)p_wtmm->wav_wtmm, p_wtmm->ord_der_wtmm,
#else
	       (double)WAVWTMM, ORDDERWTMM,
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
	       nmoms, Moms, nsc, scArray, ExtWTlis, Z, n_ext );
    
    if(
#ifdef _PARSE_FRACTAL_PARAMETERS_
       (p_wtmm->flag_methodwtmm==CANONICAL) &&
#else
       (FLAG_METHODWTMM==CANONICAL) && 
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
       (sTqLogT!=NULL) && 
       (sTq!=NULL) ) 
      canonmean1D( ExtWTlis, dimx, n_ext, nmoms, Moms, nsc, scArray, sTq, sTqLogT );
  }
  
  return OK;
} // end of wtmm1d_compute


/***************************************************************************/
int wtmm1d_estima( int dimx, int nsc, double *Moms, int nmoms, 
		   double **Z, int *n_ext, double *scArray,
		   double *h_wtmm, double *Dh_wtmm, double *errDh_wtmm,
		   double *taup, double **sTq, double**sTqLogT ) {
  /***************************************************************************
  .* Called by Dh_estima_wtmm
  ***************************************************************************/
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  if(p_wtmm->flag_methodwtmm == DIRECT) 
#else
    if(FLAG_METHODWTMM == DIRECT)   
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
      /* Approximate Legendre Transform for estimating multifractal exponents
       * with DIRECT method */
      directspec1D( Z, dimx, nmoms, Moms, nsc, scArray, taup, h_wtmm, Dh_wtmm );
  
    else if((sTqLogT!=NULL) && (sTq!=NULL)) 
      /* Ibid with CANONICAL method */
      canonspec1D( sTq, sTqLogT, dimx, nmoms, Moms, nsc, scArray, taup,
		   h_wtmm, Dh_wtmm );
  
  /*       Padding proxy errorbars and shifting spectrum       */
  for(ip=0;ip<nmoms;ip++)     errDh_wtmm[ip]=0.2; // conventionally

  return OK;
} // end of wtmm1d_estima



/***************************************************************************/
int wtmm1D_pf( double *signal, int dimx, double wav, int order,
	       int nq, double *qArray, int nsc, double *scArray, 
	       double **ExtWTlis, double **Z, int *n_ext ) {
  /***************************************************************************
   * Find the extrema of WT for each scale of analysis, by performing the  
   * following steps:
   * 1) Wavelet convolution of the signal for increasing wavelet scale.
   * 2) Locate the local maxima of the absolute value of wavelet
   *    coefficient as a function of time for each wavelet scale.
   * 3) Check whether a local maximum at a given wavelet scale is located
   *    close to a maximum at a smaller scale - if yes connect both maxima,
   *    otherwise cancel it. Generate maxima lines.
   * 4) Check that the number of maxima at larger scales is less or equal
   *    to that at a smaller scale. 
   * 5) Track maxima lines for increasing wavelet scale by choosing at each
   *    scale the supremum between all previous values at smaller scales.
   *
   * Parse:
   *    - signal : original signal,
   *    - order : variable used for the choice of the wavelet,
   *    - maxfilt : maximal size of the wavelet filter (pre-computed),
   *    - dimx : size of the original signal,
   *    - nq, qArray : resp. no of moments and list of moments,
   *    - nsc, scArray : resp. no of sceles and list of scales.
   * Outputs:
   *    - Z : 2d tabular (scale by moment) storing the factors of the 
   *      partition function,
   *    - ExtWTlis : 2d tabular (scale by list of extrema) storing the WTMM,
   *      i.e. the values of the WT over selected extrema,
   *    - n_ext : tabular storing the number of extrema at each scale,
   * Returns the number icolor of data in color representation, if computed.
   *
   * Called by estima1D_wtmm
   ***************************************************************************/
  int i, j;
  double *maxsig, *wavesig, *pt, sc;
  double *wtrans;
  int *wtrans_ind, *maxsig_ind, *ptl;
  int ind_sc, ind_q, n_max;
  int jj=0; /* used to test wether the number of maxima at larger scales is  
	     * less or equal to that at a smaller scale.*/
  
  /* Allocations of memory
   * Note: making the allocations here (once for all, instead of making
   * it for each scale) enables to reduce time expense and allows
   * comparisons between scales */
  /* Arrays of the wavelet transform */
  TrackNullAlloc( wtrans=(double*)calloc(dimx+1,sizeof(double)) );
  TrackNullAlloc( wtrans_ind=(int*)calloc(dimx+1,sizeof(int)) );
  /* Arrays to manipulate the maxima through the scales */
  TrackNullAlloc( maxsig=(double*)calloc(dimx+1,sizeof(double)) );
  TrackNullAlloc( maxsig_ind=(int*)calloc(dimx+1,sizeof(int)) );
  
  for( ind_sc=0; ind_sc<nsc; ind_sc++ )   {
    sc=scArray[ind_sc]; /* current scale of analysis */
    
    /* Wavelet transform of the signal for increasing wavelet scale. */
    wt1D_transform( signal, dimx, wav, order, sc,  wtrans );
    
    /** Find the local maxima of the wavelet coefficient for each scale **/
    n_max = wtmm1D_extrema_find( wtrans, wtrans_ind, dimx, sc); 
    /* n_max : number of extrema at current scale */
    
    /** Tracking the maxima lines: test for supremum 
     ** Generate maxima lines and compute the partition function. */
    if( ind_sc > 0)    /* i.e.: if sc > min_sc */
      n_ext[ind_sc] =
	pf1D_extrema_track( wtrans, wtrans_ind, maxsig, maxsig_ind, 
			  jj, n_max, nq, qArray, nsc, scArray, 
			  ExtWTlis[ind_sc], Z[ind_sc] );
    
    /* Update for the next scale */
    jj = n_max;

    pt = maxsig;
    maxsig = wtrans; /* store in maxsig the WT of the current scale in order
		      * to compare it with the next scale */
    wtrans = pt; /* change the variable pointed by wtrans, so that
		  * further modification of wtrans won't modify maxsig */
    ptl = maxsig_ind;
    maxsig_ind = wtrans_ind; /* ibid with the indexes */ 
    wtrans_ind = ptl;  /* ibid with the indexes */
    
  } /* end loop over the scales: "for (ind_sc=0;..." */

  /* Free memories */
  if(wtrans) free(wtrans);
  if(wtrans_ind) free(wtrans_ind);
  if(maxsig) free(maxsig);
  if(maxsig_ind) free(maxsig_ind);
  
  return OK;
} // end of wtmm1D_pf


/***************************************************************************/
int wtmm1D_extrema_find( double *wtrans, int *wtrans_ind, int dimx, 
			double sc ) {
  /***************************************************************************
   * Find the local maxima of the wavelet coefficient for each scale.
   *
   * Parse:
   *    - sc : current scale of analysis,
   *    - wtrans : 1d array storing  the values of the WT at scale sc,
   *    - dimx : common size of both  wtrans and wtrans_ind.
   * Outputs:
   *    - wtrans : it is modified to finally store the values of the extrema
   *      only,
   *    - wtrans_ind : 1d array storing the indexes of the selected 
   *      extrema.
   * Returns the number n_max of selected extrema at scale sc.
   *
   * Called by : wtmm1D_pf  
   ***************************************************************************/
  int n_max=0;
  double temp, temp1;
  int i,j,sign;
  
  /* useless: for (i=0; i<dimx+1; i++) wtrans_ind[i] = 0; */
  
  sign = SIGN(wtrans[1] - wtrans[0]);
  temp1 = 0.; 
  for( j=0, i=2; i<dimx; i++ ) {  
    if ((fabs(wtrans[i] - wtrans[i-1]) > 0.0) && 
	((sign == 1) && ((SIGN(wtrans[i] - wtrans[i-1])) == -1))) {
      temp = wtrans[i-1];
            
      n_max++;
      /*  maximum condition */
      wtrans[j]=temp;
      wtrans_ind[j]=i-1;
      temp1=temp;
      j++;
    }
    sign=SIGN(wtrans[i]-wtrans[i-1]);
  } /* end loop on the filter "for( j=0,..." */
  
  return n_max;
} // end of wtmm1D_extrema_find


/***************************************************************************/
int pf1D_extrema_track( double *wtrans, int *wtrans_ind, 
		      double *maxsig, int *maxsig_ind,
		      int jj, int  n_max,
		      int nq, double *qArray, int nsc, double *scArray, 
		      double *ExtWTlis, double *Z ) {
  /***************************************************************************
   * Check whether the local maximum is located close to a maximum at a 
   * smaller scale. if yes connect both maxima, otherwise cancel it. 
   *
   * Parse:
   *    - wtrans, wtrans_ind : 1d arrays storing resp. the values and the
   *      indexes of the WT at the current scale,
   *    - maxsig, maxsig_ind : 1d arrays storing resp. the values and the
   *      indexes of the WT at the previous scale,
   *    - n_max : number of local extrema,
   *    - jj : index of the last detected extrema in the list of extremas
   *      already selected.
   *    - nq, qArray : resp. no of moments and list of moments,
   *    - nsc, scArray : resp. no of sceles and list of scales,
   * Outputs:
   *    - Z : 1d tabular (indexed by moment) storing the factors of the 
   *      partition function at current scale,
   *    - ExtWTlis : 1d array storing the WTMM over selected extrema
   *      at current scale,
   * Returns the final number of tracked global extrema at current scale. 
   *
   * Called by : wtmm1D_pf  
   ***************************************************************************/
  int i1=0, i2=0;
  int n_ext=0,ind_q;
  
  /* With jj, check that the number of maxima at larger scales is less or 
   * equal to that at a smaller scale. */
  while (((i1-1) < n_max) && ((i2-1) < jj))   {
    
    /* Choose the supremum between all previous values at smaller scales. */
    /* Activated only if the FLAG_SUPWTMM option is given                     */
    
#ifdef _PARSE_FRACTAL_PARAMETERS_
    IF(p_wtmm->flag_sup_wtmm)  
#else
    IF(FLAG_SUPWTMM)  
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
      {
	if ((wtrans_ind[i1] - maxsig_ind[i2]) <=
	    (maxsig_ind[i2+1] - wtrans_ind[i1]))
	  wtrans[i1] = fMax( wtrans[i1], maxsig[i2] );
	else
	  wtrans[i1] = fMax( wtrans[i1], maxsig[i2+1] );
      }
    
    /* Store in ExtWTlis the list of the values of the WT over the detected 
     * maxima, i.e. the final WTMM */
    ExtWTlis[n_ext++] = wtrans[i1];
    
    /* ...and also compute directly the factors of the partition funtion
     * from the WTMM */
    for ( ind_q=0; ind_q<nq; ind_q++ )
      Z[ind_q] += pow( wtrans[i1], qArray[ind_q]);
    /* Note that in the case of several signals (NSERIES>1), we simply
     * add the factors of the partition function */
    
    i1++, i2++;
    while ((i2 < jj) && (wtrans_ind[i1] >= maxsig_ind[i2]))
      i2++;
    i2--;
  }	  
  
  return n_ext;
} // end of pf1D_extrema_track



/**
 ** METHOD I : direct computation of the density spectrum through the 
 **            exponent tau(q) and the partition function
 **/


/***************************************************************************/
int spectrum1D_direct( double **Z, int dimx, 
		       int nq, double *qArray, int nws, double *wsArray,
		       double *Tauq, double *H, double *Dh ) {
  /***************************************************************************
   * Direct method to compute the spectrum thanks to the Legendre transform
   *
   * Parse:
   *    - Z : 2d tabular (scale by moment) storing the factors of the
   *      partition function,
   *    - dimx : size of the original signal,
   *    - nq, qArray : resp. no of moments and list of moments,
   *    - nws, wsArray : resp. no of scales and list of scales,
   * Outputs:
   *    - Tauq : multifractal exponents,
   *    - H : singularity exponents,
   *    - Dh : spectrum.
   *
   * The Legendre transform is approximated through the relations:
   *            h = \frac{d\tau}{dq} 
   *            D(h) = q h - tau(q) 
   *
   * Note : not implemented in the LastWave toolbox.
   ***************************************************************************/
  double *sumpf, *sumscpf;
  int i, count=0, ind=0;
  double sumsc=0., sumsqsc=0.;
  int ind_q, ind_ws;
  double min_q=qArray[0],  /* minimum considered moment */
    max_q;
  double Lws, LZ;
  int ind_ext; /* index of extrema */
  double wtext; /* local variable to store the value of the wavelet transform
		 * over the extrema */
   
  /* Local allocations */
  TrackNullAlloc( sumpf=(double*)calloc(nq,sizeof(double)) );
  TrackNullAlloc( sumscpf=(double*)calloc(nq,sizeof(double)) );

  /* Initialization of both local tabulars */
  for(ind_q=0; ind_q<nq; ind_q++) sumpf[ind_q] = sumscpf[ind_q] = 0.;
  
  for (ind_ws=1; ind_ws<nws; ind_ws++) {
    /* log of the current scale */
    Lws = log(wsArray[ind_ws])/log(10.);

      /* number of scales considered till now */
    count ++;
    /* sum of scales */
    sumsc += Lws;
    /* sum of squared scales */
    sumsqsc += Lws*Lws;
    for( ind_q=0; ind_q<nq; ind_q++ ) {
	/* log of the partition function */
	LZ = log( Z[ind_ws][ind_q])/log(10.);
	/* sum of the partition function */
	sumpf[ind_q] += LZ;
	/* sum of partition function weightened by the scale  */
	sumscpf[ind_q] += Lws * LZ;	
    }
  } /* end loop over the scales: "for (ind_ws=0;..." */
  
   /** Compute the tauq */
  for(ind_q=0; ind_q<nq; ind_q++) {
    /** Compute the tauq exponent corresponding to qArray[ind_q]
     * */
    Tauq[ind_q] = (sumscpf[ind_q]*count - sumsc*sumpf[ind_q]) 
      / (count*sumsqsc - sumsc*sumsc);
  } /* end of the 1st loop over the moments: "for (ind_q=0;..." */

 /* Compute the spectrum (h,D(h)) */
  for(ind_q=1; ind_q<nq-1; ind_q++) {

    /** First compute the singularity exponent h:
     *        h = \frac{d\tau}{dq}  */
    H[ind_q] = (Tauq[ind_q+1]-Tauq[ind_q-1])/(qArray[ind_q+1]-qArray[ind_q-1]);
    /* (tau[iq+1]-tau[iq-1])/(qArray[iq+1]-qArray[iq-1]) is an approximation
     * of the derivative \frac{d\tau}{dq} in iq */

    /** Then compute the density spectrum:
     *        D(h) = q h -tau(q) */
    Dh[ind_q] = qArray[ind_q]*H[ind_q] - Tauq[ind_q]; 
    /* variant for numerical percision : */
    /* Dh[ind_q] = qArray[ind_q] / (qArray[ind_q+1]-qArray[ind_q-1])
     * (Tauq[ind_q+1]-Tauq[ind_q-1]) - Tauq[ind_q]; */
    
  } /* end of the 2nd loop over the moments: "for (ind_q=0;..." */
  
  return OK;
} // end of spectrum1D_direct 



/**
 ** METHOD II : canonical computation through the canonical formula for
 **             the partition function and averaged exponents 
 **/



/***************************************************************************/
int mean1D_canon( double **ExtWTlis, int dimx, int *n_ext,
		  int nq, double *qArray, int nws, double *wsArray,
		  double **sTq, double **sTqLogT 
		  /* double **logSTq */ ) {
  /***************************************************************************
   * Function used in the canonical method proposed to compute the spectrum 
   * through to the Legendre transform - part I
   * Averages of multifractal are computed for the different scales.
   *
   * Parse:
   *    - ExtWTlis : 2d tabular (scale by moment) storing the WTMM, i.e.  
   *      the values of the WT computed over selected extrema,
   *    - n_ext : number of extrema at each scale, 
   *    - dimx : size of the original signal,
   *    - nq, qArray : resp. no of moments and list of moments,
   *    - nws, wsArray : resp. no of sceles and list of scales.
   * Outputs:
   *    - sTq, sTqLogT : 2d (scale by moment) intermedary tabulars used
   *      further to compute the different multifractal exponents and
   *      parseeterized similarly to the factors of the partition function.
   *      They look like:
   *         sTq[ws,q] = \sum{x \in L(ws)} |WT(x,ws)|^q
   *         sTqLogT[ws,q] = \sum{x \in L(ws)} |WT(x,ws)|^q * log(|WT(x,ws)|)
   *         logSTq[ws,q] =  log(\sum{x \in L(ws)} |WT(x,ws)|^q)
   *
   * The Legendre transform is approximated thanks to the method inspired
   * by Chhabra and Jensen: 
   *       \tilde{T_\psi} [s](q,x,a)=
   *                         \frac{T_\psi[s](x,a)}{Z(q,a)}
   * then the averages below are computed:
   *       < h >(q,a) = \sum_{(x,a)} 
   *                    \tilde{T_\psi} [s](q,x,a) \ln |T_\psi[s](x,a)|
   *       Dh(q,a) = \sum_{(x,a)} 
   *                 \tilde{T_\psi} [s](q,x,a) \ln \tilde{T_\psi}[s](x,a)
   * and, finally, the slopes of these quantities will provide the 
   * multifractal exponents (see spectrum1D_canon below).
   *
   * Note: in LastWave toolbox, this is the goal of the function
   *                PFComputeOneScaleF (one scale only)
   * in file pf_lib.c (package package_wtmm) called by: 
   *                ComputePartFuncOnExtrep (all scales)
   * in file pf_functions (package package_wtmm).
   ***************************************************************************/
  double *tempTq,*tempLogT;
  double *tempWT; 
  int ind_ws, ind_q, i, imin;
  int maxsize=-1, size;
  double tm, q;
  
  /* Find the maximum number of extrema over scales for further useful
   * allocations */
  for( i=0; i<nws; i++ )   maxsize = Max( maxsize, n_ext[i] );
  
  /* Local allocation of space for tempTq and tempLogT */
  TrackNullAlloc( tempTq=(double*)malloc(2*maxsize*sizeof(double)) );
  tempLogT = tempTq + maxsize;
  
  /* DEBUG   for( ind_ws=0; ind_ws<nws; ind_ws++ )... */
  for( ind_ws=1; ind_ws<nws; ind_ws++ ) {
    
    /* pointer on the list of selected extrema at scale ws */
    tempWT = ExtWTlis[ind_ws];
    /* number of extrema at this scale */
    size = n_ext[ind_ws];

    /* Rearrange in increasing order the values of the WTMM over the extrema */
    qsort((void*)tempWT,size,sizeof(double),
	  (int (*)(const void*,const void*))dumcompare);
    /* check that it is the same qsort for your compilator: it may depend
     * on the libraries used... */

    /* Find the first occurence of non null WTMM */
    imin = 0;
    while(imin<size && tempWT[imin] == 0.)  imin++;
    
    /* Approximation of the Legendre transform */
    for(ind_q=0; ind_q<nq; ind_q++) {
      
      /* First determine the current q */
      q = qArray[ind_q];
      
      /* Do we want T/tm to be >= 1 or <= 1 ? */
      if(q >= 0.)	tm = tempWT[imin];
      else	tm = tempWT[size-1];
      /* DEBUG tm=1.; */
      
      /* We compute Tq and Log(T/tm) */
      for( i=imin; i<size; i++ )      {
	tempTq[i] = pow( tempWT[i], q );
	tempLogT[i] = log(tempWT[i]/tm)/log(10.);
      }
      
      /** Note that the sum below allow to compute statistical variables 
       ** over several signals (NSERIES>1) by simply adding their values **/

      /* We compute sTq and sTqLogT */
      if(q >= 0.)      
	for( i=imin; i<size; i++ )	{
	  sTq[ind_ws][ind_q] += tempTq[i]; /* pow( tempWT[i], q ); */
	  sTqLogT[ind_ws][ind_q] += tempTq[i] * tempLogT[i];
	  /* pow( tempWT[i], q ) *log( tempWT[i]/tm )/log(10.); */
	}
      
      else 
	for( i=size-1; i>=imin; i-- )	{
	  sTq[ind_ws][ind_q] += tempTq[i];/* pow( tempWT[i], q ); */
	  sTqLogT[ind_ws][ind_q] += tempTq[i] * tempLogT[i]; 
	  /* pow( tempWT[i], q ) * log( tempWT[i]/tm )/log(10.); */
	}
      
      /* sTqLogT = sTqLogT + log(tm)*sTq */
      sTqLogT[ind_ws][ind_q] += log(tm)/log(10.)*sTq[ind_ws][ind_q];
      
      /* We compute LogSTq
	 if(sTq[ind_ws][ind_q] != 0.0)	{
	 logSTq[ind_ws][ind_q] = 
	 log(sTq[ind_ws][ind_q]/((double) (size-imin)))/log(10.);
	 }    else 	{
	 logSTq[ind_ws][ind_q] = 0.;
	 } 
      */ 
      
    } /* end loop over the moments: "for (ind_q=0;..." */
  } /* end loop over the scales: "for (ind_ws=0;..." */
  
  if(tempTq) free(tempTq); /* then tempLogT is automatically free */
  
  return OK;
} // end of mean1D_canon


/***************************************************************************/
int spectrum1D_canon( double **sTq, double **sTqLogT, /* double **logSTq */
		      int dimx, 
		      int nq, double *qArray, int nws, double *wsArray,
		      double *Tauq, double *H, double *Dh ) {
  /***************************************************************************
   * Function used in the canonical method proposed to compute the spectrum 
   * through to the Legendre transform - part II
   * Estimation of multifractal exponents are realized through linear
   * regression over the scales.
   *
   * Parse:
   *    - sTq, sTqLogT : 2d (scale by moment) tabulars parseeterized
   *      similarly to the factors of the partition function,
   *    - dimx : size of the original signal,
   *    - nq, qArray : resp. no of moments and list of moments,
   *    - nws, wsArray : resp. no of scales and list of scales,
   * Outputs:
   *    - Tauq : multifractal exponents,
   *    - H : singularity exponents,
   *    - Dh : spectrum of singularity.
   *
   * Note: in LastWave toolbox, this operation is realized through the 
   * functions 
   *        tauqSpectrum and singSpectrum
   * of the script file wtmm1d.pkg (package scripts). 
   * These scripts make an implicit use of the functions:
   *        PFAccessTQFloat    PFAccessHQFloat    PFAccessDQFloat
   * of file pf_lib.c (package package_wtmm), and of the function :
   *        LineFitSig
   * of file signal_function.c (package package_signal), declared in the
   * file utl_stats.c.
   ***************************************************************************/
  int ind_ws, ind_q;
  double q, stq;
  double *tq, *hq, *dq;
  double *LwsArray;
  int ss;
  
  /* Local allocations */
  TrackNullAlloc( tq=(double*)calloc(3*nws,sizeof(double)) );
  TrackNullAlloc( LwsArray=(double*)calloc(nws,sizeof(double)) );
  
  hq = tq + nws;
  dq = tq + 2*nws;
  
  /* Array storing the log of the scales */
  for( ind_ws=0; ind_ws<nws; ind_ws++ )
    LwsArray[ind_ws] = log(wsArray[ind_ws])/log(10.);
  
  for(ind_q=0; ind_q<nq; ind_q++) { 
    
    /* consider the current moment for computing the variables */
    q=qArray[ind_q];
    
    for( ind_ws=0; ind_ws<nws; ind_ws++ ) { 
      stq = sTq[ind_ws][ind_q];
      
      /** Compute the canonical multifractal exponents tq */ 
      if(stq == 0.)	tq[ind_ws] = 0.;
      else   tq[ind_ws] = (double) log(stq)/log(10.);
      /* equivalent to the function PFAccessTQFloat */

      /** Compute the canonical singularity exponents hq:
       *      <h>(q,ws) = \sum_{(x,ws)} 
       *         \tilde{WT}(q,x,ws) \ln |WT(x,ws)|
       * where:
       *   \tilde{WT}(q,x,ws)= WT(x,ws) / Z(q,ws)
       * by using: 
       *     sTq[ws,q] = \sum{x \in L(ws)} |WT(x,ws)|^q
       *     sTqLogT[ws,q] = \sum{x \in L(ws)} |WT(x,ws)|^q * log(|WT(x,ws)|)
       */
      if(stq == 0.)	hq[ind_ws] = 0.;
      else   hq[ind_ws] = (double)(sTqLogT[ind_ws][ind_q] / stq);
      /* equivalent to the function PFAccessHQFloat in file pf_lib.c */
      
      /* Compute the canonical spectrum dq:
       *      dq(q,ws) = \sum_{(x,ws)} 
       *        \tilde{WT}(q,x,ws) \ln \tilde{WT}(x,ws)
       */
      if(stq == 0.)	dq[ind_ws] = 0.;
      else   dq[ind_ws] = (double)(q * sTqLogT[ind_ws][ind_q] / stq 
				   - log(stq)/log(10.));
      /* equivalent to the function PFAccessDQFloat */
      
    } /* end loop over the scales: "for (ind_ws=0;..." */
    
    
    /* Regression over the scales to get the different exponents.
     * For each moment, the tq, hq and dq are tabulars of nws values
     * whose slopes give the corresponding exponent Tauq, H and Dh. */
    ss = fit_line( tq, LwsArray, nws,  &(Tauq[ind_q]) );
    
    fit_line( hq, LwsArray, nws, &(H[ind_q]) );

    fit_line( dq, LwsArray, nws, &(Dh[ind_q]) );
    
  } /* end loop over the moments: "for (ind_q=0;..." */
  
  IFVERBOSE 
    WarningV("%d scales used for the estimation of multifractal exponents\n", ss);
  
  /* Free memory */ 
  if(tq) free(tq);
  if(LwsArray) free(LwsArray);
  
  return OK;
} // end of spectrum1D_canon


/***************************************************************************/
int dumcompare(const double *d1,const double *d2) {
  /***************************************************************************
   * Stupid function to compare double values and used by qsort in function
   * canonmean below.
   ***************************************************************************/
  
  if(*d1<*d2)    return -1;
  else if(*d1 == *d2)    return 0;
  else    return +1;
} // end of dumcompare
