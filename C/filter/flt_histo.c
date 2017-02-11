#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_alloc.h>      
#include <utl_stats.h>      
#include <utl_operator.h>

#include <flt_stats1d.h>
#include <flt_stats2d.h>
#include <flt_histo.h>

#ifndef BIN_MODE_INF
#define BIN_MODE_INF
#endif

#ifdef BIN_MODE_INF
#define IRBIN(ir,bin,value) ( \
	( \
         (((ir)=(int)((bin)*(value)))>=0) && ((ir)<((bin))) \
          ) ? \
	 (ir) : ((bin)-1) \
         )
#else
#define IRBIN(ir,bin,value) ((int)((bin-1) * (value)))
#endif


/* see: compute_histo, equalize_intensity and match_intensity */

/****************************************************************/
Histo* alloc_histo() {
/****************************************************************/
  Histo *hist=NULL;
  
  hist=(Histo*)malloc(sizeof(Histo));
  return hist;
}


/****************************************************************/
Histo* create_histo(int bin) {
/****************************************************************/
  
  Histo *hist=NULL;
  
  if( (hist=alloc_histo()) == NULL ||
      (hist->h=(double*)calloc(bin,sizeof(double))) == NULL ||
      (hist->cum=(double*)calloc(bin,sizeof(double))) == NULL )
    return NULL;
  
  hist->bin = bin;
  hist->nsample = 0;

  return hist;
}

/****************************************************************/
int free_histo(Histo* hist) {
/****************************************************************/

  Free(hist->h);
  Free(hist->cum);
  Free(hist);

  return OK;
}

/****************************************************************/
int copy_histo( Histo* hin, Histo* hdest ) {
/****************************************************************/

  if(hdest == NULL)
    TrackNullAlloc( hdest=create_histo(hin->hin) );
  
  if(hin->h != NULL) 
    copy1D(hin->bin,hin->h,hdest->h,NULL);
  if(hin->cum != NULL) 
    copy1D(hin->bin,hin->cum,hdest->cum,NULL);

  return OK;
}

/****************************************************************/
int compute_histo( int dimx, int dimy, double **s, 	
		   double *mms, Histo *hist, char **mask ) {
  /****************************************************************/
  register int i, j, ir, bin;
  int flag=FALSE;
  
  TrackNullAlloc( hist );  TrackNullAlloc( hist->h );

  if( mms == NULL) {
    flag = TRUE;
    TrackNullAlloc( mms=(double*)malloc(3*sizeof(double)) );
    extrema( dimx, dimy, s, mms, mask );
    mms[2] = mms[1] - mms[0];
  }
  
  /* First step: vectors are cleaned up */
  bin = hist->bin;
  for(ir=0; ir<bin; ir++)  hist->h[ir] = hist->val[ir] = 0.;
  
  /* Second step: Data are explored, taking the mask into account if required */
  if(mask == NULL) {
    for(i=0; i<dimx; i++)
      for(j=0; j<dimy; j++) {
	ir = IRBIN( ir, bin, (s[j][i]-mms[0])/mms[2]);
	hist->h[ir] += 1.;
	hist->val[ir] += s[j][i];
      }
    hist->nsample = dimx * dimy;
    
  } else
    for(i=0; i<dimx; i++)
      for(j=0; j<dimy; j++) 
	if(mask[j][i] == TRUE) {
	  ir = IRBIN( ir, bin, (s[j][i]-mms[0])/mms[2]);
	  hist->h[ir] += 1.;
	  hist->val[ir] += s[j][i];
	  hist->nsample++;
	}
  
  IF(flag) {
    free(mms);
    mms = NULL;
  }
  
  return OK;
}


/****************************************************************/
int norma_histo(Histo *hist) {
  /****************************************************************/
  
  register int ir;
  
  for (ir=0; ir<hist->bin; ir++) {
    hist->h[ir] /= (double)hist->nsample;
    hist->cum[ir] /= (double)hist->nsample;
  }

  hist->flnorma = TRUE;

  return OK;
}


/****************************************************************/
int cum_histo( Histo *hist ) {
  /****************************************************************/

  register int ir;
  double sum=0.,total;

  TrackNullAlloc( hist->h);
  
  if(hist->flnorma = TRUE)    total=1.;
  else                        total=(double)hist->nsample;

  for (ir=0; ir<hist->bin; ir++) {
    sum += hist->h[ir];
    hist->cum[ir] = sum;
    if (sum > total) break;
  }
  for( /* we start where we stopped */; ir<hist->bin; ir++ ) 
    hist->cum[ir] = total;
  
  return OK;
}


/***************************************************************************/
int scale_histogram( Histo *h0, double sc0, Histo *h1, double sc1 ) {
  /***************************************************************************
   * Adapt (for comparison) the histograms coming from data with different 
   * resolutions sc1 to a common basic resolution sc0.
   * Note that the transformation rule is general and would be valid for
   * any quantity exhibiting a general power-law scaling in the histogram.
   * In addition, notice that if sc0=sc1, Nnew=h1[ih].
   * So that, this routine does not induce any change in the data if they are
   * completely homogeneous in size. Notice also that, as the dependence in
   * sc0, sc1 is logarithmic, the induced change is very small if data do
   * not differ in size by a considerable amount.  
   ***************************************************************************/
   double hmode=0.;
   double Nnew;
   int ih;

   if((bin=h0->bin) != h1->bin) Error("Histograms should be of same size");

   /* Obtain the histogram mode */
   for (ih=0;ih<bin;ih++) hmode = fMax(hmode,h1);
   
   /* We now modify the histogram taking the following into account:
    * h1(h)= hmode * r^{d-D(h)}, and now r must change from sc1 (that
    * of data) to sc0 (the new basic resolution). The average h associated
    * to the bin box must also be changed accordingly. */
   for( ih=0; ih<bin; ih++ ) 
     if(h1[ih] > 0.5) { // there is at least one case
       
       Nnew = hmode * exp(log(h1[ih]/hmode) * log(sc0) / log(sc1));
       
     /* The values are now updated */
       h0[ih] *= (Nnew / h1[ih]);
       h1[ih] = Nnew;
     }
   
   return OK;
} // end of scale_histogram


/****************************************************************/
int equalize_intensity( double** s, int xdim, int ydim,
			double *mms, Histo *hist, char **mask) {
/****************************************************************
   * Calculate the mapping to the new intensity value.
   ****************************************************************/
 
  register int ir, i, j;
  int bin=hist->bin;
  int *intens=NULL;

  TrackNullAlloc( intens=(int*)malloc(bin*sizeof(int)) );

  if(hist->flnorma = FALSE)      
    TrackError( norma_histo(hist), "Error normalizing histogram" );

  for (ir=0; ir<bin; ir++) {
    /* scaleup and round to nearest integer */
    intens[ir] = (int) (hist->cum[ir] * bin+ 0.5);	
    if (intens[ir] > bin) intens[ir] = bin;
  }
  
  if(mask == NULL) 
    for(i=0; i<xdim; i++)
      for(j=0; j<ydim; j++) {
	ir = IRBIN( ir, bin, (s[j][i]-mms[0])/mms[2]);
	s[j][i] = (double)intens[ir];
      }
  else
    for(i=0; i<xdim; i++)
      for(j=0; j<ydim; j++) 
	if(mask[j][i] == TRUE) {
	  ir = IRBIN( ir, bin, (s[j][i]-mms[0])/mms[2]);
	  s[j][i] = (double)intens[ir];
	}

  return OK;
}




/****************************************************************/
int match_intensity( double** s, int xdim, int ydim, double *mms, 
		    Histo *hist, Histo *hfit,
		    char **mask ) {
  /****************************************************************/
  
  double scale;
  int *left=NULL, *right=NULL;
  int bin=hist->bin;
  register int ir, i, j;
  int IR=0,p;
  int hsum=0;
  
  TrackNullAlloc( left=(int*)calloc(bin,sizeof(int)) );
  TrackNullAlloc( right=(int*)calloc(bin,sizeof(int)) );
  
  /* normalize h2 to conform with number of samples of h1 */
  scale = (double)hist->nsample / hfit->nsample;
  TrackError( op_scale(bin,1,scale,&(hfit->h),NULL),
	      "Error scaling signal" ); 

  /* Evaluate remapping of all input gray levels;
   * Each input gray value maps to an interval of valid output values.
   * The endpoints of the intervals are left[] and right[] */
  for(ir=0; ir<bin; ir++) {
    left[ir] = IR;	/* left end of interval */
    hsum += hist->h[ir];	/* cumulative value for interval */
    /* compute width of interval */
    while(hsum>hfit->h[IR] && IR<bin-1) { 
      hsum -= hfit->h[IR];	/* adjust hsum as interval widens */
      IR++;		/* update */
    }
    right[ir] = IR;	/* init right end of interval */
  }

  /* Clear h1 and reuse it below */
  for(ir=0; ir<bin; ir++) hist->h[ir] = 0;
  
  /* Revisit all input pixels */
 if(mask == NULL) {
    for(i=0; i<xdim; i++)
      for(j=0; j<ydim; j++) {
	ir = IRBIN( ir, bin, (s[j][i]-mms[0])/mms[2]);
	p = left[ir];
	if(hist->h[p] >= hfit->h[p])	/* mapping doesnt satisfy hfit */
	  p = left[ir] = Min(p+1, right[ir]);
	s[j][i] = (double)p;
	hist->h[p]++;
      }
  
 }  else {
    for(i=0; i<xdim; i++)
      for(j=0; j<ydim; j++) 
	if(mask[j][i] == TRUE) {
	  ir = IRBIN( ir, bin, (s[j][i]-mms[0])/mms[2]);
	  p = left[ir];
	  if(hist->h[p] >= hfit->h[p])	/* mapping doesnt satisfy hfit */
	    p = left[ir] = Min(p+1, right[ir]);
	  s[j][i] = (double)p;
	  hist->h[p]++;
	}
}
  
  return OK;
}


/****************************************************************/
int add_histo(Histo *h1, Histo *h2) {
  /****************************************************************/
  
  register int ir;
  double nh1, nh2, newnh1; 
  
  nh1 = (h1->flnorma==TRUE) ? (double)h1->nsample : 1.;
  nh2 = (h2->flnorma==TRUE) ? (double)h2->nsample : 1.;
  
  h1->nsample += h2->nsample;
  newnh1 = (h1->flnorma==TRUE) ? (double)h1->nsample : 1.;

  for (ir=0; ir<Min(h1->bin,h2->bin); ir++) {
    h1->h[ir] =  nh1 * h1->h[ir] + nh2 * h2->h[ir];
    h1->h[ir] /= newnh1;
  }
  
  return OK;
}

/****************************************************************/
int sub_histo(Histo *h1, Histo *h2) {
  /****************************************************************/
  
  register int ir;
  double nh1, nh2, newnh1; 
  
  nh1 = (h1->flnorma==TRUE) ? (double)h1->nsample : 1.;
  nh2 = (h2->flnorma==TRUE) ? (double)h2->nsample : 1.;
  
  h1->nsample += h2->nsample;
  newnh1 = (h1->flnorma==TRUE) ? (double)h1->nsample : 1.;
  
  for (ir=0; ir<Min(h1->bin,h2->bin); ir++) {
    h1->h[ir] =  nh1 * h1->h[ir] - nh2 * h2->h[ir]; 
    h1->h[ir] /= newnh1;
  }
  
  return OK;
}

/*
unsigned char ThresholdOtsu(unsigned char** data, int nRows, int nCols)
 {
  int i;
  int threshold;
  double totMean;	// mean gray-level for the whole image 
  double* bcv;		// between-class variance 
  double* histo;	// normalized histogram 
  double* cmf;		// cumulative normalized histogram (probability mass function) 
  double* mean;		// mean gray-level 
  
  histo = CreateNormalizedHistogram(data, nRows, nCols);

  cmf = (double *) malloc(256 * sizeof(double));
  mean = (double *) malloc(256 * sizeof(double));

  cmf[0] = histo[0];
  mean[0] = 0.0;

  for(i = 1; i < 256; i++)
   {
    cmf[i] = cmf[i - 1] + histo[i];
    mean[i] = mean[i - 1] + i * histo[i];
   }

  totMean = mean[255];

  bcv = (double *) malloc(256 * sizeof(double));

  // when i = 255 the denominator of the bcv expression 
  // becomes 0.0 since cmf[255] = 1.0
  bcv[255] = 0.0;
  threshold = 255;

  // calculate the bcv at each gray-level      and find the maximum 
  for(i = 0; i < 256 - 1; i++)
   {
    if(ZERO(cmf[i]) == FALSE && ZERO(1.0 - cmf[i]) == FALSE)	
      // guard against zero denominator
     {
      bcv[i] = SQUARE(totMean * cmf[i] - mean[i]) / (cmf[i] * (1.0 - cmf[i]));
      
      if(bcv[i] > bcv[threshold])
       threshold = i;
     }
   }

  return (unsigned char) threshold;
 }

*/

