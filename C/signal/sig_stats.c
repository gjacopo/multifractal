#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_alloc.h>      
#include <utl_parse.h>
#include <utl_operator.h>	

#include <flt_stats1d.h>
#include <flt_stats2d.h>
#include <flt_histo.h>

#include <signal.h>
#include <sig_alloc.h>
#include <sig_operator.h>
#include <sig_stats.h>

#ifdef DEBUG
static int cnt=1;
#endif

/***************************************************************************/
double rmse_signal(Signal* s1, Signal* s2, char**m) {
/***************************************************************************/

  double rmse=0.;
  Signal *s=NULL;
  
  TrackNullAlloc( s1 );   TrackNullAlloc( s1->signal );
  TrackNullAlloc( s2 );   TrackNullAlloc( s2->signal );

  TrackNullAlloc( s=create_signal(s1->xdim,s1->ydim) );
  copy_signal(s1,s,m);

  substract_signal(s2, s, m);
  cuadra_signal(s,m);
  rmse=media_signal(s,m);

  free_signal(s);

  return sqrt(rmse);
} // end of rmse_signal


/***************************************************************************/
int mediadisp_signal( Signal* s, double *med, double *disp, char **m ) {
  /***************************************************************************
   * Computes med(s) over all pixels
   ***************************************************************************/

  TrackNullAlloc( s ); TrackNullAlloc( s->signal );
  mediadisp(s->xdim, s->ydim, s->signal, med, disp, m );

  return OK;
} // end of mediadisp_signal


/***************************************************************************/
double media_signal( Signal* s, char **mask ) {
  /***************************************************************************
   * Computes med(s) over all pixels
   ***************************************************************************/

  double med;
  
  TrackNullAlloc( s ); TrackNullAlloc( s->signal );
  med = media(s->xdim, s->ydim, s->signal, mask );

  return med;
} // end of media_signal


/***************************************************************************/
double dispersion_signal( Signal* s, char **mask ) {
  /***************************************************************************
   * Computes sigma(s) over all pixels
   ***************************************************************************/

  double disp;
  
  TrackNullAlloc( s ); TrackNullAlloc( s->signal );
  disp = dispersion(s->xdim, s->ydim, s->signal, mask );
  
  return disp;
} // end of dispersion_signal


/***************************************************************************/
int extrema_signal(Signal* s, char **mask ) {
/***************************************************************************/

  TrackError( extrema( s->xdim, s->ydim, s->signal, s->mms, mask ),
	      "Error while computing signal extrema");
  s->mms[2] = s->mms[1] - s->mms[0];

  return OK;
} // end of extrema_signal


/***************************************************************************/
int display_extrema_signal(char *text, Signal* s, char **mask ) {
/***************************************************************************/

  extrema_signal( s, mask );

  return OK;
} // end of display_extrema_signal


/****************************************************************/
int compute_histo_signal(Signal *sig, char **mask) {
  /****************************************************************/
  
  if(sig->hist == NULL) return ERROR;
  else {
    if(sig->mms == NULL) extrema_signal(sig, mask );
    compute_histo( sig->xdim, sig->ydim, sig->signal, 	
		   sig->mms, sig->hist, mask );
  }
  
  return OK;
} // end of compute_histo_signal


/****************************************************************/
int equalize_histo_signal( Signal *s, int bin, char **mask ) {
  /****************************************************************/

  TrackNullAlloc( s );   TrackNullAlloc( s->signal );

  if(s->hist == NULL) 
    TrackNull( s->hist=init_histo(bin), 
	       "Error allocation while creating histogram" );
  
  TrackError( compute_histo_signal(s,mask), "Error computing histogram" );
  TrackError( norma_histo(s->hist), "Error normalizing histogram" );  
  TrackError( cum_histo(s->hist), "Error computing cumulated histogram" );
  
  TrackError( equalize_intensity(s->signal,s->xdim,s->ydim,s->mms,
				s->hist,mask), 
	      "Error processing intensity equalization" ); 
  // New signal is in the range [0 bin]
  
  return OK;
} // end of equalize_histo_signal


/****************************************************************/
int histmatch_signal( Signal *s1, Signal *s2, int bin,
		      char **mask ) {
  /****************************************************************
   * Apply histogram matching to s2 by remapping intensities in s2 
   * to have histogram specified by the histo of s1. 
   * Output goes to s2.
   ****************************************************************/
  
  Histo *h1=NULL, *h2=NULL;
  
 if(s1->hist == NULL)
   TrackNull( s1->hist=init_histo(bin), 
	      "Error allocation while creating histogram" );
 if(s2->hist == NULL)
   TrackNull( s2->hist=init_histo(bin), 
	      "Error allocation while creating histogram" );
 
  /* Compute (non normalized) histogram of frequencies of first
   * image */
  TrackError( compute_histo_signal(s1,mask), "Error computing histogram"  );
  /* ...and second image */
   TrackError( compute_histo_signal(s2,mask), "Error computing histogram"  );

  /* Compute extrema of the signal to be matched */
  TrackError( extrema_signal(s2, mask ), 
	      "Error computing extrema of signal" );

  /* Old note: while the histogram is computed for pixels of mask
   * only, the intensity matching is performed for all pixels 
   * of the input image. For that reason, a NULL mask is passed
   * as an argument to the matchintensity function */
  TrackError( match_intensity(s2->signal,s2->xdim,s2->ydim,s2->mms,
			     s2->hist,s1->hist,mask), 
			     "Error processing histogram matching"  ); 
  
  return OK;
} // end of histmatch_signal


#ifdef DEBUG

/****************************************************************/
Signal* test_histmatch( Signal *s, int bin, char **mask ) {
  /****************************************************************
   * Apply histogram matching to s1 by remapping intensities in s2 to
   * have histogram specified in histo. Output goes to s2.
   ****************************************************************/
  
  Signal *seq;
  double rmse;

  TrackNullAlloc( seq=create_signal(s->xdim,s->ydim) );

  /* Compute non normalized histogram (frequencies) */
  TrackNull( s->hist=init_histo(bin), 
	     "Error allocation while creating histogram" );
  TrackError( compute_histo_signal(s,mask), "Error computing histogram" );
 
  /* Equalize */
  TrackNull( copy_signal(s,seq,mask), 
	     "Error allocation while copying signal" );
  TrackError( equalize_histo_signal(seq, bin, mask), 
	      "Error processing histogram equalization" );
  
  /* Compute new equalized histogram */
  TrackNull( seq->hist=init_histo(bin) );
  TrackError( compute_histo_signal(seq,mask), "Error computing histogram" );
  
  TrackError( match_intensity(s->signal,s->xdim,s->ydim,s->mms,
			     s->hist,seq->hist,mask), 
	      "Error processing histogram matching" );

  TrackError( compute_histo_signal(s,mask), "Error computing histogram" );

  return s;
} // end of test_histmatch

#endif
