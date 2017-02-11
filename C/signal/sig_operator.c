#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Librairies  ROUTINES */

#include <utils.h>
#include <utl_alloc.h>      
#include <utl_operator.h>	

#include <flt_stats1d.h>
#include <flt_stats2d.h>
#include <flt_histo.h>

#include <signal.h>
#include <sig_alloc.h>
#include <sig_operator.h>
#include <sig_stats.h>


/***************************************************************************/
int compare_dimmask( Mask *m1, Mask *m2) {
  /***************************************************************************/
  
  if(m1->xdim!=m2->xdim || m1->ydim!=m2->ydim)
    return ERROR;
  else 
    return OK;
} // end of compare_dimsignal


/***************************************************************************/
int compare_dimsignal( Signal *s1, Signal *s2) {
  /***************************************************************************/
  
  if(s1->xdim!=s2->xdim || s1->ydim!=s2->ydim)
    return ERROR;
  else 
    return OK;
} // end of compare_dimsignal


/***************************************************************************/
int compare_dimimage( Image *im1, Image *im2) {
  /***************************************************************************/
  
  if(im1->xdim!=im2->xdim || im1->ydim!=im2->ydim || im1->zdim!=im2->zdim)
    return ERROR;
  else 
    return OK;
} // end of compare_dimimage


/***************************************************************************/
int affine_signal( Signal* s, double scale, double shift, char **m ) {
/***************************************************************************/
  
  TrackNullAlloc( s ); 
  // TrackError( 
  op_affine(s->xdim,s->ydim,scale,shift,s->signal,m);
  //	      "Error operating affine transformation" );

  return OK;
} // end of affine_signal


/***************************************************************************/
int scale_signal( Signal* s, double scale, char **m ) {
/***************************************************************************/
  
  TrackNullAlloc( s ); 
  // TrackError( 
  op_scale(s->xdim,s->ydim,scale,s->signal,m);
  //	      "Error operating scale transformation" );

  return OK;
} // end of scale_signal


/***************************************************************************/
int fill_signal( Signal* s, double val, char **m ){
  /***************************************************************************/
  
  //  TrackError( 
  fill( s->xdim, s->ydim, val, s->signal, m );
  //      "Error filling signal" );
 
  return OK;
} // end of fill_signal


/***************************************************************************/
int substract_signal( Signal* sin, Signal* sout, char **m ) {
  /***************************************************************************
   * Compute sout-sin and stores the result in sout
   ***************************************************************************/
  int dimx, dimy;

  TrackNullAlloc( sin ); TrackNullAlloc( sin->signal );
  TrackNullAlloc( sout ); TrackNullAlloc( sout->signal );

  dimx = Min(sin->xdim,sout->xdim);
  dimy = Min(sin->ydim,sout->ydim);

  //  TrackError( 
  op_diff(dimx, dimy, sin->signal, sout->signal,m);
  //	      "Error operating diffence of signals" );

  return OK;
} // end of substract_signal


/***************************************************************************/
int add_signal( Signal* sin, Signal* sout, char **m ) {
  /***************************************************************************
   * Compute sout+sin and stores the result in sout
   ***************************************************************************/
  int dimx, dimy;

  TrackNullAlloc( sin ); TrackNullAlloc( sin->signal );
  TrackNullAlloc( sout ); TrackNullAlloc( sout->signal );

  dimx = Min(sin->xdim,sout->xdim);
  dimy = Min(sin->ydim,sout->ydim);

  //  TrackError( 
  op_add(dimx, dimy, sin->signal, sout->signal,m);
  //	      "Error operating addition of signals" );

  return OK;
} // end of add_signal


/***************************************************************************/
int divide_signal( Signal* sin, Signal* sout, char **m ) {
  /***************************************************************************
   * Compute sout/sin and stores the result in sout
   ***************************************************************************/
  int dimx, dimy;

  TrackNullAlloc( sin ); TrackNullAlloc( sin->signal );
  TrackNullAlloc( sout ); TrackNullAlloc( sout->signal );

  dimx = Min(sin->xdim,sout->xdim);
  dimy = Min(sin->ydim,sout->ydim);

  op_divide(dimx, dimy, sin->signal, sout->signal,m);
  //	      "Error operating division of signals" );

  return OK;
} // end of divide_signal


/***************************************************************************/
int cuadra_signal( Signal* s, char **m ) {
  /***************************************************************************
   * Computes sqr(s) for each pixel
   ***************************************************************************/

  TrackNullAlloc( s ); TrackNullAlloc( s->signal );
  op_square(s->xdim, s->ydim, s->signal, s->signal,m);

  return OK;
} // end of cuadra_signal


/***************************************************************************/
double fabs_signal( Signal* s, char **m ) {
  /***************************************************************************
   * Computes abs(s) over all pixels
   ***************************************************************************/

  TrackNullAlloc( s ); TrackNullAlloc( s->signal );
  op_fabs(s->xdim, s->ydim, s->signal, s->signal,m);

  return OK;
} // end of fabs_signal


/***************************************************************************/
int norma_signal( Signal* s, char **m ) {
/***************************************************************************/
  
  extrema_signal( s, m ); // keep in memory old borns of the signal
  norma( s->xdim, s->ydim, s->signal, m );
  
  return OK;
} // end of norma_signal


/***************************************************************************/
int threshold2ptr_signal( Signal* sig, int sign, double thres, double* s ) {
  /***************************************************************************/
  
  int N;
  int dimx=sig->xdim, dimy=sig->ydim;
  
  if(s == NULL) 
    TrackNullAlloc( s=(double*)malloc(dimx*dimy*sizeof(double)) );
  
  N = op2ptr_threshold( dimx, dimy, sig->signal, s, sign, thres );
  
  return N;
} // end of threshold2ptr_signal


/***************************************************************************/
int mask2ptr_signal(Signal* sig, int sign, double thres, char* s ) {
  /***************************************************************************/
  
  int N;
  int dimx=sig->xdim, dimy=sig->ydim;
  
  if(s == NULL) 
    TrackNullAlloc( s=(char*)malloc(dimx*dimy*sizeof(char)) );
  
  N = op2ptr_binary( dimx, dimy, sig->signal, s, sign, thres );
  
  return N;
} // end of mask2ptr_signal


/***************************************************************************/
int mask_signal(Signal* sig, int sign, double thres, char**m ) {
  /***************************************************************************/
  
  int N;
  int dimx=sig->xdim, dimy=sig->ydim;
  
  if(m == NULL)    
    TrackNullAlloc( m=cmatrix2D(dimy,dimx) );
  
  /* if sign>0, check:      signal>thres
   * i.e: if signal>thres, then mask=TRUE
   *      else                  mask=FALSE
   * if sign<0, check:      -signal>-thres <=> signal<thres
   * i.e. the opposite    
   */
  N = op_binary( dimx, dimy, sig->signal, m, sign, thres );
  
  return N;
} // end of mask_signal


/***************************************************************************/
int vratio_signal( Image* im, Signal *s, int ic1, int ic2, double c, 
		   char**m ) {
  /***************************************************************************
   * Computes normalized ratio of vectorial components IC1 and IC2 of the  
   * image IM.
   * Returns in signal S the ratio:
   *        S = c*IC2 / IC1
   ***************************************************************************/
  
  int dimx=im->xdim, dimy=im->ydim, dimv=im->vdim;
  Signal *s1;
  
  if(ic1>dimv || ic2>dimv) 
    Error("Index outside the range of possible vectorial components");
  
  if(s == NULL)        TrackNullAlloc( s=create_signal(dimx,dimy) );
  
  /* allocate pointer, do not allocate buffer */
  TrackNullAlloc( s1=create_psignal(dimx,dimy) );
  TrackError( read_im2sig(im,ic1,s1), 
	      "Error reading input channel" ); /* s1 points to channel ic1 */
  
  /* copy of buffer: content of channel ic2 is copied in s */
  TrackError( read_im2sig(im,ic2,s), 
	      "Error reading input channel" );  
  
  TrackError( divide_signal(s1,s,m), 
	      "Error operating channel division" ); /* s = s/s1 */ 
  
  if(c != 1.)
    TrackError( scale_signal(s,c,m),
		"Error operating scale factor" ); /* s = c*s */ 
  
  /* Make s1 point on NULL to free the memory */
  s1->signal = NULL;
  free_signal(s1);
  
  return OK;
} // end of vratio_signal


/***************************************************************************/
int vdiffratio_signal( Image* im, Signal *s, int ic1, int ic2, 
		       double c, char**m ) {
  /***************************************************************************
   * Computes normalized ratio of channels stored in components IC1 and IC2  
   * of the image IM.
   * Returns in signal S the ratio:
   *        S = (IC2 - c*IC1) / (IC1 + c*IC2)
   ***************************************************************************/
  
  int dimx=im->xdim, dimy=im->ydim, dimv=im->vdim;
  Signal *s1, *s2;
  
  if(ic1>dimv || ic2>dimv) 
    Error("Index outside the range of possible vectorial components");
  
  if(s == NULL)       TrackNullAlloc( s=create_signal(dimx,dimy) );
  
  /* decide if allocate pointers or matrix depending on the value of c */
  if(c != 1.) {/* allocate buffer */
    TrackNullAlloc( s1=create_signal(dimx,dimy) );
  } else       /* use pointer */
    TrackNullAlloc( s1=create_psignal(dimx,dimy) );
  /* use pointer or copy buffer */
  TrackError( read_im2sig(im,ic1,s1),
	      "Error reading input channel" ); 
  
  /* allocate pointer, do not allocate buffer */
  TrackNullAlloc( s2=create_psignal(dimx,dimy) );
  TrackError( read_im2sig(im,ic2,s2),
	      "Error reading input channel" ); /* s2 points to channel ic2 */
  
  /* copy of buffer: content of channel ic2 is copied in s */
  TrackError( read_im2sig(im,ic2,s),
	      "Error reading input channel" );  
  
  /* Operations */
  if(c != 1.)
   TrackError( scale_signal(s1,c,m),
		"Error operating scale factor" ); /* s1 = c*s1 */ 
 
  TrackError( substract_signal(s1,s,m),
	      "Error operating channel difference" ); /* s = s - c*s1 */ 
  TrackError( add_signal(s1,s2,m),
	      "Error operating channel addition" ); /* s2 = c*s1 + s2 */ 
  
  TrackError( divide_signal(s2,s,m),
	      "Error operating channel division" ); /* s = s/s2 */ 
  /* Final result is: (s2-c*s1) / (s2+c*s1) */
  
  /* Make them point on NULL to free the memory */
  s2->signal = NULL;
  if(c == 1.) s1->signal = NULL;
  free_signal(s1);
  free_signal(s2);
  
  return OK;
} // end of vdiffratio_signal


/***************************************************************************/
int rvdiffratio_signal( Image* im, Signal *s, int ic1, int ic2, int ic3, 
			double c, char**m ) {
  /***************************************************************************
   * Computes normalized ratio of channels stored in components IC1 and IC2  
   * of the image IM.
   * Returns in signal S the ratio:
   *        S = (IC1/IC2 - c*IC2/IC3) / (IC1/IC2 + c*IC2/IC3)
   ***************************************************************************/
  
  int dimx=im->xdim, dimy=im->ydim, dimv=im->vdim;
  Signal *s1, *s2, *s3;
  
  if(ic1>dimv || ic2>dimv) 
    Error("Index outside the range of possible vectorial components");
  
  if(s == NULL)       TrackNullAlloc( s=create_signal(dimx,dimy) );
  
  /* allocate pointer, do not allocate buffer */
  TrackNullAlloc( s1=create_psignal(dimx,dimy) );
  TrackNullAlloc( s3=create_psignal(dimx,dimy) );
  /* allocate buffer */
  TrackNullAlloc( s2=create_signal(dimx,dimy) );
  
  /* use pointers */
  TrackError( read_im2sig(im,ic1,s1), /* s1 points to channel ic1 */
	      "Error reading input channel" );
  TrackError( read_im2sig(im,ic3,s3), /* s3 points to channel ic3 */
	      "Error reading input channel" ); 

  /* copy of buffer: content of channel ic2 is copied in s2 */
  TrackError( read_im2sig(im,ic2,s2),
	      "Error reading input channel" ); 
  /* copy of buffer: content of channel ic1 (through s1) is copied 
   * in s */
  TrackError( copy_signal(s1,s,m),
	      "Error copying signal" );

  /* Operate channels division */
  TrackError( divide_signal(s2,s,m),    /* s = s/s2 = s1/s2 */
	      "Error operating channel division" );  
  TrackError( divide_signal(s3,s2,m),   /* s2 = s2/s3 */ 
	      "Error operating channel division" );  

  if(c != 1.)
    TrackError( scale_signal(s2,c,m),   /* s2 = c*s2 = c*s2/s3 */ 
		"Error operating scale factor" ); 

  TrackError( substract_signal(s2,s,m), /* s = s - s2 = s1/s2 - c*s2/s3 */ 
	      "Error operating channel difference" ); 
  TrackError( scale_signal(s2,2.,m),    /* s2 = 2*s2 = 2*c*s2/s3 */ 
	      "Error operating scale factor" ); 
  TrackError( add_signal(s,s2,m),   /* s2 = s + s2 
				     *    = s1/s2 - c*s2/s3 + 2*c*s2/s3  */
	      "Error operating channel addition" ); 
						           
  TrackError( divide_signal(s2,s,m),
	      "Error operating channel division" ); /* s = s/s2 */ 

  /* Final result is: (s1/s2 - c*s2/s3) / (s1/s2 + c*s2/s3) */

  /* Make them point on NULL to free the memory */
  s1->signal = s3->signal = NULL;
  free_signal(s1);
  free_signal(s3);
  /* free buffer */
  free_signal(s2);
  
  return OK;
} // end of rvdiffratio_signal


