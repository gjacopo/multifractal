/* ===================================
** sig_inout.c
** started on Tue Jan 30 08:52:33 2007 
** ===================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <sig_inout.h>



/***************************************************************************/
int read_im2sig( Image* im, int ic, Signal* sig ) {
/***************************************************************************/

  if(sig->signal == NULL)          sig->signal = &(im->image[ic][0]);
  else      copy(sig->xdim,sig->ydim,im->image[ic],sig->signal,NULL);

  return OK;
} // end of read_im2sig


/***************************************************************************/
int read_im2buf( Image* im, int ic, double** s ) {
/***************************************************************************/

  if(s == NULL)          s = &(im->image[ic][0]);
  else     copy(im->xdim,im->ydim,im->image[ic],s,NULL);

  return OK;
} // end of read_im2buf

/***************************************************************************/
int write_sig2im( Signal* sig, int ic, Image* im ) {
/***************************************************************************/

  if(im->image[ic] == NULL)   im->image[ic] = sig->signal;
  else     copy(sig->xdim,sig->ydim,sig->signal,im->image[ic],NULL);
  
  return OK;
} // end of write_sig2im


/***************************************************************************/
int write_buf2im( double** s, int dimx, int dimy, int ic, Image* im ) {
/***************************************************************************/

  if(im->image[ic] == NULL)   im->image[ic] = s;
  else     copy(dimx,dimy,s,im->image[ic],NULL);
  
  return OK;
} // end of write_buf2im


/***************************************************************************/
int dimension_image( char *name, int flagfoto, int flagcolor, 
		     int *bd, int *type,
		     int *xdim, int *ydim, int *vdim, int *zdim ) {
  /***************************************************************************/

  if(flagfoto == TRUE) {/* Picture image */
    *type = 0; 
    *bd = 1; 
    TrackError( read_dimensiones_foto(name,xdim,ydim),
		"Error reading PIC image header" );
    *zdim = 1;
    if(flagcolor == TRUE) *vdim = 3;
    else *vdim = 1;
  } else if(flagfoto == -TRUE)	{/* Inrimage image */
#ifdef  FLAG_INRIMAGE
    TrackError( read_hdrinr(name,bd,type,xdim,ydim,vdim,zdim),
		"Error reading INR image header" );
#else
    *type = 0; 
    TrackError( autoread_hdrinr(name,bd,xdim,ydim,vdim,zdim),
		"Error reading INR image header" );
#endif
  } else {
    TrackError( read_dimensiones_raw(name,bd,xdim,ydim,vdim,zdim),
		"Error reading RAW image header" );
    *type = 0; 
  }
  
  return OK;
} // end of dimension_image


/***************************************************************************/
Image *init_image( char *name, int flagfoto, int flagcolor ) {
  /***************************************************************************/
  
  Image *im;
  int xdim, ydim, vdim, N;
  int type, bd;
  int ic, iz;

  if(dimension_image( name, flagfoto, flagcolor, &bd, &type, 
		      &xdim, &ydim, &vdim, &N ) == ERROR)
    return NULL;
  
  if((im=create_image(xdim,ydim,vdim,N)) == NULL)    return NULL;

  // fillname_image(im, name);
  im->name = name;
  
  im->foto=flagfoto;
  im->color=flagcolor;
  im->bd = bd;
  im->type = type;
  
  for( iz=0; iz<N; iz++ ) 
    im->flag_read[iz] = (char)FALSE; /* nothing read yet */
  
  return im;
} // end of init_image


/***************************************************************************/
int read_image( Image *im, int iz ) {
  /***************************************************************************/
  
  int levels;
  int dimx=im->xdim, dimy=im->ydim, dimv=im->vdim;
  int N=im->zdim;
  double* med_cr;

  TrackNullAlloc( im->image );
  TrackNullAlloc( med_cr=(double*)calloc(dimv,sizeof(double)) );

  if(im->foto == TRUE) { /* PIC image */
    if(dimv == 3)
      levels = read_color_block( dimx, dimy, 1., im->name, 
				 im->image, med_cr ); 
    else if(dimv == 1)
      levels = read_foto_gris( dimx, dimy, im->name, im->image[0]);
    
  } else if(im->foto == -TRUE) { /* INR image */
#ifdef  FLAG_INRIMAGE
    read_inr( im->name, im->bd, dimx, dimy, dimv, N, iz, 
	      im->image );
#else 
    autoread_inr( im->name, im->bd, dimx, dimy, dimv, N, iz,
		  im->image );
#endif
  } else  /* RAW image */
    read_raw( im->name, im->bd, dimx, dimy, dimv, N, iz, 
	      im->image );
  
  im->iz = iz;
  im->flag_read[iz] = (char)TRUE; /* the frame has been loaded */

  Free(med_cr);

  return OK;
} // end of read_image

/***************************************************************************/
int write_image( char *base, int flagvideo, Image *im, int iz ) {
/***************************************************************************/
  
#ifdef FLAG_INRIMAGE
  if(im->foto == TRUE) {
#endif
    TrackError( write_foto(im->xdim, im->ydim, im->vdim, im->zdim,
			   /*im->iz*/ iz, base, im->name, im->image ),
		"Error writing final result");
    
#ifdef FLAG_INRIMAGE
  } else { /* Inrimage image */
    TrackError( write_inr(im->bd, im->xdim, im->ydim, im->vdim, im->zdim, 
			  /*im->iz*/ iz, im->name, im->image),
		"Error writing final result" );
  }
#endif
  
  return OK;
} // end of write_image

