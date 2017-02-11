#include <stdio.h>
#include <math.h>
	 	
#ifdef  FLAG_INRIMAGE
#include <inrimage/image.h>
#endif

/* Librairies  ROUTINES */

#include <utils.h>
#include <utl_alloc.h>	 	
#include <utl_char.h>	 	
#include <utl_stats.h>	 	
#include <utl_operator.h>

#include <inout.h>              
#include <io_inria.h>              
#include <io_grafic.h>

#include <flt_stats1d.h>
#include <flt_stats2d.h>

#include <signal.h>
#include <sig_alloc.h>
#include <sig_gradient.h>


/***************************************************************************/
Image *alloc_image( ) {
/***************************************************************************/
  Image *im;

  if( (im=(Image*)malloc(sizeof(Image))) == NULL)
    return NULL;
  else return im;
}  // end of alloc_image


/***************************************************************************/
Image * create_pimage( int itype, int ix, int iy, int iv, 
		       int iz ) {
/***************************************************************************/
  
  Image *im;
  int ic;
  
  if((im=alloc_image()) == NULL)        return NULL;

  /* Information about the dimension of the image itself */
  im->xdim=ix;
  im->ydim=iy;
  im->vdim=iv;
  /* Information about the dimension of the sequence */
  im->zdim=iz;

  CASETYPEDO(itype, 
	     im->image=(MYCHAR***)calloc(iv,sizeof(MYCHAR**)),
	     im->image=(char***)calloc(iv,sizeof(char**)),
	     im->image=(uchar***)calloc(iv,sizeof(uchar**)),
	     im->image=(***)calloc(iv,sizeof(**)),
	     im->image=(***)calloc(iv,sizeof(**)),
	     im->image=(***)calloc(iv,sizeof(**)),
	     im->image=(double***)calloc(iv,sizeof(double**)))
  
  if( (im->image == NULL) ||
      (im->flag_read=(char*)malloc(iz*sizeof(char))) == NULL)
    return NULL;
  
  for( ic=0; ic<iz; ic++) im->image[ic] = NULL;
  
  im->prec = im->suiv = NULL;

  return im;
} // end of create_pimage


/***************************************************************************/
Image * create_image( int itype, int ix, int iy, int iv, 
		      int iz ) {
/***************************************************************************/
  
  Image *im;
  int ic;
  
  if( (im=create_pimage(itype,ix,iy,iv,iz)) == NULL) return NULL;
  
  for( ic=0; ic<iz; ic++) 
    if( (im->image[ic] = matrix2D(iy,ix)) == NULL)    
      return NULL;

  return im;
} // end of create_image


/***************************************************************************/
int fillname_image( Image* im, char *name ) {
  /***************************************************************************/
  
  int l;
  /*  for( l=0; l<MAXNAMELENGTH; l++ ) im->name[l]='\0';
  TrackNullAlloc( strcpy(im->name,name) );
  */
  sprintf(im->name,"%s",name);
  return OK;
} // end of fillname_image



/***************************************************************************/
int free_pimage( Image *im ) {
/***************************************************************************/

  Free(im->image);
  Free(im->flag_read);
  Free(im);

  return OK;
} // end of free_pimage


/***************************************************************************/
int free_image( Image *im ) {
/***************************************************************************/

  if(im->image != NULL)
    free_matrix3D(im->image,im->vdim,im->ydim);

  Free(im->flag_read);
  Free(im);

  return OK;
} // end of free_image


/***************************************************************************/
int nullpimage( Image* im ) {
/***************************************************************************/
  int ic;

  TrackNullAlloc( im );
  TrackNullAlloc( im->image );

  for( ic=0; ic<im->vdim; ic++ )    im->image[ic] = NULL;

  return OK;
} // end of nullpimage


/***************************************************************************/
int  display_imageinfo(Image *im, char *title){
  /***************************************************************************/
  int l; ;
  char strbar[MAXTEXTLENGTH];
  
  l=fprintf(stderr,"\n %s info - %s : ", title, im->name);
  repeat_char(strbar,' ',l-1,TRUE);

  switch(im->foto) {
  case -TRUE:
    fprintf(stderr,"Type : INR"); 
    break;
  case TRUE:
    fprintf(stderr,"Type : PIC"); 
    break;
  case FALSE:
    fprintf(stderr,"Type : RAW"); 
  }
  if(im->foto != TRUE) 
    fprintf(stderr,"\n%sByte : %d",strbar, im->bd);
  fprintf(stderr,"\n%sX : %d",strbar, im->xdim);
  fprintf(stderr,"\n%sY : %d",strbar, im->ydim);
  fprintf(stderr,"\n%sV : %d",strbar, im->vdim);
  if(im->zdim > 1) 
    fprintf(stderr,"\n%sZ : %d",strbar, im->zdim);
  
  return OK;
} // end of display_imageinfo


/***************************************************************************/
Signal* alloc_signal(  ) {
  /***************************************************************************/
  Signal *sig; 
  
  if( (sig=(Signal*)malloc(sizeof(Signal))) == NULL)
    return NULL;
  else 
    return sig;
} // end of alloc_signal


/***************************************************************************/
Signal* create_psignal( int dimx, int dimy ) {
  /***************************************************************************/
  Signal* sig;

  if((sig=alloc_signal()) == NULL) return NULL;
  
  sig->xdim = dimx, sig->ydim = dimy;
  
  sig->mms = NULL;
  /* 
     if((sig->mms=(double*)malloc(3*sizeof(double))) == NULL ) {
     free_signal(sig);
     return NULL;
     }
  */
  sig->signal = NULL;
  sig->hist = NULL;
  sig->grad = NULL;
  sig->mask = NULL;

  return sig;
} // end of create_psignal


/***************************************************************************/
Signal* create_signal( int dimx, int dimy, int *flag ) {
  /***************************************************************************/
  Signal* sig;
  int error=OK;
  
  /* memory allocation */
  if((sig=alloc_signal()) == NULL)      return NULL;

  /* initialisation of the signal */
  sig->xdim = dimx, sig->ydim = dimy;  
  if( (sig->signal=matrix2D(dimy,dimx)) == NULL ||
      (sig->mms=(double*)malloc(3*sizeof(double))) == NULL )
    error = ERROR;
  else
    sig->vsignal = &(sig->signal[0][0]);
  
  /* default initialisations of other fields */
  sig->hist = NULL, sig->grad = NULL, sig->mask = NULL;
  
  if(flag != NULL) {
    sig->flag = flag;
    if(flag[IFLAG_MASK] && ((sig->mask=create_mask(dimy,dimx))==NULL)) 
      error = ERROR;
    if(flag[IFLAG_GRAD] && ((sig->grad=create_gradient(sig))==NULL)) 
      error = ERROR;
    if(flag[IFLAG_HIST] && ((sig->hist=create_histo(NBOX))==NULL)) 
      error = ERROR;
  }
  
  /*
  if( (sig=create_psignal(dimx,dimy)) == NULL)    return NULL;
  if( (sig->signal=matrix2D(dimy,dimx)) == NULL)    return NULL;
  */

  if(error == ERROR) {
    free_signal(sig);
    return NULL;
  }
  else return sig;
} // end of create_signal


/***************************************************************************/
int init_signal(Signal* sig, Signal def) {
  /***************************************************************************/
  
  TrackNullAlloc( sig=create_signal(def.dimx,def.dimy,def.flag) );
  TrackError( copy_signal(def,sig), "Error copying the default signal");
  
  return OK;
}


/***************************************************************************/
Signal* create_nullsignal( int dimx, int dimy ) {
  /***************************************************************************/
  Signal* szero;
  
  if(((szero=create_signal(dimx,dimy))==NULL) ||
     (fill_signal(szero,0.,NULL)==ERROR))
    return NULL;
  
  else
    return szero;
} // end of create_nullsignal


/***************************************************************************/
int free_psignal( Signal* sig ) {
/***************************************************************************/

  sig->signal = NULL;
  Free(sig);

  return OK;
} // end of free_psignal


/***************************************************************************/
int free_signal( Signal* sig ) {
/***************************************************************************/

  if((sig->signal==NULL) || (sig->dimx<1) || (sig->dimy<1)) return ERROR;

  if(sig->signal != NULL)   free_matrix2D(sig->signal,sig->ydim);
  if (sig->vsignal != NULL) free(sig->vsignal);
  
  if(sig->hist != NULL) Free(sig->hist);
  if(sig->grad != NULL) free_gradient(sig->grad);
  if(sig->mask != NULL) free_mask(sig->mask);
  if(sig->mms != NULL)  Free(sig->mms);

  Free(sig);
  return OK;
} // end of free_signal


/***************************************************************************/
int copy_signal( Signal *sin, Signal* sdest, char **m ) {
  /***************************************************************************/
  int error=OK;
  
  if(sdest == NULL)
    TrackNullAlloc( sdest=create_signal(sin->xdim,sin->ydim,sin->flag) );
  
  TrackError(compare_dimsignal(sin,sdest),
	     "Signals dimensions do not match");

  error = copy(sin->xdim,sin->ydim,sin->signal,sdest->signal,m);
  
  if(sin.flag[IFLAG_GRAD] && (sin->grad!=NULL))
    error = copy_gradient(sin->grad,sdest->grad,m);
  
  if(sin.flag[IFLAG_MASK] && (sin->mask!=NULL))
    error = ccopy(sin->xdim,sin->ydim,sin->mask,sdest->mask,m);
  
  if(sin.flag[IFLAG_HIST] && (sin->histo!=NULL)) 
    error = copy_histo(sin->histo,sdest->histo);
  
  return error;
}

/***************************************************************************/
Mask* alloc_mask(  ) {
  /***************************************************************************/
   Mask *mask; 
  
  if( (mask=(Mask*)malloc(sizeof(Mask))) == NULL)
    return NULL;
  else 
    return mask;
} // end of alloc_mask


/***************************************************************************/
Mask* create_mask( int dimx, int dimy ) {
  /***************************************************************************/
  Mask* mask;
  int iy;
  int error=OK;
  
  /* memory allocation */
  if((mask=alloc_mask()) == NULL)      return NULL;
  
  /* initialisation of the signal */
  mask->xdim = dimx, mask->ydim = dimy;  
  if((mask->signal=cmatrix2D(dimy,dimx)) == NULL) {
    free_mask(mask);
    return NULL;
    
  } else /* the start of the horizontal rows is redirected to the 
	  * corresponding array pointer */
      mask->vsignal = &(mask->signal[0][0]);
  
  return mask;  
} // end of create_mask

/***************************************************************************/
int free_mask( Signal*  ) {
/***************************************************************************/

  if((mask->mask==NULL) || (mask->dimx<1) || (mask->dimy<1)) return ERROR;

  if(mask->mask != NULL)  free_cmatrix2D(mask->mask,mask->ydim);
  if(mask->vmask != NULL) free(mask->vmask);
  
  Free(mask);
  return OK;
} // end of free_mask


/***************************************************************************/
int copy_mask( Signal *min, Signal* mdest, char **m ) {
  /***************************************************************************/
  
  if(mdest == NULL)
    TrackNullAlloc( mdest=create_mask(min->xdim,min->ydim) );
  
  TrackError(compare_dimmask(min,mdest),
	     "Masks dimensions do not match");

  ccopy(min->xdim,min->ydim,min->mask,sdest->mask,m);

  return OK;
} // end of copy_mask
