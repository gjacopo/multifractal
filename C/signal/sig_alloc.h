#ifndef   	_SIG_ALLOC_H_
#define   	_SIG_ALLOC_H_
 
#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  Image* alloc_image( );
  Image * create_pimage( int dimx, int dimy, int dimv, int dimz );
  int free_pimage( Image *image );
  int fillname_image( Image* im, char *name );
  Image* create_image( int dimx, int dimy, int dimv, int dimz );
  int free_image( Image *image );
  int dimension_image( char *name, int flagfoto, int flagcolor, 
		       int *bd, int *type,
		       int *xdim, int *ydim, int *vdim, int *zdim );
  Image * init_image( char *name, int flagfoto, int flagcolor );
  int read_image( Image *image, int iz );
  int nullpimage( Image* im );
  
  int write_image( char *base, int flagvideo, Image *im, int iz );
  
  int  display_imageinfo(Image *im, char *title);
  
  Signal* alloc_signal(  );
  Signal* create_psignal( int dimx, int dimy );
  int free_psignal( Signal* sig );
  Signal* create_signal( int dimx, int dimy );
  int free_signal( Signal* sig );
  Signal* create_nullsignal( int dimx, int dimy );
  
  int read_im2sig( Image* im, int ic, Signal* sig );
  int read_im2buf( Image* im, int ic, double** s );
  int write_sig2im( Signal* sig, int ic, Image* im );
  int write_buf2im( double** s, int dimx, int dimy, int ic, Image* im );
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !SIG_ALLOC_H_ */
