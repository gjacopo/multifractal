#ifndef PIX_OPERATOR_H
#define PIX_OPERATOR_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  Pixel * mask2pixel( double **m, int xdim, int ydim, double value,
		      int *count );
  Pixel * cmask2pixel( char **m, int xdim, int ydim, char value,
		       int *count );
  Pixel* file2pixel( FILE *fseed, int *count );
  
  int region2cmask( Region* R, char value, 
		    int dimx, int dimy, char **mask);
  int region2mask( Region* R, double value, 
		   int dimx, int dimy, double **mask);
  int pixel2mask( Pixel * pixlist, double value, 
		  int dimx, int dimy, double **mask );
  int pixel2cmask( Pixel * pixlist, char value, 
		   int dimx, int dimy, char **mask );
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
