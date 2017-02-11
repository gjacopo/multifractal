#ifndef      _FLT_H_
#define      _FLT_H_

#ifdef __cplusplus
extern "C" {
#endif	/* __cplusplus */
  
  int laplacian( double **input, int dimx, int dimy, 
		 double **Sxy /*laplacian*/, double **Pxy );
  int canny( double **input, int dimx, int dimy, 
	     double **canny, double **grad );
  int zerocross( double **input, double **grad, double Thresh,
		 int dimx, int dimy, char **zeros);
  
  int filt_gaussian(double **g, double **f, int dimx, int dimy, double sigma);
  int filt_conv2D( double **g,double **f,double *gkr,double *gkc, 
		   int dimx, int dimy, int nk );
  int convR(double **f, double **g, double *gk, int nk, int n, int ir);
  int convC(double **f, double **g, double *gk, int nk, int n, int jc);
  int gkernel(double *gk, double sigma, int n);
  
  int filtbw_3x3isolate( char**mask, int dimx, int dimy );
  int filtbw_3x3erode( char** mask, int dimx, int dimy );

  double matrix_mltscale( double **weight, int wsize, 
			  double exp, int flag_norma );
  double matrix_gauss( double **weight, int wsize, 
		       double sigma, int flag_norma );

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
