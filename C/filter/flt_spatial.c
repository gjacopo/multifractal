#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_alloc.h>      
#include <utl_stats.h>	
#include <utl_operator.h>	

#include <flt_stats1d.h>	
#include <flt_stats2d.h>	
#include <filter.h>


/*****************************************************************/
int laplacian( double **input, int dimx, int dimy, 
	       double **Sxy /*laplacian*/, double **Pxy ) {
  /*****************************************************************/
  /* This routine calculates the Laplacian of vector image.     */
  /*****************************************************************/
  
  int ix, iy, iv;
  double Dx, Dy, Dxx, Dyy;
  double min, max;

  /* vector case
     double *Dx, *Dy, *Dxx, *Dyy;
     if( (Dx=(double*)calloc(dimv,sizeof(double))) == NULL ||
     (Dy=(double*)calloc(dimv,sizeof(double))) == NULL ||
     (Dxx=(double*)calloc(dimv,sizeof(double))) == NULL ||
     (Dyy=(double*)calloc(dimv,sizeof(double))) == NULL )
     return ERROR;
  */
  
  //  normalize( input, dimx, dimy, &min, &max );

  /* Process image */
  for( iy=1; iy<dimy-1; iy++ )
    for( ix=1; ix<dimx-1; ix++ ) {
      
      /* Use cubic splines to estimate derivatives */
      /* vector case: for (iv = 0; iv< dimv; iv++)	 */
      Dx/*[iv]...*/ = (input/*[iv]...*/[iy+1][ix+1] + 2.*input[iy][ix+1]
		    + input[iy-1][ix+1] - input[iy+1][ix-1]
		    - 2.*input[iy][ix-1] - input[iy-1][ix-1]) / 8.;
      
      Dy = (input[iy+1][ix+1] + 2.*input[iy+1][ix]
		    + input[iy+1][ix-1] - input[iy-1][ix+1]
		    - 2.*input[iy-1][ix] - input[iy-1][ix-1]) / 8.;
      
      Dxx = (input[iy+1][ix+1] + 4.*input[iy][ix+1]
		     + input[iy-1][ix+1] - 2.*input[iy+1][ix]
		     - 8.*input[iy][ix] - 2.*input[iy-1][ix]
		     + input[iy+1][ix-1] + 4.*input[iy][ix-1]
		     + input[iy-1][ix-1]) / 6.;
      
      Dyy = (input[iy+1][ix+1] - 2.*input[iy][ix+1]
		     + input[iy-1][ix+1] + 4.*input[iy+1][ix]
		     - 8.*input[iy][ix] + 4.*input[iy-1][ix]
		     + input[iy+1][ix-1] - 2.*input[iy][ix-1]
		     + input[iy-1][ix-1]) / 6.;
      
      /* Calculate directional derivative */
      Pxy[iy][ix] /*vector case:+*/= Dx/*[iv]...*/ * Dx + Dy * Dy; 
      Sxy[iy][ix] /*vector case:+*/= Dxx/*[iv]...*/ + Dyy;
     
      Pxy[iy][ix] = sqrt( Pxy[iy][ix] );
    }
  
  /* Handle boundary rows and columns */
  for( ix=0; ix<dimx; ix++ ) {
    Sxy[0][ix] = Sxy[1][ix];
    Sxy[dimy-1][ix] = Sxy[dimy-2][ix];
  }
  
  for (iy = 0; iy < dimy; iy++) {
    Sxy[iy][0] = Sxy[iy][1];
    Sxy[iy][dimx-1] = Sxy[iy][dimx-2];
  }

  /* vector case:     free memories 
     if( Dx ) free( Dx );
     if( Dy ) free( Dy );
     if( Dxx ) free( Dxx );
     if( Dyy ) free( Dyy );
  */
  return 0;
}

/*****************************************************************/
int canny( double **input, int dimx, int dimy, 
	   double **canny, double **grad ) {
  /*****************************************************************/
  /* Calculate directional derivative for an image.              */
  /*****************************************************************/
  
  int ix, iy, iv;
  double Dx, Dy;
  double **Cx, **Cy; 
  double Cxx, Cxy, Cyx, Cyy;
  
  if( (Cx=matrix2D( dimy, dimx )) == NULL ||  
      (Cy=matrix2D( dimy, dimx )) == NULL )
    return ERROR;
  /* cas vectoriel
     double *Dx, *Dy;
     (Dx=(double*)calloc(dimv,sizeof(double))) == NULL ||
     (Dy=(double*)calloc(dimv,sizeof(double))) == NULL 
  */
  
  /* Process image */
  for( iy=1; iy<dimy-1; iy++ )
    for( ix=1; ix<dimx-1; ix++ ) {
      
      /* Use cubic splines to estimate derivatives */
      /* vector case: for (iv = 0; iv< dimv; iv++)	 { */
      Dx/*[iv]...*/ = (input/*[iv]...*/[iy+1][ix+1] + 2. * input[iy][ix+1]
		       + input[iy-1][ix+1] - input[iy+1][ix-1]
		       - 2. * input[iy][ix-1] - input[iy-1][ix-1]) / 8.;
      
      Dy = (input[iy+1][ix+1] + 2. * input[iy+1][ix]
	    + input[iy+1][ix-1] - input[iy-1][ix+1]
	    - 2. * input[iy-1][ix] - input[iy-1][ix-1]) / 8.;
      
      /* Calculate gradient */
      grad[iy][ix] += Dx/*[iv]...*/ * Dx + Dy * Dy;
      
      /* Calculate first vector derivatives */
      Cx[iy][ix] += Dx * Dx; 
      Cy[iy][ix] += Dy * Dy;
      
      grad[iy][ix] = sqrt( grad[iy][ix] );
      Cx[iy][ix] = sqrt( Cx[iy][ix] ); 
      Cy[iy][ix] = sqrt( Cy[iy][ix] );  
    }
  
  /* Process image */
  for( iy=1; iy<dimy-1; iy++ )
    for( ix=1; ix<dimx-1; ix++ ) {
      
      /* Calculate second vector derivative */
      Cxx = (Cx[iy+1][ix+1] + 2. * Cx[iy][ix+1] + Cx[iy-1][ix+1]
	     - Cx[iy+1][ix-1] - 2. * Cx[iy][ix-1] - Cx[iy-1][ix-1]) / 8.;
      
      Cxy = (Cx[iy+1][ix+1] + 2 * Cx[iy+1][ix] + Cx[iy+1][ix-1]
	     - Cx[iy-1][ix+1] - 2 * Cx[iy-1][ix] - Cx[iy-1][ix-1]) / 8.;
      
      Cyx = (Cy[iy+1][ix+1] + 2 * Cy[iy][ix+1] + Cy[iy-1][ix+1]
	     - Cy[iy+1][ix-1] - 2 * Cy[iy][ix-1] - Cy[iy-1][ix-1]) / 8.;
      
      Cyy = (Cy[iy+1][ix+1] + 2 * Cy[iy+1][ix] + Cy[iy+1][ix-1]
	     - Cy[iy-1][ix+1] - 2 * Cy[iy-1][ix] - Cy[iy-1][ix-1]) / 8.;
      
      /* Calculate directional derivative */
      canny[iy][ix] = (Cx[iy][ix] * Cx[iy][ix] * Cxx + Cx[iy][ix] * Cy[iy][ix] * Cxy
		       + Cx[iy][ix] * Cy[iy][ix] * Cyx + Cy[iy][ix] * Cy[iy][ix] * Cyy)
	/ (0.01 + grad[iy][ix]);
    }
  
  /* Handle boundary rows and columns */
  for( ix=0; ix<dimx; ix++ ) {
    canny[0][ix] = canny[1][ix];
    canny[dimy-1][ix] = canny[dimy-2][ix];
  }
  
  for (iy = 0; iy < dimy; iy++) {
    canny[iy][0] = canny[iy][1];
    canny[iy][dimx-1] = canny[iy][dimx-2];
  }
  
  if( Cx ) free_matrix2D( Cx, dimy );
  if( Cy ) free_matrix2D( Cy, dimy );
  /* vector case:
     if( Dx ) free( Dx );
     if( Dy ) free( Dy );
  */

  return OK;
}

/*****************************************************************/
int zerocross( double **input, double **grad, double Thresh,
	       int dimx, int dimy, char **zeros) {
  /*****************************************************************/
  /* Locate zero crossings in image.                              */
  /*****************************************************************/
  
  int ix, iy, iv;
  double average;
  
  /* Loop through rows looking for zero crossings */
  for( iy=1; iy<dimy-1; iy++ )
    for( ix=1; ix<dimx-1; ix++ )  {
      
      /* Mark pixel location closest to zero */
      if( input[iy][ix-1]<0. && input[iy][ix]>=0. )
	if(-input[iy][ix-1] < input[iy][ix] )
	  zeros[iy][ix-1] = (char)MARKED;
	else
	  zeros[iy][ix] = (char)MARKED;
      
      /* Mark pixel location closest to zero */
      if( input[iy][ix]<0. && input[iy][ix-1]>=0. )
	if (-input[iy][ix] < input[iy][ix-1])
	  zeros[iy][ix] = (char)MARKED;
	else
	  zeros[iy][ix-1] = (char)MARKED;
    }
  
  /* Loop through columns looking for zero crossings */
  for( iy=1; iy<dimy-1; iy++ )
    for( ix=1; ix<dimx-1; ix++ )     { 
      
      /* Mark pixel location closest to zero */
      if(input[iy-1][ix]<0. && input[iy][ix]>=0.)
	if (-input[iy-1][ix] < input[iy][ix])
	  zeros[iy-1][ix] = (char)MARKED;
	else
	  zeros[iy][ix] = (char)MARKED;
      
      /* Mark pixel location closest to zero */
      if(input[iy][ix]<0. && input[iy-1][ix]>=0.)
	if(-input[iy][ix] < input[iy-1][ix])
	  zeros[iy][ix] = (char)MARKED;
	else
	  zeros[iy-1][ix] = (char)MARKED;
    }
  
  /* Remove phantom edges */
  for( iy=1; iy<dimy-1; iy++ )
    for (ix=1; ix <dimx-1; ix++ )
      if (zeros[iy][ix] == (char)MARKED )	{
 	average = (grad[iy+1][ix+1] + grad[iy][ix+1] + grad[iy-1][ix+1]
		   + grad[iy+1][ix] + grad[iy][ix] + grad[iy-1][ix]
		   + grad[iy+1][ix-1] + grad[iy][ix-1] + grad[iy-1][ix-1]) / 9.;
	if (grad[iy][ix]<average || grad[iy][ix]<Thresh)
	  zeros[iy][ix] = (char)NMARKED;
      }
  
  return OK;
}


#define CNT 3

/************************************************************************/
double matrix_gauss( double **weight, int wsize, 
		     double sigma, int flag_norma ) {
  /************************************************************************/
  
  int dx, dy;
  double x, y, norma= 0., fac;
  
  if(sigma < 1.e-30) sigma = 1.;
  fac = 1. / (sqrt(2.*M_PI) * sigma);
  
  for( dy=0; dy<2*wsize+1; dy++ ) {
    y = (double)(dy - wsize);
    for( dx=0; dx<2*wsize+1; dx++ ) {
      x = (double)(dx-wsize);
      weight[dy][dx] = fac * exp( -(x*x+y*y) / (2.*sigma*sigma) );
      norma += weight[dy][dx];
    } 
  }
  
  if(flag_norma)
    /* Normalisation des facteurs de pondération */
    op_scale( 2*wsize+1, 2*wsize+1, 1./norma, weight, NULL); 
  
  return norma;
}

/************************************************************************/
double matrix_mltscale( double **weight, int wsize, 
			double exp, int flag_norma ) {
  /************************************************************************/
  
  int dx, dy;
  double x, y, norma= 0., r;
  
  weight[wsize][wsize] = 0.; /* pixel central */
  
  for( dy=0; dy<2*wsize+1; dy++ ) {
    y = (double)(dy - wsize);
    for( dx=0; dx<2*wsize+1; dx++ ) {
      x = (double)(dx-wsize);
      r = sqrt( x*x + y*y );
      /* Pondération en 1/r^(-sigma) = r^sigma 
       * Remarque : l'argument sigma est negatif */
      if(r > 0.)  
	/* weight[dy][dx] = pow( r, sigma ); */
	weight[dy][dx] = 1. / pow( r, exp );      
      norma += weight[dy][dx];
    }
  }

  if(flag_norma)
    /* Normalisation des facteurs de pondération */
   op_scale( 2*wsize+1, 2*wsize+1, 1./norma, weight, NULL); 
  
  return norma;
}


/***************************************************************************/
int filt_gaussian(double **g, double **f, int dimx, int dimy, double sigma){
  /***************************************************************************
   * Normalized convolution with two 1D Gaussian kernel 
   * Called by mdl
   ***************************************************************************/

  int    nk  = SC2SIZE(sigma);
  double  *gk = centvector(nk, nk+1);
  static double **u, **ug;
  static int x, y;

  //WarnMsg("");
  if (u == NULL){
    TrackNullAlloc( u=matrix2D(dimy,dimx) );
    TrackNullAlloc( ug=matrix2D(dimy,dimx) );
    fill( dimx, dimy, 1.0, u, NULL );
  }
  
  /* 1D Gaussian kernel */
  TrackError( gkernel(gk,sigma,nk), "Error computing gaussian kernel" );

  /* separable kernels */
  TrackError( filt_conv2D(u,ug,gk,gk,dimx,dimy,nk), 
	      "Error processing convolution" );  /*  1*G(sigma)  */
  TrackError( filt_conv2D(g,f,gk,gk,dimx,dimy,nk), 
	      "Error processing convolution" );  /*  I*G(sigma)  */

  /* Normalisation */
   TrackError( op_divide( dimx, dimy, ug, f, NULL), 
	       "Error applying buffer division" ); 
  /* for (y=0; y<dimy; y++)
     for (x=0; x<dimx; x++)
     f[y][x] /= ug[y][x];
  */

  /*  free_matrix2D(u);
      free_matrix2D(ug);
  */

  return OK;
}


/***************************************************************************/
int filt_conv2D( double **g,double **f, double *gkr,double *gkc, 
	    int dimx, int dimy, int nk ){
  /***************************************************************************
   * 2D convolution with separable kernels 
   ***************************************************************************/

  static int i;
  static double **a;
  
  if(a==NULL && (a=matrix2D(dimy, dimx))==NULL)
      return ERROR;

  for( i=0; i<dimy; i++)
    convR( g, a, gkr, nk, dimx, i );
  for( i=0; i<dimx; i++ )
    convC( a, f, gkc, nk, dimy, i );

  return OK;
}


/***************************************************************************/
int convR(double **f, double **g, double *gk, int nk, int n, int ir){
/***************************************************************************/
  static int i, j, j0, j1 ;
  static double s ;

  for( i=0; i<n; i++ ){
    s = 0 ;
    j0 = i - nk; if (j0 < 0) j0 = 0;
    j1 = i + nk; if (j1 > n-1) j1 = n-1;
    for( j=j0; j<j1+1; j++ )
      s += f[ir][j] * gk[i-j];
    g[ir][i] = s ;
  }
  return OK;
}


/***************************************************************************/
int convC(double **f, double **g, double *gk, int nk, int n, int jc){
  /***************************************************************************/
  
  static int i, j, j0, j1 ;
  static double s ;

  for( i=0 ; i<n ; i++ ){
    s = 0 ;
    j0 = i - nk; if (j0 < 0) j0 = 0;
    j1 = i + nk;  if (j1 > n-1) j1 = n-1;
    for( j=j0; j<j1+1 ; j++ )
      s += f[j][jc]*gk[i-j];
    g[i][jc] = s;
  }
  return OK;
}


/***************************************************************************/
int gkernel(double *gk, double sigma, int n){
/***************************************************************************/
  static int i;
  static double s, s2;

  s = 0; s2 = 2.0*sigma*sigma ;
  gk[0] = 1 ;
 
  for( i=1 ; i<n+1 ; i++ )
    s += (gk[i] = exp(-(double)(i*i)/s2)) ;
  s2 = 1.0 + 2.0*s ; /* gk[0] + sum(gk[i]) + sum(gk[-i]) */
  gk[0] /= s2 ;

  for( i=1 ; i<n+1 ; i++ )
    gk[-i] = (gk[i] /= s2) ;

  return OK;
}


/***************************************************************************/
int filtbw_3x3isolate( char** mask, int dimx, int dimy ) {
  /***************************************************************************/
  /* Keep only pixels whose neighborhood is at least of the form: 
   *        1 1 1        1 0 0        0 0 0        0 0 1
   *        0 1 0   or   1 1 0   or   0 1 0   or   0 1 1  
   *        0 0 0        1 0 0        1 1 1        0 0 1
   *              0 0 0        0 1 0 
   *         or   1 1 1   or   0 1 0
   *              0 0 0        0 1 0
   * or with more 1 inside the neighborhood.
   */
  
  int wsize=1;
  int ix, iy, dx, dy;
  int x, y;
  int nneighb[6];
  int **iix, *iiy; /* the neighbours */
  char **tmp=cmatrix2D(dimy,dimx);
  
  TrackNullAlloc( iix=imatrix2D(dimx,2*wsize+1) );
  TrackNullAlloc( iiy=(int*)calloc(2*wsize+1,sizeof(int)) );
  
  for( ix=0; ix<dimx; ix++ ) 
    vecindex( iix[ix]/* vector */, ix/* scalar */, 
	      wsize, dimx, FALSE /* no period */ );     
  
  for( iy=0; iy<dimy; iy++ ) {
    vecindex( iiy, iy, wsize, dimy, FALSE ); 
    
    for( ix=0; ix<dimx; ix++ ) 
      
      if((tmp[iy][ix]=mask[iy][ix]) == TRUE) {
	
	for( dx=0; dx<6; dx++ ) nneighb[dx] = 0;
	
	for( dx=0; dx<=2*wsize; dx++ ) 
	  if((x=iix[ix][dx])>=0) {
	    if((y=iiy[0]) >= 0) nneighb[0] += (int)mask[y][x];
	    if((y=iiy[2]) >= 0) nneighb[2] += (int)mask[y][x];	  
	    nneighb[4] += (int)mask[iiy[1]][x];
	  }
	
	for( dy=0; dy<=2*wsize; dy++ ) 
	  if((y=iiy[dy])>=0) {
	    if((x=iix[ix][0]) >= 0) nneighb[1] += (int)mask[y][x];
	    if((x=iix[ix][2]) >= 0) nneighb[3] += (int)mask[y][x];	  
	    nneighb[5] += (int)mask[y][iix[ix][1]];
	  }
	
	/* TEST */
	if( nneighb[0]!=CNT && nneighb[2]!=CNT &&
	    nneighb[1]!=CNT && nneighb[3]!=CNT &&
	    nneighb[4]!=CNT && nneighb[5]!=CNT ) 
	  tmp[iy][ix] = FALSE;
	
      } /* else tmp[iy][ix]=mask[iy][ix]) already equal to FALSE */
   
  }
  
  TrackError( ccopy(dimx,dimy,tmp,mask,NULL), "Error copying buffer" );
  
  free_cmatrix2D( tmp, dimy );
  free_imatrix2D( iix, dimx );
  Free( iiy );

  return OK;
}



/***************************************************************************/
int filtbw_3x3erode( char** mask, int dimx, int dimy ) {
  /***************************************************************************/
  /* Naive erosion of an image with a 3x3 square structuring element with 
   * 8-connection */
  
  int wsize=1;
  int ix, iy, dd;
  int x, y;
  int nneighb;
  int *ii; /* the neighbours */
  char **tmp=cmatrix2D(dimy,dimx);
  
  TrackNullAlloc( ii=(int*)calloc(2*wsize+1,sizeof(int)) );
  
  for( iy=0; iy<dimy; iy++ ) {
    vecindex( ii, iy, wsize, dimy, FALSE ); 
    for( ix=0; ix<dimx; ix++ ) 
      if((tmp[iy][ix]=mask[iy][ix]) == TRUE) {
	nneighb = 0;
	for( dd=0; dd<=2*wsize; dd++ ) 
	  if((y=ii[dd])>=0) nneighb += (int)mask[y][ix];
	if(nneighb != CNT) tmp[iy][ix] = FALSE;
      }
  }
  
  for( ix=0; ix<dimx; ix++ ) {
    vecindex( ii, ix, wsize, dimx, FALSE ); 
    for( iy=0; iy<dimy; iy++ ) 
      if((mask[iy][ix]=tmp[iy][ix]) == TRUE) {
	nneighb = 0;
	for( dd=0; dd<=2*wsize; dd++ ) 
	  if((x=ii[dd]) >= 0) nneighb += (int)tmp[iy][x];
	if(nneighb != CNT) mask[iy][ix] = FALSE;
      }
  }
  
  free_cmatrix2D( tmp, dimy );
  Free(ii);

  return OK;
}

