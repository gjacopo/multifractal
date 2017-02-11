/* ===================================
** mf_source.c
** ===================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*	My libraries 		*/
#include <utils.h>
#include <utl_alloc.h>
#include <utl_stats.h>
#include <utl_op.h>	
	
#include <flt_.h>
#include <flt_stats.h>
#include <flt_fft2d.h>		
#include <flt_deriva2d.h>		
#include <flt_wavelet2d.h>


#include <mfractal.h>
#include <mf_source.h>

/* Global external variables */

#ifdef _PARSE_FRACTAL_PARAMETERS_
#include <mf_parse.h>
extern ParMSM *p_msm;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */

#ifndef OPTIMAL_PATCHEA
#define OPTIMAL_PATCHEA FALSE
#endif

/***************************************************************************
 * reconstruct2D_FFT called by:    deriva_orient_limits
 *               compute_UPM
 *               recupera_dummy
 *               recupera_gradiente_esencial
 *               determine_unitary 
 *               kernelea
 *               deskernelea
 ***************************************************************************/


/***************************************************************************/
int determine_MSM( int dimx, int dimy, int dimv, int dimz, int iz, 
		   int xeff, int yeff, int levels, 
		   double ***signal, char ***msm) {
  /***************************************************************************
   * Computes the MSM of a greyscale image.
   ***************************************************************************/
  double **gx,**gy;
  int ic;
  
  TrackNullAlloc( gx=matrix2D(yeff,xeff) );
  TrackNullAlloc( gy=matrix2D(yeff,xeff) );
  
  for(ic=0;ic<dimv;ic++) {

    /* fill0(xeff,yeff,gx,NULL); */
    copy(dimx,dimy,signal[ic],gx,NULL);

    /* Calcul du signal derive */
    deriva2D_bak(xeff,yeff,gx,gy);

    /* Initialisation de la MSM */
    cfill(dimx,dimy,C0,msm[ic],NULL);
    
    compute_upm( dimx,dimy,xeff,yeff,levels,signal[ic],gx,gy,msm[ic] );    
    
    singularize_sign(dimx,dimy,xeff,yeff,0,msm[ic],gx,gy);
  }
  
  /* Computing/displaying the multifractal statistics if needed	
     IF(MULTIFRACTAL) 
     write_multifractal(dimx,dimy,dimv,dimz,iz, ext,msm,signal);
     write_foto_RGB(VIDEO,BLOUT,dimx,dimy,dimv,dimz,iz,"contour",ext,msm,msm,msm);
  */
  
  free_matrix2D(gx,yeff);
  free_matrix2D(gy,yeff);
  
  return OK;
} // end of determine_msm


/***************************************************************************/
int determine_unitary( int dimx, int dimy, int dimv,  int xeff, int yeff,
		       double ***gx, double ***gy, double ***dummy ) {
  /***************************************************************************
   * Computes the chromatically reduced image reconstructed using an unitary 
   * gradient instead of the original one (but with the same direction and 
   * orientation).
   ***************************************************************************/
  double **ax,**ay;
  double mod;
  int ix,iy,ic;
  
  TrackNullAlloc( ax=matrix2D(yeff,xeff) );
  TrackNullAlloc( ay=matrix2D(yeff,xeff) );
  
  for(ic=0;ic<dimv;ic++)	{
    for(iy=0;iy<dimy;iy++) {
      for(ix=0;ix<dimx;ix++)		{
	/* Calcul du  module du gradient */
	mod = sqrt(gx[ic][iy][ix]*gx[ic][iy][ix]
		   + gy[ic][iy][ix]*gy[ic][iy][ix]);
	/* Gradient normalise (unitaire) mais dans la meme direction 
	 * que le gradient de l'image */
	if(mod > 1.e-30)  {
	  ax[iy][ix] = gx[ic][iy][ix] / mod;
	  ay[iy][ix] = gy[ic][iy][ix] / mod;
	} else		  
	  ax[iy][ix] = ay[iy][ix] = 0.;
	
      }
    }
    /* Reconstruction a partir du gradient unitaire d'une image
     * chromatiquement reduite */
    reconstruct2D_FFT(xeff,yeff,ax,ay);
    copy(dimx,dimy,ax,dummy[ic],NULL);
  }
  
  free_matrix2D(ax,yeff);
  free_matrix2D(ay,yeff);

  return OK;
} // end of determine_unitary


/***************************************************************************/
double compute_upm( int dimx, int dimy, int xeff, int yeff, int levels, 
		    double **signal, 
		    double **gx, double **gy, char **msm) {
  /***************************************************************************/
  /* Called by determine_msm */
  /***************************************************************************/
  double **disp,**errorx,**errory;
  double umbral,dens;
  double norsignal,norerror;
  double psnr;

  double upm_dens, upm_thr;
   
  TrackNullAlloc( disp=matrix2D(dimy,dimx) );
  
  TrackNullAlloc( norsignal=dispersion(dimx,dimy,signal,NULL) );
  TrackNullAlloc( errorx=matrix2D(yeff,xeff) );
  TrackNullAlloc( errory=matrix2D(yeff,xeff) ); 
  
  /*** Operations pour la determination de la MSM */
#ifdef _PARSE_FRACTAL_PARAMETERS_
  switch(p_msm->mode_upm)    
#else
    switch(MODE_UPM)  
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
      {
      case MODUPM:
	genera_upm(dimx,dimy,signal,disp);
	break;
      case MODCORRUPM:
	genera_upm_corr(dimx,dimy,signal,gx,gy,disp); 
	break;
      case MODGLOBUPM:
	genera_upm_glob(dimx,dimy,signal,gx,gy,disp);	
	break;
      case MODRELUPM:
      default:
	genera_upm_rel(dimx,dimy,signal,gx,gy,disp);
      }
  
  /* Calcul du seuil de determination des valeurs des exposants de la MSM 
   * a partir de leur distribution */
#ifdef _PARSE_FRACTAL_PARAMETERS_
  upm_dens = p_msm->upm_dens;
  upm_thr = p_msm->upm_thr;
#else
  upm_dens = UPM_DENS;
  upm_thr = UPM_THR;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
  if(upm_dens > 0.) 
    umbral=cuantil_2D(dimx,dimy,disp,upm_dens) / media(dimx,dimy,disp,NULL);
  else 
    umbral=upm_thr;
  
  /* Determination des points de la MSM ou du complementaire */
  dens=mascara_comp(dimx,dimy,
#ifdef _PARSE_FRACTAL_PARAMETERS_
		    ((p_msm->flag_msmcomp) ? -1 : 1),
#else
		    ((FLAG_MSMCOMP) ? - 1 : 1),
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
		    umbral,disp,msm);
  
  /*** Operations pour le calcul de la PSNR */
  
  copy(xeff,yeff,gx,errorx,NULL);
  copy(xeff,yeff,gy,errory,NULL);
  /* errorx et erroy contiennent les gradients de l'image */
  
  /* Mise a 0 des valeurs du gradient en dehors de la MSM */
  singularize(dimx,dimy,xeff,yeff,msm,errorx,errory); 
  /* Reconstruction a partir du gradient sur la MSM
   * errorx contient l'image reconstruite */
  reconstruct2D_FFT(xeff,yeff,errorx,errory);
  /* Represent reconstruction: */
  
  IFDEBUG  { 
    WarnMsg("create kernel.gif...");
    write_foto_block(dimx,dimy,BLOCKOUT,"kernel.gif",errorx); 
  }
  
  /* Signal-to-noise (SNR) measures estimates the quality of the reconstructed 
   * image compared with the original image.
   * Computation of the mean squared error MSE of the reconstructed image IR 
   * as follows:
   *       MSE = \frac{\sum_x (IR(x) - I(x))^2} {N^2}
   * where I is the original image and N is the size of the image
   */
  /* 1) Compute the difference between original and reconstructed images */
  op_diff(dimx,dimy,signal,errorx,NULL);
  /* 2) Compute the root mean squared error RMSE which is the square root of MSE:
   *       RMSE = \sqrt{MSE}                                               */ 
  norerror = dispersion(dimx,dimy,errorx,NULL) / 
    ((double)levels)*256.;
  /* 3) Compute the PSNR as follows:
   *       PSNR = 20 \log_10 (\frac{255.}{RMSE})
   * Note that error metrics are computed on the luminance signal only so 
   * the pixel values I and IR range between black (0) and white (255).
   */
  psnr = 20. * log(255./norerror) / log(10.);
  
  IFVERBOSE 
    WarningVV("Density of UPM : %f%% at PSNR %f dB\n",
#ifdef _PARSE_FRACTAL_PARAMETERS_
	      ((p_msm->flag_msmcomp) ? (100.*(1.-dens)): (100.*dens)),
#else
	      ((FLAG_MSMCOMP) ?  (100.*(1.-dens)): (100.*dens)),
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
	      psnr ); 
  
  /*	Freeing memory associated to the routine	*/
  free_matrix2D(errorx,yeff);
  free_matrix2D(errory,yeff);
  free_matrix2D(disp,dimy);
  
  return dens;
} // end of compute_upm

  
/****************************************************************************/
int genera_upm( int dimx, int dimy, double **signal, double **disp) {
  /****************************************************************************/
  double ***Fourier;
  double *parche;
  int ix,iy;
  
  TrackNullAlloc( parche=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Fourier=matrix3D(4,5,5) );
  
  init_cross(Fourier);
  
  for(iy=0;iy<dimy;iy++) 
    for(ix=0;ix<dimx;ix++) {
      recorta(dimx,dimy,ix,iy,signal,parche);
      disp[iy][ix] = patchea(parche,Fourier);
    }
  
  free(parche);
  free_matrix3D(Fourier,4,5);

  return OK;
} // end of genera_upm


/****************************************************************************/
int genera_upm_corr( int dimx, int dimy, double **signal, 
		     double **gx, double **gy, double **disp) {
  /****************************************************************************/
  double ***Fourier;
  double **errx,**erry;
  double *parche;
  double sigproy;
  int *iix,*iiy;
  int ix,iy,dx,dy,width=1;

  TrackNullAlloc( errx=matrix2D(dimy,dimx) );
  TrackNullAlloc( erry=matrix2D(dimy,dimx) );
  TrackNullAlloc( parche=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( iix=(int*)calloc(2*width+1,sizeof(int)) );
  TrackNullAlloc( iiy=(int*)calloc(2*width+1,sizeof(int)) );
  
  TrackNullAlloc( Fourier=matrix3D(4,5,5) );
  init_cross(Fourier);

  for(iy=0;iy<dimy;iy++) 
    for(ix=0;ix<dimx;ix++) {
      recorta(dimx,dimy,ix,iy,signal,parche);
      patchea_vec(parche,Fourier,&errx[iy][ix],&erry[iy][ix]);
    }
  
  for(iy=0;iy<dimy;iy++)    {
    for(dy=0;dy<2*width+1;dy++) iiy[dy]=Mod(iy+dy-width,dimy);
    for(ix=0;ix<dimx;ix++)      {
      for(dx=0;dx<2*width+1;dx++)
	iix[dx]=Mod(ix+dx-width,dimx);
      
      disp[iy][ix]= 0.;
      for(dy=0;dy<2*width+1;dy++)	   
	for(dx=0;dx<2*width+1;dx++)	  {
	  sigproy=gx[iy][ix]*gx[iiy[dy]][iix[dx]]+
	    gy[iy][ix]*gy[iiy[dy]][iix[dx]];
	  if(sigproy>=0) sigproy=1.;
	  else sigproy=-1.;
	  
	  disp[iy][ix]+=
	    errx[iy][ix]*errx[iiy[dy]][iix[dx]]+
	    erry[iy][ix]*erry[iiy[dy]][iix[dx]];
	}
      disp[iy][ix]=sqrt(fabs(disp[iy][ix])); 
    }
  }
  
  free_matrix2D(errx,dimy);
  free_matrix2D(erry,dimy);
  free(parche);
  free(iix);  free(iiy);
  free_matrix3D(Fourier,4,5);

  return OK;
} // end of genera_upm_corr


/***************************************************************************/
int genera_upm_rel( int dimx, int dimy, double **signal, 
			    double **gx, double **gy, double **disp) {
  /***************************************************************************/
  /* Computes the values and the distribution of singularity exponents
   * Returns the matrice disp of singularity exponents for exch pixel of the image.
   * ========================================================================== 
   * Called by compute_upm */
  /***************************************************************************/
  double ***Fourier;
  double **errx,**erry;
  double *parche;
  double proy,mod;
  int *iix,*iiy;
  int ix,iy,dx,dy,width=1;
  
  TrackNullAlloc( Fourier=matrix3D(4,5,5) );
  init_cross(Fourier);
  /* Matrix for Cross Fourier Transform: Fourier[0], Fourier[1]
   * Inversa matrix for Cross Fourier Transform: Fourier[2], Fourier[3]  */
  
  TrackNullAlloc( errx=matrix2D(dimy,dimx) );
  TrackNullAlloc( erry=matrix2D(dimy,dimx) );
  TrackNullAlloc( parche=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( iix=(int*)calloc(2*width+1,sizeof(int)) ) ;
  TrackNullAlloc( iiy=(int*)calloc(2*width+1,sizeof(int)) );
  
  for(iy=0;iy<dimy;iy++)  
    for(ix=0;ix<dimx;ix++)	{
      /* Extraction d'un patch (voisinage en croix) autour du pixel d'interet 
       * du signal */
      recorta(dimx,dimy,ix,iy,signal,parche);
      /* Convolution with the micro-wavelet */
      IF(OPTIMAL_PATCHEA) 
	patchea_vec_optimal(parche,Fourier,&errx[iy][ix],&erry[iy][ix]);
      ELSE
	patchea_vec(parche,Fourier,&errx[iy][ix],&erry[iy][ix]);
    }
  
  for(iy=0;iy<dimy;iy++) {
    for(dy=0;dy<2*width+1;dy++) iiy[dy]=Mod(iy+dy-width,dimy);
    for(ix=0;ix<dimx;ix++)	{
      for(dx=0;dx<2*width+1;dx++)
	iix[dx]=Mod(ix+dx-width,dimx);
      
      mod = proy = 0.;
      for(dy=0;dy<2*width+1;dy++)	{
	for(dx=0;dx<2*width+1;dx++)	  {
	  /* Somme des modules du gradient de l'image dans le voisinage considere */
	  proy +=
	    gx[iiy[dy]][iix[dx]]* gx[iiy[dy]][iix[dx]]
	    + gy[iiy[dy]][iix[dx]]* gy[iiy[dy]][iix[dx]];
	  
	  /* Somme des modules de err dans le meme voisinage */ 
	  mod +=
	    errx[iy][ix]*errx[iiy[dy]][iix[dx]] 
	    + erry[iy][ix]*erry[iiy[dy]][iix[dx]];
	}
      }

      if(proy > 1.e-30) disp[iy][ix] = sqrt(fabs(mod/proy));
      else disp[iy][ix]=1.;
    }
  }

  free_matrix2D(errx,dimy);
  free_matrix2D(erry,dimy);
  free(parche);
  free(iix);  free(iiy);
  free_matrix3D(Fourier,4,5);

  return OK;
} // end of genera_upm_rel


/***************************************************************************/
int genera_upm_glob( int dimx, int dimy, double **signal, 
			     double **gx, double **gy, double **disp) {
  /***************************************************************************/
  
  double **errx,**erry;
  double mod,proy; 
  int *iix,*iiy;          
  int ix,iy;
  int dx,dy,width=1;
  
  TrackNullAlloc( iix=(int*)calloc(2*width+1,sizeof(int)) );
  TrackNullAlloc( iiy=(int*)calloc(2*width+1,sizeof(int)) );
  TrackNullAlloc( errx=matrix2D(dimy,dimx) );
  TrackNullAlloc( erry=matrix2D(dimy,dimx) );
  
  for(iy=0;iy<dimy;iy++) {
    for(dy=0;dy<2*width+1;dy++) iiy[dy] = Mod(iy+dy-width,dimy);
    for(ix=0;ix<dimx;ix++)      {
      for(dx=0;dx<2*width+1;dx++)
	iix[dx] = Mod(ix+dx-width,dimx);    
      errx[iy][ix] = 0.25*(gx[iiy[2]][ix]+gx[iiy[0]][ix]+
			   gx[iy][iix[2]]+gx[iy][iix[0]]);
      erry[iy][ix] = 0.25*(gy[iiy[2]][ix]+gy[iiy[0]][ix]+
			   gy[iy][iix[2]]+gy[iy][iix[0]]);
    }
  }
  
  for(iy=0;iy<dimy;iy++)    {
    for(dy=0;dy<2*width+1;dy++) iiy[dy]=Mod(iy+dy-width,dimy);
    for(ix=0;ix<dimx;ix++)     {
      for(dx=0;dx<2*width+1;dx++)
	iix[dx]=Mod(ix+dx-width,dimx);
      mod = proy = 0.;
      for(dy=0;dy<2*width+1;dy++)	{
	for(dx=0;dx<2*width+1;dx++)	  {
	  proy += gx[iiy[dy]][iix[dx]]* gx[iiy[dy]][iix[dx]]
	    + gy[iiy[dy]][iix[dx]]* gy[iiy[dy]][iix[dx]];
	  
	  mod += errx[iy][ix]*errx[iiy[dy]][iix[dx]]
	    + erry[iy][ix]*erry[iiy[dy]][iix[dx]];
	}
      }
      if(proy > 1.e-30) disp[iy][ix]=sqrt(fabs(mod/proy));
      else disp[iy][ix]=1.;
      
    }       
  }   
            
  free_matrix2D(errx,dimy);
  free_matrix2D(erry,dimy);

  free(iix);
  free(iiy);  

  return OK;
} // end of genera_upm_glob


/***************************************************************************/
void init_cross( double ***Fou) {
  /***************************************************************************/
  /* Initialisation de la micro-ondelette d'analyse-reconstruction
   * =======================================================================
   *	Matrix for Cross Fourier Transform (real and imaginary parts)	
   * ------------------------------------|----------------------------------|
   *           Fou[0] (real)             |            Fou[1] (imag.)        | 
   * ------------------------------------|----------------------------------|
   *                                     |                                  |
   *    1/3   1/3    1/3    1/3    1/3   |     a   0    0    0    0         |
   *    1/3  -1/6   -1/6    1/3    1/3   |     0   a   -a    0    0         |
   *    1/3  -1/6   -1/6    1/3    1/3   |     0  -a    a    0    0         |
   *    1/3   1/3    1/3   -1/6   -1/6   |     0   0    0    a   -a         |
   *    1/3   1/3    1/3   -1/6   -1/6   |     0   0    0   -a    a         |  
   *                                     |                                  |
   *                                     | a=sqrt(3.)/6.                    |
   * ------------------------------------|----------------------------------|
   * ------------------------------------|----------------------------------|
   *  Inversa matrix for Cross Fourier Transform (real and imaginary parts)
   * ------------------------------------|----------------------------------|
   *           Fou[2] (real)             |           Fou[3] (imag.)             
   * ------------------------------------|----------------------------------|
   *                                     |                                  |
   *   -1    1      1      1      1      |     0   0    0    0    0         |
   *    1   -1/2   -1/2    0      0      |     0  -b    b    0    0         |
   *    1   -1/2   -1/2    0      0      |     0   b   -b    0    0         | 
   *    1    0      0     -1/2   -1/2    |     0   0    0   -b    b         |
   *    1    0      0     -1/2   -1/2    |     0   0    0    b   -b         |
   *                                     |                                  |
   *                                     | b=sqrt(3.)/2.                    |
   * ------------------------------------|----------------------------------|
   * =======================================================================
   * Ces matrices servent pour la convolution locale avec le signal original
   * (voir Fourier_cross). Noter comment les coeffcients de la matrice 
   * interviennent dans la convolution.
   * =======================================================================
   * Called by genera_medida_UPM_rel
   *           genera_medida_UPM 
   *           genera_medida_UPM_corr
   */
  /***************************************************************************/

  int ipx,ipy;
  
  /*	Matrix for Cross Fourier Transform	*/
  for(ipx=0; ipx<5; ipx++)  
    for(ipy=0; ipy<5; ipy++)	{
      Fou[0][ipx][ipy]=1./3.;
      Fou[1][ipx][ipy]= Fou[2][ipx][ipy]= Fou[3][ipx][ipy]=0.;
    }
  
  for(ipx=1; ipx<5; ipx++)  {	
    Fou[0][ipx][ipx]=-1./6.;
    Fou[1][ipx][ipx]= sqrt(3.)/6.;
  }
  
  Fou[0][1][2]= Fou[0][2][1]= Fou[0][3][4]= Fou[0][4][3]= -1./6.;
  Fou[1][1][2]= Fou[1][2][1]= Fou[1][3][4]= Fou[1][4][3]= -sqrt(3.)/6.;
  
  /*	Inversa matrix for Cross Fourier Transform	*/
  
  /*	This line is the only depending on the dimensionality 	*/
  Fou[2][0][0]=-1.;
  /*	and it is irrelevant for zero mean signals		*/
  
  for(ipx=1; ipx<5; ipx++) {
    Fou[2][ipx][0]= Fou[2][0][ipx]= 1.;
    Fou[2][ipx][ipx] = -1./2.;
    Fou[3][ipx][ipx] = -sqrt(3.)/2.;
  }
  /* Fou[2][1][1]= Fou[2][2][2]= Fou[2][3][3]= Fou[2][4][4]= -1./2.;
     Fou[3][1][1]= Fou[3][2][2]= Fou[3][3][3]= Fou[3][4][4]= -sqrt(3.)/2.;
  */
  
  Fou[2][1][2]= Fou[2][2][1]= Fou[2][3][4]= Fou[2][4][3]= -1./2.;
  Fou[3][1][2]= Fou[3][2][1]= Fou[3][3][4]= Fou[3][4][3]= sqrt(3.)/2.;

} // end of init_cross


/****************************************************************************/
void deriva_cross( double *Rc, double *Ic, double *Rgx, double *Igx, 
		  double *Rgy, double *Igy, double ***Fourier) {
  /****************************************************************************/
  /* Reconstruction du patch local a l'aide de la micro-ondelette.
   * - Rc: patch des valeurs reelles du signal dans un voisinage du pixel 
   *   d'interet;
   * - Ic: valeurs imaginaires de ce signal;
   * =========================================================================
   * Called by patchea_vec.   */
  /****************************************************************************/
  int ip;
  double aux,diff;

  /*	Definition of the differential element	*/
  /*	Naif element:	diff=2.*PI/3.;  */
  diff=sqrt(3.); /* id est: 2.*sin(PI/3) */

  IFDEBUG {
    WarnMsg("BEFORE Cross Fourier Transform of (Rgx,Igx)");
    for(ip=0; ip<5; ip++)  
      fprintf(stderr,"\nip=%d => Rc=%f, Ic=%f, Rgx=%f, Igx=%f, Rgy=%f, Igy=%f",
	      ip,Rc[ip],Ic[ip],Rgx[ip],Igx[ip],Rgy[ip],Igy[ip]);
  }

  /*	Initialization				*/
  for(ip=0; ip<5; ip++)    {
    /* Assignation des valeurs du signal d'entree a Rgx et Igx */
    Rgx[ip]=Rc[ip];
    Igx[ip]=Ic[ip];
    /* Initialisation de Rgy et Igy */
    Rgy[ip]= Igy[ip]= 0.;
  }
  
  /* Cross Fourier Transform of (Rgx,Igx) */
  Fourier_cross(Fourier, Rgx, Igx, -1);

  IFDEBUG {
    WarnMsg("AFTER Cross Fourier Transform of (Rgx,Igx)");
    for(ip=0; ip<5; ip++)  
      fprintf(stderr,"\nip=%d => Rc=%f, Ic=%f, Rgx=%f, Igx=%f, Rgy=%f, Igy=%f",
	      ip,Rc[ip],Ic[ip],Rgx[ip],Igx[ip],Rgy[ip],Igy[ip]);
  }
  
  /* Rappel sur les indices de la clique ("patch") 
   *                         [iy-1][ix]  <- 3 
   *     1 -> [iy][ix-1]      [iy][ix]  <- 0       [iy][ix+1]  <- 2
   *                         [iy+1][ix]  <- 4
   */
  
  Rgx[0]= Igx[0]= 0.;	
  /* ici, comme si: Rgy[0] = Igy[0]= 0.; */
  
  /* Multiplication du nombre complexe avec l'element differentiel 
   * et avec i=sqrt(-1) (d'ou la rotation) */
  aux=Rgx[1];
  Rgx[1] = -diff*Igx[1];
  Igx[1] = diff*aux;
  /* ici, comme si: Rgy[1] = Igy[1]= 0.; */

  aux=Rgx[2];
  Rgx[2] = diff*Igx[2];
  Igx[2] = -diff*aux;
  /* ici, comme si: Rgy[2] = Igy[2]= 0.; */

  Rgy[3] = -diff*Igx[3];
  Igy[3] = diff*Rgx[3];
  Rgx[3] = Igx[3]= 0.;

  Rgy[4] = diff*Igx[4];
  Igy[4] = -diff*Rgx[4];
  Rgx[4] = Igx[4]= 0.;

  /* Pour Rgx, Igx: indices 0, 3 et 4 => nul   
   * Pour Rgy, Igy: indices 0, 1 et 2 => nul  */

  IFDEBUG {
    WarnMsg("BEFORE Inverse Cross Fourier Transform of (Rgx,Igx) and (Rgy,Igy)");
    for(ip=0; ip<5; ip++)  
      fprintf(stderr,"\nip=%d => Rc=%f, Ic=%f, Rgx=%f, Igx=%f, Rgy=%f, Igy=%f",
	      ip,Rc[ip],Ic[ip],Rgx[ip],Igx[ip],Rgy[ip],Igy[ip]);
  }

  /* Inverse Cross Fourier Transform */
  Fourier_cross(Fourier, Rgx, Igx, 1); 
  Fourier_cross(Fourier, Rgy, Igy, 1); 
  /* (Rgx,Igx) et (Rgy,Igy) contiennent l'information de gradient
   * dans les directions x et y respectivement */
  
  IFDEBUG {
    WarnMsg("AFTER  Inverse Cross Fourier Transform of (Rgx,Igx) and (Rgy,Igy)");
    for(ip=0; ip<5; ip++)   
      fprintf(stderr,"\nip=%d => Rc=%f, Ic=%f, Rgx=%f, Igx=%f, Rgy=%f, Igy=%f",
	      ip,Rc[ip],Ic[ip],Rgx[ip],Igx[ip],Rgy[ip],Igy[ip]);
  }

}  // end of deriva_cross


/****************************************************************************/
void reconstruct_cross(double *Rgx, double *Igx, double *Rgy, double *Igy,
		      double *Rc, double *Ic, double ***Fourier) {
  /****************************************************************************/
  /* Reconstruction of a local patch with the micro-wavelet defined by init_cross.
   * =========================================================================
   * Called by patchea_vec. */
  /****************************************************************************/
  double diff;
  int ip;
    
  /*	Definition of the differential element	*/
  /*	Naif element:	diff=2.*PI/3.; */
  diff=sqrt(3.); /* id est: 2.*sin(PI/3) */
  
  /*	Processing				*/
  
  IFDEBUG {
    WarnMsg("BEFORE Cross Fourier Transform of (Rgx,Igx) and (Rgy,Igy)");
    for(ip=0; ip<5; ip++)    
      fprintf(stderr,"\nip=%d => Rc=%f, Ic=%f, Rgx=%f, Igx=%f, Rgy=%f, Igy=%f",
	      ip,Rc[ip],Ic[ip],Rgx[ip],Igx[ip],Rgy[ip],Igy[ip]);
  }

  /* Cross Fourier Transform of (Rgx,Igx) */
  Fourier_cross(Fourier, Rgx, Igx, -1);
  /* Cross Fourier Transform of (Rgy,Igy) */
  Fourier_cross(Fourier, Rgy, Igy, -1);
  
  IFDEBUG {
    WarnMsg("AFTER Cross Fourier Transform of (Rgx,Igx) and (Rgy,Igy)");
    for(ip=0; ip<5; ip++)    
      fprintf(stderr,"\nip=%d => Rc=%f, Ic=%f, Rgx=%f, Igx=%f, Rgy=%f, Igy=%f",
	      ip,Rc[ip],Ic[ip],Rgx[ip],Igx[ip],Rgy[ip],Igy[ip]);
  }
  
  Rc[0]= Ic[0]= 0.;
    
  Rc[1]=Igx[1]/diff;
  Ic[1]=-Rgx[1]/diff;

  Rc[2]=-Igx[2]/diff;
  Ic[2]=Rgx[2]/diff;

  Rc[3]=Igy[3]/diff;
  Ic[3]=-Rgy[3]/diff;

  Rc[4]=-Igy[4]/diff;
  Ic[4]=Rgy[4]/diff;
  
  IFDEBUG {
    WarnMsg("BEFORE Inverse Cross Fourier Transform of (Rc,Ic)");
    for(ip=0; ip<5; ip++)    
      fprintf(stderr,"\nip=%d => Rc=%f, Ic=%f, Rgx=%f, Igx=%f, Rgy=%f, Igy=%f",
	      ip,Rc[ip],Ic[ip],Rgx[ip],Igx[ip],Rgy[ip],Igy[ip]);
  }

  /* Inverse Cross Fourier Transform of Rgx */
  Fourier_cross(Fourier, Rc, Ic, 1);

  IFDEBUG {
    WarnMsg("AFTER Inverse Cross Fourier Transform of (Rc,Ic)");
    for(ip=0; ip<5; ip++)    
      fprintf(stderr,"\nip=%d => Rc=%f, Ic=%f, Rgx=%f, Igx=%f, Rgy=%f, Igy=%f",
	      ip,Rc[ip],Ic[ip],Rgx[ip],Igx[ip],Rgy[ip],Igy[ip]);
  }

} // end of reconstruct_cross


/****************************************************************************/
void Fourier_cross( double ***Fourier, double *Rf, double *If, int sign) {
  /****************************************************************************/
  /* Cross Fourier Transform or Inverse Cross Fourier Transform.
   * =========================================================================
   * Called by reconstruct_cross 
   *           deriva_cross  */
  /****************************************************************************/
  
  int ipx,ipy,iFou;
  double Raux[5], Iaux[5];
  double Ra,Ia;
  
  iFou = (int)((sign<0) ? 0 : 2);
  /*  if(sign<0) iFou=0; else iFou=2; */
  
  for(ipx=0; ipx<5; ipx++) Raux[ipx]= Iaux[ipx]= 0.;

  for(ipx=0; ipx<5; ipx++)    {
    for(ipy=0; ipy<5; ipy++)      {
      /* Multiplication complexe pour la convolution de (Rf+j*If) avec 
       * Fourier[ifou] + j*Fourier[ifou+1]    */
      C_mult( Rf[ipy], If[ipy], Fourier[iFou][ipy][ipx], 
	      Fourier[iFou+1][ipy][ipx], &Ra, &Ia);
      Raux[ipx]+=Ra;
      Iaux[ipx]+=Ia;
    }
  }
  
  for(ipx=0; ipx<5; ipx++)  {
    Rf[ipx]=Raux[ipx];
    If[ipx]=Iaux[ipx];
  }
} // end of Fourier_cross


/***************************************************************************/
void recorta( int dimx, int dimy, int ix, int iy, double **signal, 
	      double *parche) {
  /****************************************************************************/
  /* Extraction d'un "patch" (clique ou voisinage) en croix du signal:
   *                         
   *                           [iy-1][ix]  <- 3 
   *     1 -> [iy][ix-1]        [iy][ix]  <- 0       [iy][ix+1]  <- 2
   *                           [iy+1][ix]  <- 4
   *                  
   * en tenant compte des bords et de la periodicite du signal  
   * ========================================================================
   * Called by genera_medida_UPM_rel
   *           genera_medida_UPM
   *           genera_medida_UPM_corr */
  /***************************************************************************/
  
  int iix,isx,iiy,isy;
  
  /* if(ix>0) iix=ix-1; else iix=dimx-1;
     if(ix<dimx-1) isx=ix+1; else isx=0;
     if(iy>0) iiy=iy-1; else iiy=dimy-1;
     if(iy<dimy-1) isy=iy+1; else isy=0; */
  
  parche[0] = signal[iy][ix] ;
  parche[1] = signal[iy][(ix>0) ? (ix-1) : (dimx-1)];
  parche[2] = signal[iy][(ix<dimx-1) ? (ix+1) : 0];
  parche[3] = signal[(iy>0) ? (iy-1) : (dimy-1)][ix];
  parche[4] = signal[(iy<dimy-1) ? (iy+1) : 0][ix];
  
  /* parche[0]=signal[iy][ix];
     parche[1]=signal[iy][iix];
     parche[2]=signal[iy][isx];
     parche[3]=signal[iiy][ix];
     parche[4]=signal[isy][ix]; */
} // end of recorta


/****************************************************************************/
double patchea( double *parche, double ***Fourier) {
  /****************************************************************************/
  /* Called by genera_medida_UPM  */
  /****************************************************************************/
  
  double *Rgx,*Igx,*Rgy,*Igy,*Iparche;
  double *Rax,*Iax,*Ray,*Iay;
  double mediap=0.,disp=0.;
  int ip;

  TrackNullAlloc( Rgx=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Igx=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Rgy=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Igy=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Iparche=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Rax=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Iax=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Ray=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Iay=(double*)calloc(5,sizeof(double)) );

  /* somme du signal dans le voisinage considere */
  for(ip=0; ip<5; ip++) mediap += parche[ip];

  mediap = mediap/3.;
  parche[0] += mediap;
  for(ip=1; ip<5; ip++) parche[ip] -= mediap;

  /* Calcul du signal derive avec la micro-ondelette */
  deriva_cross(parche,Iparche,Rgx,Igx,Rgy,Igy,Fourier);
	
  for(ip=0; ip<5; ip++) {
    Rax[ip] = Rgx[ip];
    Iax[ip] = Igx[ip];
    Ray[ip] = Rgy[ip];
    Iay[ip] = Igy[ip];
  }
  Rgx[0]= Igx[0]= Rgy[0]= Igy[0]=0.;

  /* Reconstruction du signal avec la micro-ondelette */
  reconstruct_cross(Rgx,Igx,Rgy,Igy,parche,Iparche,Fourier);
  deriva_cross(parche,Iparche,Rgx,Igx,Rgy,Igy,Fourier);

  for(ip=0; ip<5; ip++)  {
    Rax[ip] -= Rgx[ip];
    Iax[ip] -= Igx[ip];
    Ray[ip] -= Rgy[ip];
    Iay[ip] -= Igy[ip];
  }
  /* on retourne le module */
  disp = sqrt(Rax[0]*Rax[0] + Ray[0]*Ray[0]);
  
  free(Rgx); free(Igx);
  free(Rgy);  free(Igy);
  free(Rax);  free(Iax);
  free(Ray);  free(Iay);
  free(Iparche);
	
  return disp;
} // end of patchea


/****************************************************************************/
int patchea_vec( double *parche, double ***Fourier, double *errx,
		  double *erry) {
  /****************************************************************************/
  /* Convolution with the micro-wavelet defined by init_cross.
   * Parameters:
   *     - patch (neighborhood) extracted from the signal,
   *     - Fourier : micro-wavelet defined in Direct and in Fourier space 
   *       for decomposition/reconstruction 
   * ========================================================================== 
   * Called by genera_medida_UPM_rel 
   *           genera_medida_UPM_corr
   * See also patchea_vec_optimal */
  /****************************************************************************/
  
  double *Rgx,*Igx,*Rgy,*Igy,*Iparche;
  double *Rax,*Iax,*Ray,*Iay;
  double mediap=0.;
  int ip;

  TrackNullAlloc( Rgx=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Igx=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Rgy=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Igy=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Iparche=(double*)calloc(5,sizeof(double)) );

  TrackNullAlloc( Rax=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Iax=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Ray=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Iay=(double*)calloc(5,sizeof(double)) );

  /* Somme du signal dans le voisinage considere */
  for(ip=0; ip<5; ip++) mediap+=parche[ip];

  /* Ponderation locale du voisinage */
  mediap = mediap/3.;
  parche[0] += mediap;
  for(ip=1; ip<5; ip++) parche[ip]-=mediap;
  /* On a alors :  
   *      parche[0] = - sum_{ip=1}^4 parche[ip] 
   *  =>  sum_{ip=0}^4 parche[ip] = 0                 */
  
  /* Calcul du signal derive avec la micro-ondelette */
  deriva_cross(parche,Iparche,Rgx,Igx,Rgy,Igy,Fourier);
  
  for(ip=0; ip<5; ip++) {
    Rax[ip] = Rgx[ip];
    Iax[ip] = Igx[ip];
    Ray[ip] = Rgy[ip];
    Iay[ip] = Igy[ip];
  }

  Rgx[0]= Igx[0]= Rgy[0]= Igy[0]=0.;

  IFDEBUG {
    WarnMsg("\nintermedary results");
    for(ip=0; ip<5; ip++)    
      fprintf(stderr,"\nip=%d => Rax=%f, Iax=%f, Ray=%f, Iay=%f",
	      ip,Rax[ip],Iax[ip],Ray[ip],Iay[ip]);
  }
 
  /* Reconstruction du signal avec la micro-ondelette */
  reconstruct_cross(Rgx,Igx,Rgy,Igy,parche,Iparche,Fourier);
  deriva_cross(parche,Iparche,Rgx,Igx,Rgy,Igy,Fourier);
  
  IFDEBUG {
    WarnMsg("\nresults");
    for(ip=0; ip<5; ip++)    
      fprintf(stderr,"\nip=%d => Rgx=%f, Igx=%f, Rgy=%f, Igy=%f",
	      ip,Rgx[ip],Igx[ip],Rgy[ip],Igy[ip]);
  }
  
  for(ip=0; ip<5; ip++) {
    Rax[ip] -= Rgx[ip];
    Iax[ip] -= Igx[ip];
    Ray[ip] -= Rgy[ip];
    Iay[ip] -= Igy[ip];
  }
  
  *errx=Rax[0];
  *erry=Ray[0];
  
  free(Rgx);  free(Igx);
  free(Rgy);  free(Igy);
  free(Rax);  free(Iax);
  free(Ray);  free(Iay);
  free(Iparche);	

  return OK;
} // end of patchea_vec



/***************************************************************************/
int patchea_vec_optimal( double *parche, double ***Fourier, double *errx,
			  double *erry) {
  /***************************************************************************/
  /* Version patchea_vec minimale: tous les calculs ne sont pas effectues */
  /***************************************************************************/
  double *Rgy,*Igy,*Iparche;
  double mediap=0.;
  int ip;
  double aux;
  double diff=sqrt(3.); /* id est: 2.*sin(M_PI/3) */
  
  TrackNullAlloc( Rgy=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Igy=(double*)calloc(5,sizeof(double)) );
  TrackNullAlloc( Iparche=(double*)calloc(5,sizeof(double)) );
  
  /* Somme du signal dans le voisinage considere */
  for(ip=0; ip<5; ip++) mediap+=parche[ip];
  
  /* Ponderation locale du voisinage */
  mediap = mediap/3.;
  parche[0] += mediap;
  for(ip=1; ip<5; ip++) parche[ip]-=mediap;
  
  IFDEBUG  {
    WarnMsg("BEFORE Cross Fourier Transform of (parche,Iparche)");
    for(ip=0; ip<5; ip++)  
      fprintf(stderr,"\nip=%d => parche=%f, Iparche=%f",ip,parche[ip], Iparche[ip]);
  }  
  
  /* Calcul du signal derive avec la micro-ondelette */
  Fourier_cross(Fourier, parche, Iparche, -1);
  
  IFDEBUG  {
    WarnMsg("AFTER Cross Fourier Transform of (parche,Iparche)");
    for(ip=0; ip<5; ip++)  
      fprintf(stderr,"\nip=%d => parche=%f, Iparche=%f",ip,parche[ip], Iparche[ip]);
  }  
  
  /* Composante en x */
  parche[0] = Iparche[0] = 0.;	
  aux=parche[1];
  parche[1] = -diff*Iparche[1]; 
  Iparche[1] = diff*aux;
  aux=parche[2];
  parche[2] = diff*Iparche[2]; 
  Iparche[2] = -diff*aux;
  
  /* Composante en y */
  for(ip=0; ip<3; ip++) Rgy[ip] = Igy[ip] = 0.;
  Rgy[3] = -diff*Iparche[3]; 
  Igy[3] = diff*parche[3];
  Rgy[4] = diff*Iparche[4];  
  Igy[4] = -diff*parche[4];
  
  parche[3] = Iparche[3]= 0.;
  parche[4] = Iparche[4]= 0.; 
  
  *errx = parche[1]+parche[2];
  *erry = Rgy[3] + Rgy[4];
  
  IFDEBUG  {
    WarnMsg("BEFORE Inverse Cross Fourier Transform of (parche,Iparche) and (Rgy,Igy)");
    for(ip=0; ip<5; ip++)  
      fprintf(stderr,"\nip=%d => parche=%f, Iparche=%f, Rgy=%f, Igy=%f",
	      ip,parche[ip], Iparche[ip],Rgy[ip],Igy[ip]);
    WarningVV("\nintermedary results: errx = %f, erry = %f",*errx,*erry);
  }
  
  Fourier_cross(Fourier, parche, Iparche, 1); 
  Fourier_cross(Fourier, Rgy, Igy, 1);  
  
  IFDEBUG  {
    WarnMsg("AFTER Inverse Cross Fourier Transform of (parche,Iparche) and (Rgy,Igy)");
    for(ip=0; ip<5; ip++)  
      fprintf(stderr,"\nip=%d => parche=%f, Iparche=%f, Rgy=%f, Igy=%f",
	      ip,parche[ip], Iparche[ip],Rgy[ip],Igy[ip]);
  }
  
  parche[0]= Iparche[0]= Rgy[0]= Igy[0]=0.;
  Fourier_cross(Fourier, parche, Iparche, -1);
  Fourier_cross(Fourier, Rgy, Igy, -1);
  
  IFDEBUG  {
    WarnMsg("AFTER last Cross Fourier Transform of (parche,Iparche) and (Rgy,Igy)");
    for(ip=0; ip<5; ip++)  
      fprintf(stderr,"\nip=%d => parche=%f, Iparche=%f, Rgy=%f, Igy=%f",
	      ip,parche[ip], Iparche[ip],Rgy[ip],Igy[ip]);
  }
  
  *errx += parche[0];
  *erry += Rgy[0];
  
  free(Rgy);  free(Igy);
  free(Iparche);	
}


/***************************************************************************/
double mascara( int dimx, int dimy, int sign, double umbral, 
		double **disp, char **msm ) {
  /***************************************************************************/
  /* Computes the MSM and its density according to the distribution of the 
   * singularity exponents.
   * =========================================================================
   * Called by compute_UPM */
  /****************************************************************************/
  
  int ix,iy;
  double buff,dens=0.;
  double media_tot=0.;
  
  for(iy=0; iy<dimy; iy++) 
    for(ix=0; ix<dimx; ix++) {
      msm[iy][ix] = C0;
      /* C0 = grey level associated to 0 */
      media_tot += disp[iy][ix];
    }
  
  media_tot /= (dimy*dimx);
  
  for(iy=0; iy<dimy; iy++)
    for(ix=0; ix<dimx; ix++)
      if(sign*disp[iy][ix] >= sign*umbral*media_tot)	{
	msm[iy][ix] = CP;
	/* CP = grey level associated to 1 */
	dens+=1.;
     }
  
  dens /= ((double)dimx*dimy);
  
  return dens;
} // end of mascara


/***************************************************************************/
void marca_limits( int dimx, int dimy, double **signal, double *min_c,
		    double *max_c) {
  /***************************************************************************/
  
  double suma=0.;
  int ix,iy;

  *max_c=-1e30;
  *min_c=1e30;
  for(ix=0;ix<dimx;ix++) 
    for(iy=0;iy<dimy;iy++) {
      suma += signal[iy][ix];
      *max_c = fMax(*max_c,signal[iy][ix]);
      *min_c = fMin(*min_c,signal[iy][ix]);
    }
  
  suma = suma / ((double) dimx*dimy);
  *max_c -= suma;
  *min_c -= suma;
} // end of marca_limits



/***************************************************************************/
void singularize_sign( int dimx, int dimy, int xeff, int yeff, int modo, 
		       char **msm, double **gx, double **gy) {
  /***************************************************************************/
  /* Computes a signed MSM or not (according to the orientation of
   * the gradient over MSM) and signed or not gradient's vectors 
   * according to this orientation.
   * =========================================================================
   * Called by determine_MSM (where modo=0)
   *           deriva_orient_limits (where modo=1)
   * =========================================================================
   * Remarks:
   *   1) in compute_dummy_msm: [gx,gy]=deriva(signal): gradient de l'image  
   *   originale;
   *   2) in deriva_orient_limits: [gx,gy]=deriva_limits(MSM);
   */   
  /***************************************************************************/
  
  double **ax,**ay;
  double sigproy;
  int ix,iy;
  
  TrackNullAlloc( ax=matrix2D(yeff,xeff) );
  TrackNullAlloc( ay=matrix2D(yeff,xeff) );
  
  deriva_limits(dimx,dimy,xeff,yeff,msm,ax,ay);
  
  /* Calcul du produit scalaire des vecteurs (ax,ay) et (gx,gy) sur l'image */
  for(ix=0;ix<dimx;ix++)     
    for(iy=0;iy<dimy;iy++) 
      if(msm[iy][ix] == C0) // C0 = grey level associated to 0
	gx[iy][ix]= gy[iy][ix]= 0.;
      else {
	sigproy = ax[iy][ix]*gx[iy][ix] +  
	  ay[iy][ix]*gy[iy][ix];
	/*  sigproy renseigne sur l'orientation du gradient (gx,gy) par rapport
	 * au gradient unitaire (ax,ay) sur les pixels de la MSM */
	if(modo == 0) { 
	  /* Give a sign to the MSM according to the orientation of the 
	   * gradient */
	  if(sigproy<0) msm[iy][ix] = CM;
	  /* CM = grey level associated to -1 */
	}  else if( ((sigproy>=0.) && (msm[iy][ix]==CM))||
		    ((sigproy<0.) && (msm[iy][ix]==CP)) )  {
	  /* Orientate the gradient and leave the MSM unchanged */
	  gx[iy][ix] = -gx[iy][ix];
	  gy[iy][ix] = -gy[iy][ix];
	  /* CP = grey level associated to 1 */
	}
      }
  
  if(xeff>dimx)   
    for(ix=dimx;ix<xeff;ix++)	
      for(iy=0;iy<yeff;iy++)	    
	gx[iy][ix]= gy[iy][ix]= 0.;
  if(yeff>dimy)  
    for(ix=0;ix<xeff;ix++)     
      for(iy=dimy;iy<yeff;iy++)
	gx[iy][ix]= gy[iy][ix]= 0.;
  
  free_matrix2D(ax,yeff);
  free_matrix2D(ay,yeff);
  
  return OK;
} // end of singularize_sign


/***************************************************************************/
void singularize(int dimx, int dimy, int xeff, int yeff, char **msm, 
		 double **gx, double **gy) {
  /***************************************************************************
  * Applies the mask of MSM on the gradient images. 
   * gx=gy=0 on pixels outside the MSM (ie. when, MSM=C0) and outside the 
   * dimension of the MSM.
   * =========================================================================
   * Called by compute_UPM 
   *           deriva_orient_limits 
   *           kernelea
   *           deskernelea
   ***************************************************************************/

  int ix,iy;
  
  for(iy=0;iy<dimy;iy++) 
    for(ix=0;ix<dimx;ix++)      
      if(msm[iy][ix] == C0)	
	/* C0 = grey level associated to 0 */
	gx[iy][ix]= gy[iy][ix]= 0.;
  
  if(xeff>dimx)  
    for(iy=0;iy<yeff;iy++)	
      for(ix=dimx;ix<xeff;ix++)
	gx[iy][ix]= gy[iy][ix]= 0.;
  if(yeff>dimy)
    for(iy=dimy;iy<yeff;iy++)
      for(ix=0;ix<xeff;ix++)
	gx[iy][ix]= gy[iy][ix]= 0.;
  
} // end of singularize


/****************************************************************************/
int genera_sources( int dimx, int dimy, int tipo, int sign,  double **gx,
		     double **gy ) {
  /****************************************************************************/
  
  double **dR,**dI;
  double norm;
  double x4,y4;
  double bgx,bgy;
  double modd;
  double x,y;
  double conj;
  int ix,iy;

  TrackNullAlloc( dR=matrix2D(dimy,dimx) );
  TrackNullAlloc( dI=matrix2D(dimy,dimx) );
  
  if(tipo==2) conj=1.; else conj=-1.;
  
  /*             Definicion a medio pixel          */
  norm=sqrt((double)(dimx*dimy));
  for(iy=0;iy<dimy;iy++)    {
    y=((double)iy+0.5)/((double)dimy);
    if(iy>=dimy/2) y-=1.;
    
    for(ix=0;ix<dimx;ix++)	{
      
      x=((double)ix+0.5)/((double)dimx);
      if(ix>=dimx/2) x-=1.;
      
      x4 = x*x*x*x + y*y*y*y - 6.*x*x*y*y - 1.;
      y4 = 2.*x*y*(x*x-y*y);
      
      modd = (x*x+y*y) * (x4*x4+y4*y4);
      if(modd>1e-30)	    {
	bgx = -(x*x4-y*y4)/modd;
	bgy = (x*y4+y*x4)/modd;
	dR[iy][ix] = bgx*(x4+5./4.)-bgy*y4;
	dI[iy][ix] = conj*(bgx*y4+bgy*(x4+5./4.));
      }      else	
	dR[iy][ix]= dI[iy][ix]= 0.;
      dR[iy][ix] = norm*dR[iy][ix];
      dI[iy][ix] = norm*dI[iy][ix];
    }
  }
  
  Fourier2D(dimx,dimy,dR,dI,-1);
  Fourier2D(dimx,dimy,gx,gy,-1);

  for(iy=0;iy<dimy;iy++)  
    for(ix=0;ix<dimx;ix++)	{      
      modd = dR[iy][ix]*dR[iy][ix]+dI[iy][ix]*dI[iy][ix];
      if(sign>0)	{
	if(modd>1e-30)	  {
	  dR[iy][ix]=dR[iy][ix]/modd;
	  dI[iy][ix]=-dI[iy][ix]/modd;
	} else 
	  dR[iy][ix]= dI[iy][ix]= 0.;
      }
      if(tipo==0)	{
	dR[iy][ix]=sqrt(dR[iy][ix]*dR[iy][ix]+
			dI[iy][ix]*dI[iy][ix]);
	dI[iy][ix]=0.;
      }
      bgx = gx[iy][ix];
      bgy = gy[iy][ix];
      
      gx[iy][ix]=bgx*dR[iy][ix]-bgy*dI[iy][ix];
      gy[iy][ix]=bgy*dR[iy][ix]+bgx*dI[iy][ix];   
    }
  
  Fourier2D(dimx,dimy,gx,gy,1);
  
  free_matrix2D(dR,dimy);
  free_matrix2D(dI,dimy);

  return OK;
} // end of genera_sources

  
/***************************************************************************/
double compute_norma_filtro(int dimx, int dimy, double expon) {
  /***************************************************************************/
  int ix,iy;
  double **aux;
  double mini,norma;
  double f1,f2,cambio;

  TrackNullAlloc( aux=matrix2D(dimy,dimx) );
  /* fill0(dimx,dimy,aux,NULL); */
  mini = aux[0][0] = 1.;
  filtro2D_bak(dimx,dimy,0.,expon,aux);
  
  /* en espagnol dans le texte : Calculo por minimo */
  for(iy=0;iy<dimy;iy++)  
    for(ix=0;ix<dimx;ix++)
      mini = fMin(mini,aux[iy][ix]);
  
  norma = -mini * ((double)dimx*dimy);
  
  /* en espagnol dans le texte : Calculo por relacion 
     f1=aux[dimy/2][dimx/2];
     f2=aux[dimy/4][dimx/4];
     if(fabs(2.-expon)>1e-5)  {
     cambio=pow(2.,2.-expon);
     norma=-(cambio*f1/(cambio-1.)-f2/(cambio-1.))*
     ((double)dimx*dimy);
     }    else  {
     cambio=sqrt((double)(dimx*dimx+dimy*dimy))/2.;
     cambio=log(cambio)/log(2.);
     norma=-(f1*(1.-cambio)+f2*cambio)*((double)dimx*dimy);
     }
  */
  
  return norma;
} // end of compute_norma_filtro


/***************************************************************************/
int  kernelea( int dimx, int dimy, int xeff, int yeff, 
	       char **msm, double **signal) {
  /***************************************************************************/
  double **gx,**gy;
  
  TrackNullAlloc( gx=matrix2D(yeff,xeff) );
  TrackNullAlloc( gy=matrix2D(yeff,xeff) );

  copy(dimx,dimy,signal,gx,NULL);
  deriva2D_bak(xeff,yeff,gx,gy);
  singularize(dimx,dimy,xeff,yeff,msm,gx,gy);
  reconstruct2D_FFT(xeff,yeff,gx,gy);
  copy(dimx,dimy,gx,signal,NULL);
	
  free_matrix2D(gx,yeff);
  free_matrix2D(gy,yeff);

  return OK;
} // end of kernelea


/***************************************************************************/
int deskernelea( int dimx, int dimy, int xeff, int yeff, 
		  char **msm, double **signal) {
  /***************************************************************************/
 
  double **gx,**gy;

  TrackNullAlloc( gx=matrix2D(yeff,xeff) );
  TrackNullAlloc( gy=matrix2D(yeff,xeff) );
  
  copy(dimx,dimy,signal,gx,NULL);
  deriva2D_bak(xeff,yeff,gx,gy);
  singularize(dimx,dimy,xeff,yeff,msm,gx,gy);
  reconstruct2D_FFT(xeff,yeff,gx,gy);
  op_diff(dimx,dimy,gx,signal,NULL);
	
  free_matrix2D(gx,yeff);
  free_matrix2D(gy,yeff);
  
  return OK;
} // end of deskernelea


/***************************************************************************/
int compute_sources( int dimx, int dimy, int xeff, int yeff, double **dummy, 
		      double **gx, double **gy, double *meang ) {
  /***************************************************************************/
  /* Computes the sources of an image, starting from the estimation of the MSM 
   * (namely from the chromatically reduced image) and of the derivatives of the
   * reconstructed image */
  /***************************************************************************/
  
  double **ax,**ay;
  double mmgx[2],mmgy[2];
  double normax,normay,norma;  
  int ix,iy;
  
  TrackNullAlloc( ax=matrix2D(yeff,xeff) );
  TrackNullAlloc( ay=matrix2D(yeff,xeff) );
  
  /*	Procesado				*/
  copy(dimx,dimy,dummy,ax,NULL);
  filtro2D_bak(xeff,yeff,0.,-EXP_MU,ax);
  deriva2D_bak(xeff,yeff,ax,ay);

  IFDEBUG  {
    WarnMsgVV("\ncreate filtro_dummy.gif");
    write_foto_block(dimx,dimy,1.,"filtro_dummy.gif",ax);
  }
  
  filtro2D_bak(xeff,yeff,0.,-EXP_MU,gx);
  filtro2D_bak(xeff,yeff,0.,-EXP_MU,gy);
  opvec_divide(xeff,yeff,ax,ay,gx,gy,NULL);
  
  free_matrix2D(ax,yeff);
  free_matrix2D(ay,yeff);
  
  meang[0] = media(xeff,yeff,gx,NULL);
  meang[1] = media(xeff,yeff,gy,NULL);
  IFVERBOSE WarningVV("Source average: (%f,%f):\n",meang[0],meang[1]);
  
  normax = dispersion(xeff,yeff,gx,NULL);
  normay = dispersion(xeff,yeff,gy,NULL);
  IFVERBOSE WarningV("Source dispersion: %f\n", sqrt(normax*normax+normay*normay));

  extrema(xeff,yeff,gx,mmgx,NULL);
  extrema(xeff,yeff,gy,mmgy,NULL);
  IFVERBOSE WarningVV("Component x: Min.: %f  Max. %f\n",mmgx[0],mmgx[1]);
  IFVERBOSE WarningVV("Component y: Min.: %f  Max. %f\n",mmgy[0],mmgy[1]);
  
  return OK;
} // end of 



/***************************************************************************/
int retrieve_dummy( int dimx, int dimy,  int xeff, int yeff,
		     char **msm, double **dummy) /*recupera_dummy*/{
  /***************************************************************************/
  /* Computes the reconstruction of the image starting from the reduced 
   * unitary essential vectorial field, which :
   *   - is unitary on the MSM, null otherwise
   *   - is perpendicular to the MSM on each pixel of the MSM,
   *   - has identical orientation as the original gradient.
   * =========================================================================
   * Called by the main              */
 /***************************************************************************/

  double **ax,**ay;
  int ix,iy;

  TrackNullAlloc( ax=matrix2D(yeff,xeff) );
  TrackNullAlloc( ay=matrix2D(yeff,xeff) );
  
  /*	Procesado				*/
  deriva_orient_limits(dimx,dimy,xeff,yeff,msm,ax,ay);
  reconstruct2D_FFT(xeff,yeff,ax,ay);  
  copy(dimx,dimy,ax,dummy,NULL);
  
  free_matrix2D(ax,yeff);
  free_matrix2D(ay,yeff);
  
  return OK;
} // end of retrieve_dummy


/***************************************************************************/
int retrieve_essential_gradient( int dimx, int dimy, int xeff, int yeff, 
				 double *expon, double **dummy, 
				 double **gx, double **gy, double *meang)
     /*recupera_gradiente_esencial*/ {
  /***************************************************************************/
  /* Computes the essential gradient field given the chromatically reduced 
   * image. 
   * The essential gradient field is the gradient restricted to the MSM.
   * =========================================================================
   *  Parameters:
   *   - dummy :chromatically reduced image (computed by compute_dummy_msm
   *     or determine_unitary),
   *   - gx, gy: gradient images computed by compute_sources.
   * =========================================================================
   * Called by the main                                                      */ 
  /***************************************************************************/
  
  double **ax,**ay;
  double expon;
  int ix,iy;
  //  expon=EXP_MU;
  
  TrackNullAlloc( ax=matrix2D(yeff,xeff) );
  TrackNullAlloc( ay=matrix2D(yeff,xeff) );
  
  /*	Procesado				*/
  copy(dimx,dimy,dummy,ax,NULL);
  filtro2D_bak(xeff,yeff,0.,-expon,ax);
  deriva2D_bak(xeff,yeff,ax,ay);
  
  opvec_multiply(xeff,yeff,ax,ay,gx,gy,NULL);
  reconstruct2D_FFT(xeff,yeff,gx,gy);

  filtro2D_bak(xeff,yeff,0.,expon,gx);
  deriva2D_bak(xeff,yeff,gx,gy);

  free_matrix2D(ax,yeff);
  free_matrix2D(ay,yeff);

  return OK;
} // end of retrieve_essential_gradient


/***************************************************************************/
int process_sources( int dimx, int dimy, int xeff, int yeff, 
		     double **gx, double **gy) /*procesa_sources*/{
  /***************************************************************************/

  double **mod;
  int ix,iy;
  
  TrackNullAlloc( mod=matrix2D(yeff,xeff) );
  
  for(iy=0;iy<dimy;iy++)  {
    for(ix=0;ix<dimx;ix++)	{
      mod[iy][ix]=sqrt(gx[iy][ix]*gx[iy][ix]+gy[iy][ix]*gy[iy][ix]);
      if(mod[iy][ix]>1e-30) mod[iy][ix]=log(mod[iy][ix]);
      else mod[iy][ix]=-30.*log(10.);
    }
  }
  filtro2D_bak(xeff,yeff,0.,2.,mod);
  // write_foto_4(xeff,yeff,"prueba.gif",mod);
  
  free_matrix2D(mod,yeff);

  return OK;
} // end of procesa_sources


/***************************************************************************/
int haz_histograma( int dimx, int dimy, double *mm, double **data ) {
  /***************************************************************************/

  double pos;
  double *prob;
  int ix,iy,ip;

#ifdef _PARSE_FRACTAL_PARAMETERS_
  nbox = p_msm->nbox;
#else
  nbox = NBOX;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/

  TrackNullAlloc( prob=(double*)calloc(nbox+1,sizeof(double)) );

  extrema(dimx,dimy,data,mm,NULL);
  
  for(iy=0;iy<dimy;iy++)    
    for(ix=0;ix<dimx;ix++)     {
      pos =(data[iy][ix]-min) / (max-min);
      ip = (int)(nbox*pos);
      prob[ip] += 1.;
    }
  
  for(ip=0;ip<=nbox;ip++)
    prob[ip] = prob[ip]*nbox/((double)(dimx*dimy)*(max-min));
  
  /*  strcat(name,".hst");
      write_histograma(name,mm[0],mm[1],prob);
  */
  
  free(prob);
  
  return OK;
} // end of haz_histograma


/***************************************************************************/
void pinta( int dimx, int dimy, char **msm, double **data ) {
  /***************************************************************************/
  int ix,iy;
  
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)
      if(msm[iy][ix] != C0)  data[iy][ix]=1.;
      else                       data[iy][ix]=0.;
  
} // end of pinta


/***************************************************************************/
void pinta_suave( int dimx, int dimy, int xeff, int yeff, 
		  char **msm, double **data ) {
  /***************************************************************************/
  
  int ix,iy;
  
  fill0(xeff,yeff,data,NULL);
  
  /* pinta(dimx,dimy,msm,data); */
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)
      if(msm[iy][ix] != C0)  data[iy][ix] = 1.;
      else                   data[iy][ix] = 0.;
  
  /* Filtrage  dans l'espace des frequences:
   *    - multiplication par 0 pour la composante de frequence nulle,
   *    - multiplication par f^(-1), ou f designe la norme du vecteur 
   *      frequence, pour les autres composantes. */
  filtro2D_bak(xeff,yeff,0.,-1.,data);
  
} // end of pinta_suave


/***************************************************************************/
void deriva_limits( int dimx, int dimy, int xeff, int yeff, char **msm, 
		    double **ax, double **ay) {
  /***************************************************************************/
  /* Transforms the Hausdorff measure of the MSM in the frequency space.
   * Namely, computes:
   *     IFFT( FFT(\delta_{msm} * f^(-1)) ) 
   * where \delta_{msm} stands for the Hausdorff measure restricted to the MSM 
   * and f stands for the norm of the frequency vector.
   * Normalizes the result (unitary vector).
   * =========================================================================
   * Called by:    deriva_orient_limits
   *               singularize_sign         */
  /***************************************************************************/

  double norma;
  
  /*	Calculo de direcccion		*/
  fill0(xeff,yeff,ax,NULL);
  pinta_suave(dimx,dimy,xeff,yeff,msm,ax);

  IFDEBUG  { /*  Represent filtered indicatrice image:  */
    WarnMsg("create filtered image filtro.gif...");
    write_foto_block(dimx,dimy,BLOCKOUT,"filtro.gif",ax);
  }
  
  deriva2D_bak(xeff,yeff,ax,ay);
  opvec_norma(xeff,yeff,ax,ay,NULL);
  
  /* No es preciso singularizar este vector (especializarlo a la MSM)
   * porque todas las rutinas a la fecha lo singularizan o usan solo
   * sobre el soporte de la MSM.
   */

  IFDEBUG {
    WarnMsg("create der_limits.[x,y].gif...");
    write_foto_block(dimx,dimy,BLOCKOUT,"der_limits.x.gif",ax);
    write_foto_block(dimx,dimy,BLOCKOUT,"der_limits.y.gif",ay);
  }
  
} // end of deriva_limits


/***************************************************************************/
void deriva_orient_limits( int dimx, int dimy, int xeff, int yeff, 
			    char **msm, double **ax, double **ay) {
  /***************************************************************************/
  /* Computes the reduced unitary essential vectorial field which:
   *   - is unitary on the MSM, null otherwise,
   *   - is perpendicular to the MSM on each pixel of the MSM,
   *   - has identical orientation as the original gradient.
   * =========================================================================
   * Called by  recupera_dummy   */ 
  /***************************************************************************/

  double norma;
  int ix,iy;
  
  deriva_limits(dimx,dimy,xeff,yeff,msm,ax,ay);
  singularize_sign(dimx,dimy,xeff,yeff,1,msm,ax,ay);

  /* Las siguientes lineas suavizan algo la determinacion de la
   * derivada de la funcion de limits, al imponer compatibilidad
   * con la reconstructibilidad del vector obtenido y suavizacion
   * de los limits que dieron origen a la imagen */
  reconstruct2D_FFT(xeff,yeff,ax,ay);

  deriva2D_bak(xeff,yeff,ax,ay);
  singularize(dimx,dimy,xeff,yeff,msm,ax,ay);
} // end of deriva_orient_limits


/***************************************************************************/
int trunc_reconstruct2D_FFT(int dimx, int dimy, double **gxR, double **gyR) {
  /***************************************************************************
   * Reconstruction en partie seulement du signal: le but est de verifier que 
   * l'expression A^2(f) de la relation suivante dans l'espace des frequences:
   *          S(\vec{f}) = g^2(f) . A^2(vec{f})
   * (\vec{f} est le vecteur frequence, f son module) exprimee dans 
   * l'equation (21) de l'article:
   *    "Reconstructing images from their most singular fractal manifold"
   *     Turiel et del Pozo
   * varie bien en f^\eta. On rappelle (cf. meme article) que:
   *          A(f) = |\vec{v}_\infty \cdot \vec{f}| / f
   * ou \vec{v}_\infty est le gradient essentiel et \cdot designe le produit
   * scalaire.
   * ========================================================================
   * Des parties de la fonction reconstruct2D_FFT ont simplement ete tronquees de 
   * maniere a retourner le scalaire A^2(f).
   ***************************************************************************/
  
  int ix,iy;
  double x,y,agxR,agxI,agyR,agyI;
  double dx,dy,prefm;
  double **gxI, **gyI;
  
  /* (gxR,gyR):  bidimensionnal vector field */
  TrackNullAlloc( gxI=matrix2D(dimy,dimx);
  TrackNullAlloc( gyI=matrix2D(dimy,dimx);
  
  /* ((gxR,gxI) (gyR,gyI)):  Fourier transform of (gxR,gyR) 
   * ie., complex bidimensionnal vector field */
  Fourier2D(dimx,dimy,gxR,gxI,-1);
  Fourier2D(dimx,dimy,gyR,gyI,-1);
  
  for(iy=0;iy<dimy;iy++)  {
    y=((double)iy)/((double) dimy);
    if(iy>=dimy/2) y-=1.;
    if(iy==dimy/2) y=0.;
    
    for(ix=0;ix<dimx;ix++) {
      x=((double)ix)/((double) dimx);
      if(ix>=dimx/2) x-=1.;
      if(ix==dimx/2) x=0.;
      
      /* frequency vector: f=(dx,dy) */
      dx=sin(M_PI*x);      dy=sin(M_PI*y);
      /*  dx=x; dy=y; */
      prefm=dx*dx+dy*dy;
     
      agxR=gxR[iy][ix];
      agxI=gxI[iy][ix];
      agyR=gyR[iy][ix];
      agyI=gyI[iy][ix];
      
      /* The scalar product with the frequency vector f is computed */
      if(prefm>1e-30) {
	
	/* MODIF:
	 * La multiplication par le nombre imaginaire i=sqrt{-1} et la
	 * division par f^2 sont supprimees
	 gxR[iy][ix] = (dx*agxI+dy*agyI) / prefm;
	 gxI[iy][ix] = -(dx*agxR+dy*agyR) / prefm;
	*/      
	gxR[iy][ix] = (dx*agxR+dy*agyR);
	gxI[iy][ix] = (dx*agxI+dy*agyI);
	/* La norme du vecteur complexe ainsi obtenue est sauvee dans gxR */
	gxR[iy][ix] *= gxI[iy][ix];
	/* END_MODIF */
      } else {
	gxR[iy][ix] = gxI[iy][ix] = 0;
      }
      /* MODIF:
       * plus aucun interet...
       gyR[iy][ix] = gyI[iy][ix] = 0.;
       * END_MODIF */
      
    }
  }
  
  /* MODIF:
   * Inverse Fourier transform 
   Fourier2D(dimx,dimy,gxR,gxI,1);
   * La transformee inverse est supprimee
   * END_MODIF */
  
  free_matrix2D(gxI,dimy);
  free_matrix2D(gyI,dimy);
  
  return OK;
} // end of trunc_reconstruct2D_FFT





