#include <stdio.h>
#include <math.h>

#include <utl_alloc.h>
#include <utl_stats.h>
#include <utl_op.h>
		
#include <stats.h>
#include <fft2d.h>		
#include <deriva2d.h>		
#include <wavelet.h>

#include <io_grafic.h>  	

#include <mf_manifold.h>
#include <mf_source.h>
#include <mf_operate.h>


extern int ENTROPY;
extern double EXP_MU;
extern int NBOX;
extern int NPUNTOS;
extern int MULTIFRACTAL;

extern double UPM_DENS;
extern double UPM_THR;
extern double ENTR_THR;

extern int VIDEO;
extern double BLOUT;

extern char C0;
extern char CP;
extern char CM;

/***************************************************************************/
double moda_redef(int leff, double **signal ) {
  /***************************************************************************/

#ifdef _PARSE_FRACTAL_PARAMETERS_
  if(p_frac->dim_space == DIM1D) 
#else
    if(DSPACE == DIM1D) 
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
      return moda1D(leff,signal[0]);
    else 
      return moda(leff,leff,signal);
} // end of moda_redef


/***************************************************************************/
void singularize_sign_deconv( int dimx, int dimy, int xeff, int yeff,
			      char **msm, char **mentropy, 
			      double **gx, double **gy, double **ax, double **ay) {
  /***************************************************************************/
  /* Computes a signed MSM or not (according to the orientation of
   * the gradient over MSM) and signed or not gradient's vectors 
   * according to this orientation.
   */   
  /***************************************************************************/
  
  double sigproy;
  int ix,iy;
  double mod;

  for(ix=0;ix<dimx;ix++)     
    for(iy=0;iy<dimy;iy++) 
      if(msm[iy][ix]==C0)
	/* C0 = grey level associated to 0 */
	gx[iy][ix]= gy[iy][ix]= 0.;
      else  if(mentropy[iy][ix]==C0){
	gx[iy][ix] /= ax[iy][ix];
	gy[iy][ix] /= ay[iy][ix];
      } /* sinon il reste inchange */
  
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
}




/***************************************************************************/
void determine_MSM_deconv( int dimx, int dimy, int dimv, int dimz, int iz, 
			   int xeff, int yeff, int levels, char *ext, 
			   double ***signal, double ***entropy, char ***msm) {
/***************************************************************************/

  double **gx,**gy;
  int ic;

  gx=matrix2D(yeff,xeff);
  gy=matrix2D(yeff,xeff);
  
  for(ic=0;ic<dimv;ic++) {
    /* copy(xeff,yeff,gx,NULL); */
    copy(dimx,dimy,signal[ic],gx,NULL);
    /* Calcul du signal derive */
    deriva(xeff,yeff,gx,gy);
    /* Initialisation de la MSM */
    fill_msm(dimx,dimy,msm[ic]);
    compute_UPM_deconv( dimx,dimy,xeff,yeff,levels,ext,
			signal[ic],entropy[ic],gx,gy,msm[ic] );    
    
    singularize_sign(dimx,dimy,xeff,yeff,0,msm[ic],gx,gy);
  }
  
  /*	Computing the multifractal statistics if needed		*/
  if(MULTIFRACTAL) 
    write_multifractal(dimx,dimy,dimv,dimz,iz,/* xeff,yeff, */ ext,msm,signal);

  
  write_foto_RGB(VIDEO,BLOUT,dimx,dimy,dimv,dimz,iz,"contour",ext,msm,msm,msm);

  free_matrix2D(gx,yeff);
  free_matrix2D(gy,yeff);
}

/***************************************************************************/
double compute_UPM_deconv( int dimx, int dimy, int xeff, int yeff, int levels, 
			   char *ext, double **signal, double **entropy, 
			   double **gx, double **gy, char **msm) {
  /***************************************************************************/
  /* Called by determine_MSM */
  /***************************************************************************/

  double **disp,**errorx,**errory;
  double umbral,entr_umbral;
  double norsignal,norerror;
  double psnr;
  double entr_ave[2],dens;

  disp=matrix2D(dimy,dimx);
  
  norsignal=dispersion(dimx,dimy,signal,NULL);
  errorx=matrix2D(yeff,xeff);
  errory=matrix2D(yeff,xeff); 
  fill_msm(dimx,dimy,msm);
  
  /*** Operations pour la determination de la MSM */
  /*    genera_medida_UPM(dimx,dimy,signal,disp);*/
  /*	genera_medida_UPM_corr(dimx,dimy,signal,gx,gy,disp); */
  /*	genera_medida_UPM_glob(dimx,dimy,signal,gx,gy,disp);	*/
  genera_medida_UPM_rel(dimx,dimy,signal,gx,gy,disp);
  
  /* Calcul du seuil de determination des valeurs des exposants de la MSM 
   * a partir de leur distribution */
  if(UPM_DENS>0.) {
    entr_umbral=ENTR_THR; /* Modifier ... */
    umbral=cuantil_2D(dimx,dimy,disp,UPM_DENS) / media(dimx,dimy,disp,NULL);
  } 
  else {
    entr_umbral=ENTR_THR;
    umbral=UPM_THR;
  }

  /* Determination des points de la MSM */
  entr_ave[0]=entr_ave[1]=dens =0.;
  mascara_deconv( dimx, dimy, umbral, disp, entr_umbral, entropy,
		  msm, &dens, entr_ave);
  
  /*** Operations pour le calcul de la PSNR */
  
  copy(xeff,yeff,gx,errorx,NULL);
  copy(xeff,yeff,gy,errory,NULL);
  /* errorx et erroy contiennent les gradients de l'image */
  
  /* Mise a 0 des valeurs du gradient en dehors de la MSM */
  singularize(dimx,dimy,xeff,yeff,msm,errorx,errory); 
  /* Reconstruction a partir du gradient sur la MSM
   * errorx contient l'image reconstruite */
  reconstruct(xeff,yeff,errorx,errory);
  /* Represent reconstruction: */
  
  /* Signal-to-noise (SNR) measures estimates the quality of the reconstructed 
     /* 1) Compute the difference between original and reconstructed images: */
  op_diff(dimx,dimy,signal,errorx, NULL);
  /* 2) Compute the root mean squared error RMSE:            */ 
  norerror = dispersion(dimx,dimy,errorx,NULL) / 
    ((double)levels)*256.;
  /* 3) Compute the PSNR:  */
  psnr = 20. * log(255./norerror) / log(10.);
  printf("Density of UPM at PSNR %f dB: %f%%\n", psnr,100.*dens); 
  printf("Mean entropy on the UPM: %f bpp\n", entr_ave[0]); 
  printf("             on the complementary subset: %f bpp\n", entr_ave[1]); 
  
  /*	Freeing memory associated to the routine	*/
  free_matrix2D(errorx,yeff);
  free_matrix2D(errory,yeff);
  free_matrix2D(disp,dimy);
  
  return dens;
}


/***************************************************************************/
int mascara_deconv( int dimx, int dimy, double umbral, double **disp, 
		    double entr_umbral, double **entropy, char **msm,
		    double *dens, double *entr_ave) {
  /***************************************************************************/
  /*  */
  /****************************************************************************/
  
  int ix,iy,i;
  double buff;
  double media_tot=0., entr_media_tot=0.;
  int count[2];
  
  count[0]=count[1]=0;
  
  for(iy=0; iy<dimy; iy++) 
    for(ix=0; ix<dimx; ix++) {
      msm[iy][ix]=C0;
      /* C0 = grey level associated to 0 */
      media_tot += disp[iy][ix];
      entr_media_tot += entropy[iy][ix];
    }
  
  media_tot /= (dimy*dimx);
  entr_media_tot /= (dimy*dimx);
  
  for(iy=0; iy<dimy; iy++)
    for(ix=0; ix<dimx; ix++)
      if( disp[iy][ix]>=umbral*media_tot && 
	  entropy[iy][ix]>=entr_umbral*entr_media_tot)	{
	msm[iy][ix]=CP;
	/* CP = grey level associated to 1 */
	(*dens)+=1.;
	entr_ave[0] += entropy[iy][ix];
	count[0]++;
      } else {
	entr_ave[1] += entropy[iy][ix];
	count[1]++;
      }
  
  (*dens) /= ((double)dimx*dimy);
  for (i=0; i<2;i++)
    if(count[i]>0)  (entr_ave[i]) /= ((double)count[i]);
  
  return 0;
}
