#ifndef _MF_OPERATE_H
#define _MF_OPERATE_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

void singularize_sign_deconv( int dimx, int dimy, int xeff, int yeff,
			      char **msm, char **mentropy, 
			      double **gx, double **gy, double **ax, double **ay) ;
double compute_UPM_deconv( int dimx, int dimy, int xeff, int yeff, int levels, 
			   char *ext, double **signal, double **entropy, 
			   double **gx, double **gy, char **msm);
int mascara_deconv( int dimx, int dimy, double umbral, double **disp, 
		    double entr_umbral, double **entropy, char **msm,
		    double *dens, double *entr_ave);
void determine_MSM_deconv( int dimx, int dimy, int dimv, int dimz, int iz, 
			   int xeff, int yeff, int levels, char *ext, 
			   double ***signal, double ***entropy, char ***msm);

void truncate_reconstruct(int dimx, int dimy, double **gxR, double **gyR);
void multimf_histo( int dimx, int dimy, /* int xeff, int yeff, */
			 FILE *canal, char **msm, double **signal, double **expon);

void save_sources( int dimx, int dimy, int dimv, int dimz, int iz, char *ext, 
		   double ***gx, double ***gy);
void represent_sources_angulo( int dimx, int dimy, int dimv, int dimz, int iz, char *ext, 
			      double ***gx, double ***gy);
void optimal_patchea_vec( double *parche, double ***Fourier, double *errx,
			double *erry);
void prepare_sources_angulo( int dimx, int dimy, double block, 
			     double **mux, double **muy, double **ang );

double mascara_complementaria( int dimx, int dimy, double umbral, 
			       double **disp, char **msm);

#ifdef __cplusplus
}
#endif		/* __cplusplus */


#endif
