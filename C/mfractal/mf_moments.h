/* ===================================
** mf_moments.h
** started on Wed Jan 31 18:33:23 2007 
** ===================================
*/

#ifndef   	_MF_MOMENTS_H_
#define   	_MF_MOMENTS_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  int Dh_estima_moments( char *name_in, int ndata,
			 int dimx, int dimy,
			 double *Moms, int nmoms, int *Dist, int ndists, 
			 double *h_m, double *Dh_m, double *errDh_m );
  int moments_accumula( int dimx, int dimy, double **signal,
			double *Mom, int nmoms, int *Dist, int ndists,
			double **moments );
  void moments_taup( double **moments, int nmoms, int *Dist, int ndists,
		     double *taup );
  int legendre_transform( double *moms, int N_m, 
			  double *taup, double **h_m, double **Dh_m,
			  double **errDh_m );
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !_MF_MOMENTS_H_ */

