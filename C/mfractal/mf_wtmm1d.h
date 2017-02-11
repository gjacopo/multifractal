#ifndef   	_FRACTALWTMM1D_H_
# define   	_FRACTALWTMM1D_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  int Dh_estima_wtmm( char *name_in, int dimx, 
		      double *Moms, int nmoms, 
		      double *h_wtmm, double *Dh_wtmm, double *errDh_wtmm );
  int wtmm1d_compute( int dimx, int dimy, double **signal, 
		      int nsc, double *Moms, int nmoms, 
		      double **ExtWTlis, double **Z, int *n_ext, 
		      double *scArray, double **sTq, double**sTqLogT );
  int wtmm1d_estima( int dimx, int nsc, double *Moms, int nmoms, 
		     double **Z, int *n_ext, double *scArray,
		     double *h_wtmm, double *Dh_wtmm, double *errDh_wtmm,
		     double *taup, double **sTq, double**sTqLogT );
  
  int wtmm1D_pf( double *signal, int dimx, int wav, int order,
		 int nq, double *qArray, int nsc, double *scArray, 
		 double **ExtWTlis, double **Z, int *n_ext );
  
  int wtmm1D_extrema_find( double *wtrans, int *wtrans_ind, int dimx, 
			   double sc );
  int pf1D_extrema_track( double *wtrans, int *wtrans_ind, 
			double *maxsig, int *maxsig_ind,
			int jj, int  n_max,
			int nq, double *qArray, int nsc, double *scArray, 
			double *ExtWTlis, double *Z );

  int spectrum1D_direct( double **Z, int dimx, 
			  int nq, double *qArray, int nws, double *wsArray,
			  double *Tauq, double *H, double *Dh );
  int mean1D_canon( double **ExtWTlis, int dimx, int *n_ext,
			 int nq, double *qArray, int nws, double *wsArray,
			 double **sTq, double **sTqLogT 
			 /* double **logSTq */ );
  int spectrum1D_canon( double **sTq, double **sTqLogT, /* double **logSTq */
			 int dimx, 
			 int nq, double *qArray, int nws, double *wsArray,
			 double *Tauq, double *H, double *Dh );

  int dumcompare(const double *d1,const double *d2);

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !_MF_WTMM1D_H_ */
