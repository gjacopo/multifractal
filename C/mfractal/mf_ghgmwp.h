/* ===================================
** mf_ghgmwp.h
** started on Wed Jan 31 18:31:15 2007 
** ===================================
*/

#ifndef   	_MF_GHGMWP_H_
#define   	_MF_GHGMWP_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

  int Dh_estima_gh( /*inputs*/char *name_in, int dimx, 
		    /*outputs*/double *shift_g, 
		    double *h_g, double *Dh_g, double *errDh_g );
  int Dh_estima_gmwp( /*inputs*/char *name_in, int dimx, 
		      /*outputs*/double *shift_w,
		      double **h_w, double **Dh_w, double **errDh_w ); 
  int Dh_filter_weight( /*input*/ double *histo, double *mh, 
			/*outputs*/ double *histo_r, double *h_r, int *width );
  int Dh_filter( /*input*/ double *histo, double *mh, 
		 /*outputs*/ double *histo_r, double *h_r );
  int Dh_errorbar_weight( double sc0, int Nh,
			  /*inputs*/ double *histo_r, double *h_r, int *width, 
			  /*outputs*/ double *h, double *Dh, double *errDh );
  int Dh_histo_weight( double sc0, double *mh, double *histo );
  int Dh_histo( double sc0, double *mh, double *histo );
  void histogram_accumula( int dimx, int dimy, double *mh, double **expon, 
			   double *histo);

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !_MF_GHGMWP_H_ */

