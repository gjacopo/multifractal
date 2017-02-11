/* ===================================
** histo.h
** ===================================
*/
#ifndef   	_HISTO_H_
#define   	_HISTO_H_

#ifndef _HISTO_TYPE_
#define _HISTO_TYPE_
typedef struct myhisto {
  int bin;
  double *h;
  double *cum;
  double *val;  // it contributes to the mean binned val
  int nsample;
  int flnorma;
} Histo;
#endif

#ifdef __cplusplus
extern "C" {
#endif	/* __cplusplus */
  
#ifdef _HISTO_TYPE_
  Histo* alloc_histo() ;
  Histo* init_histo(int nbbin);
  int free_histo(Histo* hist);
  int copy_histo( Histo* hin, Histo* hdest );
  
  int compute_histo( int dimx, int dimy, double **s, 	
		    double *mms, Histo *hist, char **mask );
  int cum_histo( Histo *hist );
  int norma_histo( Histo *hist );
  
  int equalize_intensity( double** s, int xdim, int ydim,
			 double *mms, Histo *hist, char **mask);
  
  int match_intensity( double** s, int xdim, int ydim, double *mms, 
		      Histo *hist, Histo *hfit, char **mask );
  
  int add_histo(Histo *h1, Histo *h2);
  int sub_histo(Histo *h1, Histo *h2);
#endif
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !_HISTO_H_ */
