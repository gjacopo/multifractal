/********************************************************
 * mf_analysis.h               
 ********************************************************/

#ifndef   	_MF_ANALYSIS_H
# define   	_MF_ANALYSIS_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  

  int create_multifractal( int nseries, char *base );
  int analiza_multifractal( /*inputs*/char *name_in, int dimx, 
			    double *Moms, int nmoms, int *Dist, int ndists );
  
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !MF_ANALYSIS_H */
