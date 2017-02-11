/* ===================================
** mf_distributon.h
** ===================================
*/

#ifndef   	_MF_DISTRIBUTON_H_
# define   	_MF_DISTRIBUTON_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  double genera_h( double* prob_levi );
  double theoretical_Dh( double h );
  int theoretical_Deltah( double *hmin, double *hmax );

  double poisson( double lambda );
  double normal_standard( int nbox, double tch );
  int genera_levi( double*prob_levi, double alpha, int nbox, double tch );
  double levi_standard( double *prob_levi, int nbox, double tch );
  double monofractal( double h0, double prob );
  double binomial( double h0, double h1 );

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !_MF_DISTRIBUTON_H_ */
