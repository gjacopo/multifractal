#ifndef   	SIG_STATS_H_
#define   	SIG_STATS_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  int computehisto_signal(Signal *sig, char **mask);
  int histequalize_signal( Signal *s, int bin, char **mask );
  
  int histmatch_signal(Signal *s1, Signal *s2, int bin, char **mask);
  
#ifdef DEBUG
  Signal *test_histmatch(Signal *s, int bin);
#endif 
  
  double rmse_signal(Signal* s1, Signal* s2, char**m);
  int mediadisp_signal( Signal* s, double *med, double *disp, char**m );
  double media_signal( Signal* s, char**m );
  double dispersion_signal( Signal* s, char**m );
  int extrema_signal( Signal* s, char**m );
  int display_extrema_signal(char *text, Signal* s, char **m );
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !SIG_STATS_H_ */
