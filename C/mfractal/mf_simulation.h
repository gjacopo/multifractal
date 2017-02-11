/********************************************************/
/*                                                      */
/*              multimf_generator.h                */
/*            Version: 18 de Octubre, 2004              */
/*                                                      */
/********************************************************/

#ifndef MULTIFRACTAL_GENERATOR_H
#define MULTIFRACTAL_GENERATOR_H


#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  int prepara_multifractal(  );
  
  void genera_multifractal( int leff, double sc0, int wav, int ord_der,
			    double *prob_exp,  double **signal);
  int genera1D_multifractal( int leff, double sc0, int wav, int ord_der, 
			     double *prob_exp, double *serie);
  int genera2D_multifractal( int leff, double sc0, int wav, int ord_der, 
			     double *prob_exp, double **image);
  
  void genera_alpha( int dimeta, int Nwav, double *prob_exp, double *alpha );
  int genera1D_alpha( int dimeta, int Nwav, double *prob_exp, double *alpha0 );
  int genera2D_alpha( int dimeta, int Nwav, double *prob_exp, double *alpha0 );
  
  void genera_base( int leff, char *base );

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
