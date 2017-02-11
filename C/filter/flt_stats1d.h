#ifndef STATS1D_H
#define STATS1D_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

  double cuantil1D( int dimx, double *data, double prob0);
  int quicksort1D( int low, int high, double *data );
  int quicksort1D_ref( int low, int high, int *ref, double *data );
  int partition1D( int low, int high, double *data);
  int partition1D_ref( int low, int high, int *ref, double *data);
  
  double moda1D(int dimx, double *datos)/*moda_1D*/;
  double moda1D_histo( int Nhisto, double *mm, double *histo );
  double dispersion1D( int dim, double *cont, char*m )/*dispersion_lista*/;
  double media1D( int dim, double *cont, char*m );
  double covariance1D( int dimx, double *x, double *y, char*m)/*covarianza_lista*/;
  int extrema1D( int dimx, double *cont, double *mm, char *m )/*extrema_lista*/;

  double anorma1D( int dimx, int dimy, double *cont, char*m )/*anorma_lista*/;
  double anorma11D( int dimx, double *cont, char*m )/*anorma1_lista*/;
  double anorma21D( int dimx, double *cont, char*m );
  int denorma1D( int dimx, double norma, double *cont, char*m ) /*denorma_lista*/;
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
