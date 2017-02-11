#ifndef   	_MF_MANIFOLD_H_
#define   	_MF_MANIFOLD_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

void expon_density( int dimx, int dimy,  char *ext, double minimo,
		    double maximo, double **expon, 
		    double *h_inf, double *delta_h, 
		    const char *nombre_wv );
void expon_density_dif( int dimx, int dimy,  char *ext,	double **expon,
			const char *nombre_wv,int i,
			double *valor, double *disp);
void expon_comparison( int dimx, int dimy,  char *ext, double ***expon);

void multimf_histo( int dimx, int dimy, /* int xeff, int yeff, */
			 FILE *canal, char **msm, double **signal, 
			 double **expon );

void vuelcaresult( int dimx, int dimy, double sc0, double qs,
		   double minimo, double maximo, int N, 
		   const char *nombre_wv, double **expon, char *ext );

void visualize_gris( int dimx, int dimy, char *dext, double **expon, const char *nombre_wv);
double visualize_gris_pieces( int dimx, int dimy, double h_inf, double delta_h, char *dext,
		     double **expon, const char *nombre_wv);
void visualize_color( int dimx, int dimy, char *dext, double **expon, const char *nombre_wv);

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !_MF_MANIFOLD_H_ */
