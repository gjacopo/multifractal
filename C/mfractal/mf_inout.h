/* ===================================
** mf_inout.h
** started on Wed Jan 24 19:28:49 2007 
** ===================================
*/

#ifndef   	_MF_INOUT_H_
#define   	_MF_INOUT_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

  int Dh_read( char *name_in, double *h, double *Dh, double *errDh );
  int line_serie_temp( FILE *canal );
  int column_serie_temp( FILE *canal );
  int read_serie_temp( FILE *canal, int dims, int dimt, double **series );
  int Dh_write( char *name, int Nr, double *h, double *Dh, double *errDh );
  
  int write_geomap ( char *name, double **Q );
  int write_typemap ( char *name, double **Q );
    
  /**/
  int visualise_gris_pieces( int dimx, int dimy, char *dext,
			     double m, double M, double delta_h, 
			     double **data );
  
  int lecture_image( int dimx, int dimy, int dimxinrim, int dimyinrim, 
		     int dimv, int dimz, int iz, int bd, char *name_in, 
		     char *ext, double ***signal, double *med_cr);
  void write_corte( int dimx, int dimy, int dimv, int xeff, int yeff, char *ext, int *levels, 
		    double *med_cr, double **meang, char ***msm, double ***gx, 
		    double ***gy);
  void read_corte( int dimx, int dimy, int dimv, int xeff, int yeff, char *ext, int *levels, 
		   double *med_cr, double **meang, char ***msm, double ***gx, double ***gy);
  
  void write_multifractal( int dimx, int dimy, int dimv, int dimz, 
			   int iz, /* int xeff, int yeff, */ char*ext, 
			   char ***msm, double ***signal );
  void vuelcaresult( int dimx, int dimy,  double sc0, double qs,
		     double minimo, double maximo, int N, 
		     const char *nombre_wv, double **expon, char *ext ); 
  void phase_color( double mx, double my, char *Red, char *Green, char *Blue );
  
  int write_expon_density( double *prob, double *mm, int N,
			   char *ext, char *name_in );
  void process_sources( int dimx, int dimy, int xeff, int yeff, 
			double **gx, double **gy);
  void represent_sources( int dimx, int dimy, int dimv, int dimz, int iz, char *ext, 
			  double ***gx, double ***gy );
void represent_sources_mod( int dimx, int dimy, char *ext, int silog,
			    double **mux, double **muy );
  void prepare_sources_mod( int dimx, int dimy, int silog, 
			    double **mux, double **muy, double **mod );
  void represent_sources_fase( int dimx, int dimy, char *ext,
			     double **mux, double **muy);
  void prepare_sources_fase_block( int dimx, int dimy, double block, double **mux, 
				   double **muy, char **Red, char **Green, char **Blue);
  void save_sources( int dimx, int dimy, int dimv, int dimz, int iz, char *ext, 
		     double ***gx, double ***gy);
void represent_sources_angulo( int dimx, int dimy, int dimv, int dimz, int iz, char *ext, 
			       double ***gx, double ***gy);
  void prepare_sources_angulo( int dimx, int dimy, double block, 
			       double **mux, double **muy, double **ang );
  int write_multifractal_histo( int dimx, int dimy, /* int xeff, int yeff, */
				FILE *canal, char **msm, double **signal, 
				double **expon );
  
  void write_histogram( char *nombre, double minimo, double maximo, 
			double *prob );
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !_MF_INOUT_H_ */

