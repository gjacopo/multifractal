#ifndef SOURCE_FUNCION_H 
#define SOURCE_FUNCION_H 

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */


int determine_msm( int dimx, int dimy, int dimv, int dimz, int iz, int xeff, int yeff,
        int levels, char *ext, double ***signal,  char ***msm);
int determine_unitary( int dimx, int dimy, int dimv, int xeff, int yeff,
			double ***gx, double ***gy, double ***dummy);


void init_cross( double ***Fourier);
void marca_limits( int dimx, int dimy, double **signal, double *min_c,
	double *max_c);

void recorta( int dimx, int dimy, int ix, int iy, double **cont, 
	double *parche);
double patchea( double *parche, double ***Fourier);
void patchea_vec( double *parche, double ***Fourier, double *errx,
	double *erry);

double compute_UPM( int dimx, int dimy, int xeff, int yeff, int levels, char *ext, 
		    double **signal, double **gx, double **gy, char **msm);
void genera_medida_UPM( int dimx, int dimy, double **signal, double **disp);
void genera_medida_UPM_corr( int dimx, int dimy, double **signal, 
	double **gx, double **gy, double **disp);
void genera_medida_UPM_rel( int dimx, int dimy, double **signal, double **gx,
	double **gy, double **disp);
void genera_medida_UPM_glob( int dimx, int dimy, double **signal, 
	double **gx, double **gy, double **disp);
double mascara( int dimx, int dimy, double umbral, double **disp, 
	char **msm);

void derive_cross( double *Rc, double *Ic, double *Rgx, double *Igx, 
	double *Rgy, double *Igy, double ***Fourier);
void reconstruct_cross(double *Rgx, double *Igx, double *Rgy, double *Igy,
	double *Rc, double *Ic, double ***Fourier);
void Fourier_cross( double ***Fourier, double *Rf, double *If, int sign);
void singularize_sign(int dimx, int dimy, int xeff, int yeff, int modo, char **msm, 
	double **gx, double **gy);
void singularize(int dimx, int dimy, int xeff, int yeff, char **msm, 
		 double **gx, double **gy);

void reconstruct(int dimx, int dimy, double **gxR, double **gyR);
void genera_sources(int dimx, int dimy, int tipo, int sign, 
	double **gxR, double  **gyR);
void derive_complex(int dimx, int dimy, int tipo, int sign, 
	double **gxR, double  **gyR);

double compute_norma_filtro(int dimx, int dimy, double expon);

void kernelea( int dimx, int dimy, int xeff, int yeff, char **msm, double **signal);
void deskernelea( int dimx, int dimy, int xeff, int yeff, char **msm, double **signal);


void compute_sources( int dimx, int dimy, int xeff, int yeff, double **dummy,
	double **gx, double **gy, double *meang);
void recupera_dummy( int dimx, int dimy,  int xeff, int yeff, char **msm, double **dummy);
void recupera_gradiente_esencial( int dimx, int dimy, int xeff, int yeff, double **dummy,
	double **gx, double **gy, double *meang);
void process_sources( int dimx, int dimy, int xeff, int yeff, double **gx, double **gy);

int haz_histograma( int dimx, int dimy, char *nombre, double **datos);

void represent_sources( int dimx, int dimy, int dimv, int dimz, int iz, char *ext, 
	double ***gx, double ***gy);
void represent_sources_mod( int dimx, int dimy, char *ext, int silog,
	double **mux, double **muy);
void prepare_sources_mod( int dimx, int dimy, int silog, double **mux, double **muy, 
			  double **mod);
void represent_sources_fase( int dimx, int dimy, char *ext, 
	double **mux, double **muy);
void prepare_sources_fase_block( int dimx, int dimy, double block, double **mux, 
	double **muy, char **Red, char **Green, char **Blue);
void phase_color( double mx, double my, char *Red, char *Green, char *Blue);

void pinta(int dimx, int dimy, char **msm, double **data);
void pinta_suave( int dimx, int dimy, int xeff, int yeff, char **msm, double **data);
void deriva_limits( int dimx, int dimy, int xeff, int yeff, char **msm, 
		    double **ax, double **ay);
void deriva_orient_limits( int dimx, int dimy, int xeff, int yeff, char **msm, 
			    double **ax, double **ay);


#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
