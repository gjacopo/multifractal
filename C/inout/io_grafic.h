#ifndef _IO_GRAFIC_H_
#define _IO_GRAFIC_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
    
  int read_dimensiones_foto( char *nombre, int *dimx, int *dimy );
  int read_dimensiones_ppm( char *nombre, int *dimx, int *dimy );
  int read_dimensiones_gif( char *nombre, int *dimx, int *dimy );
  
  int ext_grafic( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		  char *nombre, Read2D *p_lee);
  int ext_gif( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	       char *nombre, Read2D *p_lee);
  int ext_ppm( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
	       char *nombre, Read2D *p_lee);
  
  int check_grafic( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		    char *nombre, Read2D *p_lee);
  int check_gif( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		 char *nombre, Read2D *p_lee);
  int check_ppm( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		 char *nombre, Read2D *p_lee);
  
  int read_foto_gris( int dimx, int dimy, char *nombre, double **data );
  int read_foto_color( int dimx, int dimy, char *nombre, char **Red, 
		       char **Green, char **Blue );
  int read_color_block( int dimx, int dimy, int block, char *nombre,
			double ***data, double *med_cr );
  
  int read_ppm_gris( int dimx, int dimy, char *nombre, double **data );
  int read_ppm_rgb( int dimx, int dimy, char* nombre_in, 
		    char **Red, char **Green,	char **Blue );
  int read_ppm_raw( int dimx, int dimy, char* nombre, double **data);
  int read_ppm( int dimx, int dimy, int dimv, int dimz, int iz, int bd,
		char *nombre, double ***datos);
  int read_en_ascii( FILE* canal );

  int read_pgm( int dimx, int dimy, char* nombre_in, double **cont );

  int read_gif_gris( int dimx, int dimy, char *nombre, double **data );
  int read_gif_rgb( int dimx, int dimy, int nimag, char *nombre, char **Red,
		    char **Green, char **Blue );
  int read_gif( int dimx, int dimy, int dimv, int dimz, int iz, int bd,
		char *nombre, double ***data);

  void interpreta_extension_gif( FILE *canal, int verbose );

  int read_graphic_ce_gif( FILE *canal, int verbose );
  int read_comment_gif( FILE *canal, int verbose );
  int read_application_gif(FILE *canal, int verbose );
  int interpreta_stream_gif( FILE *canal, unsigned char **encoded );
  int write_stream_gif( int l_code, unsigned char *encoded, FILE *canal );
  int annade_decoded_gif(int w0, int letra, int *tab0, int *tab1, int id0, 
			 int *decoded );
  int primer_car_gif( int w0, int letra, int *tab0 );
  
  int write_gif_rgb( int dimx, int dimy, char *nombre,
		 char **Red, char **Green, char **Blue );
  int write_gif_rgb_animado( int dimx, int dimy, int modo, int colort, 
			 char *nombre, 
			 char **Red, char **Green, char **Blue );
  
  int create_table_color_gif_rgb_icy( int dimx, int dimy, char **Red, char **Green, 
				    char **Blue, int *dimI, int *dimC, int *dimY,
				    int **table_ICY );
  int create_table_color_gif( int dimx, int dimy, 
			    char **Red, char **Green, char **Blue, 
			    int **table_c, int **freq_c );

  int RGB_ICY( int dimI, int dimC, int dimY, int mode, char *Red, char *Green,
	       char *Blue, int *iI, int *iC, int *iY );
  int reduccion_cromatica_gif(int dimx, int dimy, int factor, 
			       char **Red, char **Green, char ** Blue );
  void codify_imagen_con_table_ICY( int dimx, int dimy, int dimI, int dimC,
				      int dimY, char **Red, char **Green, 
				      char **Blue, int color_size, 
				      int *table_ICY, int *decoded );
  void codify_imagen_con_table_pos_gif( int dimx, int dimy,
					  char **Red, char **Green, char **Blue,
					  int color_size, int *table_c, 
					  int *freq_c, int *decoded );
  void codify_imagen_con_table_color_gif( int dimx, int dimy,
					    char **Red, char **Green, char **Blue,
					    int color_size, 
					    char *Rg, char *Gg, char *Bg, 
					    int *freq_c, int *decoded );
  void codify_imagen_gris_gif( int dimx, int dimy, 
				 char **Red, char **Green, char **Blue, 
				 int *decoded );
  int table_pos_table_color_gif( int color_size, int *table_c,
				  char **Rg, char **Gg, char **Bg, 
				  int dimx, int dimy,
				  char **Red, char **Green, char **Blue );
  
  
  int write_mask_block( int dimx, int dimy, double block, char* nombre,
			 char **mask );
  int write_foto_block( int dimx, int dimy, double block, char* nombre,
			 double **datos );
  int write_foto_block_limites( int dimx, int dimy, char* nombre,
				 double min_c, double max_c, double **datos );
  int write_foto_vec_block( int dimx, int dimy, char* nombre,
			     double **vx, double **vy);
  int write_foto_4( int dimx, int dimy, char* nombre, double **datos );
  int write_binary_block( int dimx, int dimy, double block, char* nombre,
			  char **mask );
  int write_binary_foto( int dimx, int dimy, int dimz, int iz, 
			 char *base, char *ext, char **bin );
  int write_RGB_block( int dimx, int dimy, char* nombre,
		       char **Red, char **Green, char **Blue );
  int write_foto_color_block( int dimx, int dimy, int n_cr, double block,
			       char *nombre, double ***datos );
  int write_video_block( int dimx, int dimy, int n_cr, double block, 
			  int modo, char *nombre, double ***datos );
  int write_video_RGB_block( int dimx, int dimy, double block, 
			      int modo, char *nombre, 
			      char **Red, char **Green, char **Blue );
  
  
  int prepare_foto_block( int dimx, int dimy, double block, double **datos,
			   char **foto );
  int prepare_foto_fijo_block( int dimx, int dimy, int levels,
			       double **datos, char **foto );
  int prepare_foto_block_limites( int dimx, int dimy, double min_c, 
				   double max_c, double block, 
				   double **datos, char **foto );
  int prepare_foto_log( int dimx, int dimy, double **datos, char **foto );
  int prepare_foto_4( int dimx, int dimy, double **datos, char **foto );
  int prepare_foto_log_4( int dimx, int dimy, double **datos, char **foto );
  int prepare_char_block( int dimx, int dimy, double block, char **grayin,
			  char **grayout );
  int write_pgm( int dimx, int dimy, int bin, char* nombre_out, char **cont );
  int write_pgm_foto( int dimx, int dimy, int dimz, int iz, 
		      int bin, char *base, char *ext, char **cont );
  int write_ppm( int dimx, int dimy, int bin, char *nombre_out, char **Red, 
		 char **Green, char **Blue );
  
  int write_gray_foto( int dimx, int dimy, int dimz, int dimv, int iz, 
		       int bin, char *base, char *ext, char ***cont );
  int write_foto_RGB( int dimx, int dimy, int dimv, int dimz, int iz, 
		      char *base, char *ext, char ***Red, char ***Green, char ***Blue );
  int write_foto( int dimx, int dimy, int dimv, int dimz, int iz, 
		  char *base, char *ext, double ***data );
  int write_foto_color( int dimx, int dimy, int dimz, int iz, 
			char *base, char *ext, double ***data );
  
  void paleta( double rep, char *Red, char *Green, char *Blue );
  
  int write_image_reference( int dimx, int dimy, int dimv, int dimz, int dimt,
			     int iz, char *base, char **nombre_tipo, char *ext,
			     double ***data );
  int write_image_foto( int leff, int dimy, char *nombre_in, double **datos);

  int visualise_color( int dimx, int dimy, char *name_in, 
		       double m, double M, double **data );
  int visualise_gris( int dimx, int dimy, char *name_in, 
		      double m, double M, double **data );
  int visualize_gris_pieces( int dimx, int dimy, char *dext,
			     double m, double M, double delta_h, double **data );
  
#ifdef __cplusplus
}
#endif	/* __cplusplus */

#endif
