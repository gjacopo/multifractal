#ifndef _IO_RAW_H_
#define _IO_RAW_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
    
  int read_winraw( char *nombre, 
		   int bd, int dimx, int dimy, int dimv, int dimz, int iz, 
		   int ix0, int iy0, int ndimx, int ndimy, 
		   double ***datos);
  int read_raw( char *nombre, 
		int bd, int dimx, int dimy, int dimv, int dimz, int iz, 
		double ***data);
  int read_dimensiones_raw( char *name,
			    int *bd, int *dimx, int *dimy, int *dimv, int *dimz );
  int read_flow( int dimx, int dimy, char* nombre_out, double **data );
  int read_unformatted( int dimx, int dimy, int dimv, int dimz, int iz, int bd,
			char *nombre, double ***data);

  int read_en_double_int( int xmax, int ymax, int ix0, int iy0, 
			  int dimx, int dimy,
			  int littleindian, char* nombre, double **data );
  
  int ext_unformatted( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
		       char *nombre, Read2D *p_lee);
  
  int check_unformatted( int *dimx, int *dimy, int *dimv, int *dimz, int *bd,
			 char *nombre, Read2D *p_lee);
  
  int read_data( int leff, int dimy, char *nombre_in, double **signal);
  int write_data( int leff, int dimy, char *nombre_in, double **signal);
  int read_data_infloat( int leff, int dimy, char *nombre_in, double **signal);
  int write_data_infloat( int leff, int dimy, char *nombre_in, double **signal);
  int write_serie( int leff, char *nombre_in, double *datos);

  
#ifdef __cplusplus
}
#endif	/* __cplusplus */

#endif
