/* ===================================
** io_parse.h
** ===================================
*/
#ifndef   	_IO_PARSE_H_
#define   	_IO_PARSE_H_

#ifndef _PARSE_IO_PARAMETERS_
#define _PARSE_IO_PARAMETERS_
#endif

/** Definition of the structures for parsing the parameters
 ** for reading input/writing output data  **/

#ifndef _PARIO_TYPE_
#define _PARIO_TYPE_
typedef struct pario{
  int dim_space;
  char in[MAXNAMELENGTH];
  char out[MAXNAMELENGTH], ext[MAXNAMELENGTH];
  int flag_window;
  int x0, y0;
  int x, y;
  int flag_bin;
  int flag_res;
  double block_in, block_out;
  int flag_foto, flag_inr, flag_color;
  int flag_video, flag_visu;
} ParIO;
#endif

/** Declaration of the prototypes of the functions parsing the
 ** input/ouput parameters  **/

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
#ifdef _PARIO_TYPE_
  ParIO* ioparse_alloc();
  int ioparse_free( ParIO *par );
  int ioparse_default( ParIO *par );
  ParIO* ioparse_create( );
  int ioparse_init( int in, int siflag, int *deflag, 
		    char **olarg, char **olval, char **olexp,
		    double **ptrvar, double **ptrval, 
		    int **iptrvar, int **iptrval, char **cptrvar,
		    int **ptrflag, int *type, int*olnumb );
  int ioparse_update( );
  int ioparse_select( const int choice[] );
  
#ifdef Debug
  int ioparse_display( );
  int ioparse_test (int argc, char *argv[]);
#endif
#endif
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

/** Useful constants when parsing arguments **/

#ifndef NARG
#define NARG 0
#endif

enum { ipar_dspace=NARG, ipar_in, ipar_out, ipar_ext,
       ipar_dim0, ipar_dim, ipar_res, 
       ipar_fot, ipar_inr, ipar_vis, ipar_col, ipar_vid, 
       ipar_ver, 
       io_npar };

#ifdef IO_NPAR
#undef IO_NPAR
#endif
#define IO_NPAR io_npar;

#endif 	    /* !_IO_PARSE_H_ */
