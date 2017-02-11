/* ===================================
** flt_parse.h
** started on Thu Jan 18 18:06:02 2007 
** ===================================
*/

#ifndef   	_FLT_PARSE_H_
#define   	_FLT_PARSE_H_

#ifndef _PARSE_FILTER_PARAMETERS_
#define _PARSE_FILTER_PARAMETERS_
#endif

/** Definition of the structures for storing the parameters
 ** of the different available filters 
 **/

#ifndef _PARFILT_TYPE_
#define _PARFILT_TYPE_
typedef struct parfil {
  /** Variables related to FFT analysis */
  int flag_memory; /* Exchanges memory efficiency with speed, by using FFT or FFFT
		    * Defaults to the use of FFT (consumes more memory, it 
		    * requires the use of dimensions being power of 2 but it
		    * is faster) with respect to FFFT */
  /** Variable related to derivative filters */
  int mode_deriva;  /* Derivation mode.
		     * 0: Half-pixel by Fourier transform. Generates aliasing in
		     *    wide areas of low gradient
		     * 1: Forward 1-pixel increment. Simple and less sensitive
		     *    to artifacts, but reconstruction from it is of less
		     *    quality (unknown reason) */
  int mode_freq; /* Choose the definition of the frequency vector : 
		  * if TRUE: f_x=sin(PI*x), f_y=sin(PI*y)
		  * else: f_x=x,f_y=y
		  * with (x,y) coordinates in Fourier space */
} ParFILT;
#endif

/** Variables related to WT analysis */
#ifndef _PARWAV_TYPE_
#define _PARWAV_TYPE_
typedef struct parwav {
  int wav;
  int ord_der;
  double wav_range;
  double thetau;
  double minscale, maxscale;
  double scale0;/* Minimum scale */
  int nvoices, sctime;
  double scratio;
} ParWAV;
#endif


/** Declaration of the prototypes of the functions parsing the
 ** different parameters of the filters
 **/

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  

#ifdef _PARFILT_TYPE_
  ParFILT* filparse_alloc();
  int filparse_free( ParFILT *par );
  ParFILT* filparse_create( );
  int filparse_default( ParFILT *par );
#endif
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

/** Useful constants when parsing arguments **/

#ifndef NARG
#define NARG 0
#endif

enum { ipar_mem=NARG, ipar_dermode, ipar_freq,
       ipar_wav, ipar_der, ipar_theta, ipar_scrange, 
       ipar_momq, ipar_shift, ipar_tauq, 
       filt_npar };

#ifdef FILT_NPAR
#undef FILT_NPAR
#endif
#define FILT_NPAR filt_npar

#endif 	    /* !_FLT_PARSE_H_ */

