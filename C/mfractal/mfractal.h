/* ===================================
** mfractal.h
** started on Fri Jan 19 17:49:29 2007 
** ===================================
*/

#ifndef   	_MFRACTAL_H_
#define   	_MFRACTAL_H_

/** Useful variables **/

#ifndef SIGCONFIDENCE
#define SIGCONFIDENCE 3. // 3 sigmas corresponds to a 99.7% confidence level
#endif

#ifndef RAND_PER 
#define RAND_PER 2147483647 // Substitute a RAND_MAX for the random routine;
#endif

#ifndef MINEV
#define MINEV 60. // 30.
#endif

/** Variables dedicated to fractal simulation **/

/* You can modify the variables below 
 * Note: some pseudo constants have a "DEF" prefix because they are
 * susceptible to be modified in the code (see the end of this file) */

/* Characteristics of the simulated data */
#ifndef DSPACE  // Dimension of the embedding space
#define DSPACE DIM1D  
#endif
#ifndef NSERIES  // Number of output series
#define NSERIES 1  
#endif
#ifndef LMAX  // Size of the generated series 
#define LMAX 512    
#endif
#ifndef OUTRES  // Number of resolutions to be eliminated to smoothen
#define OUTRES 0   
#endif

/* Defining the type of fractal */
#ifndef TYPLOGPOISSON
#define TYPLOGPOISSON 0
#endif
#ifndef TYPLOGNORMAL
#define TYPLOGNORMAL 1
#endif
#ifndef TYPLOGLEVI
#define TYPLOGLEVI 2
#endif
#ifndef TYPBINOMIAL
#define TYPBINOMIAL 3
#endif
#ifndef TYPMONOFRACTAL
#define TYPMONOFRACTAL 4
#endif
#ifndef TYPE_MFSIM /* Type of multifractal to be generated:
		    * 0: Log-Poisson   2: Log-Levi   4: Monofractal
		    * 1: Log-Normal    3: Binomial    */
#define TYPE_MFSIM TYPLOGPOISSON
#endif

#ifndef DEFFLAG_INVTRANS  /* Default initial flag for explicit translational
			   * invariance. See pseudo-constant FLAG_INVTRANS at
			   * the end of this file */
#define DEFFLAG_INVTRANS FALSE 
#endif

/* Defining the parameters of the genrated fractals */
#ifndef DEFHINF  /* Default initial Most singular exponent (Types Log-Poisson,
		  * Binomial, and Monofractal). See pseudo-constant HINF at the
		  * end of this file */
#define DEFHINF  -0.5
#endif
#ifndef CODINF  // Most singular co-dimension (Types 0 and 3)
#define CODINF 1.
#endif
#ifndef DH  // Accepted conventional error
#define DH 0.1
#endif

#ifndef H0  // Lower singularity bound in binomial model
#define H0   -1.
#endif
#ifndef DEFH1  /* Default initial upper singularity bound in binomial 
		* model. See pseudo-constant H1 at the end of this file */
#define DEFH1   0.5 // should be >= DEFHINF (typically: DEFHINF+DH)
#endif

#ifndef DEFMU  /* Default initial average singularity (for log-Normal and
		* log-Levi). See pseudo-constant MU at the end of this file */
#define DEFMU  0.5
#endif
#ifndef SIGMA  // Singularity dispersion in log-Normal and log-Levi
#define SIGMA 1.  
#endif
#ifndef DEFALPHA  /* Default initial exponent defining log-Levi.
		   * See pseudo-constant ALPHA at the end of this file */
#define DEFALPHA  1.5  
#endif 

#ifndef DENS  // Density (monofractals only)
#define DENS 0.1   
#endif

#ifndef DISPETA     // Log-Normal processing variables
#define DISPETA 0.5
#endif

#ifndef TCH  // Maximum effective dispersion (in sigmas)
#define TCH  3.    
#endif
#ifndef WAVBASE   /* Choice of the wavelet basis for representation:
		    * -1: Gaussian   0: Haar   1: Lorentian */
#define WAVBASE
#endif
#ifndef DERWAVBASE  // Derivative order of the wavelet of choice
#define DERWAVBASE 2.
#endif

#ifndef NBOX  // Number of bins in the histogram
#define NBOX 256 // 128 // 400
#endif

#ifndef FLAG_ANALYSIS  /* Enables direct obtention of singularity spectra from the
		       * currently generated data by two histogram methods */
#define FLAG_ANALYSIS FALSE
#endif


/** Variables for fractal analysis */

#ifndef TYPMOM
#define TYPMOM 0
#endif
#ifndef TYPGH
#define TYPGH 1
#endif
#ifndef TYPGMWP
#define TYPGMWP 2
#endif
#ifndef TYPWTMM
#define TYPWTMM 3
#endif

#ifndef TYPE_MFANA
#define TYPE_MFANA TYPGMWP
#endif

#ifndef NPOINTS_SPECTRUM // Number of points to be solved in the spectrum
#define NPOINTS_SPECTRUM 20
#endif

#ifndef NPOINTS  /* Number of scales in the multifractal evaluation
		  * -1 + number of points in the log-log linear regression 
		  * used to obtain the singularity exponents at every point */
#define NPOINTS 6
#endif
#ifndef FLAG_HOLDER // Analizes de measure (FALSE) or the own function
#define FLAG_HOLDER FALSE
#endif

#ifndef FLAG_SAVEDOUBLE  /* Flag for saving the data series in double format
		      * Default: DISABLE flag, the data are saved in double
		      * format */
#define FLAG_SAVEDOUBLE FALSE
#endif
#ifndef FLAG_SAVEFLOAT  // Ibid in float format
#define FLAG_SAVEFLOAT FALSE
#endif

/** Variables for MSM analysis **/

/* Boundary quantils used in the determination of the MSM */
#ifndef MSMMIN
#define MSMMIN 0.01
#endif
#ifndef MSMMAX
#define MSMMAX 0.05
#endif

/* Boundaries for possible exponent in the histograms */
#ifndef HMIN
#define HMIN -1.
#endif
#ifndef HMAX
#define HMAX 2.
#endif
/*
  #ifndef HMIN
  #define HMIN H0
  #endif
  #ifndef HMAX
  #define HMAX H1
  #endif
*/

/* Variables for UPM computation */
#ifndef MODUPM 
#define MODUPM 0
#endif
#ifndef MODGLOBUPM
#define MODGLOBUPM 1
#endif
#ifndef MODRELUPM
#define MODRELUPM 2
#endif
#ifndef MODCORRUPM
#define MODCORRUPM 3
#endif
#ifndef MODE_UPM
#define MODE_UPM MODRELUPM
#endif

/* Variables specific to the computation of sources */
#ifndef UPM_THR
#define UPM_THR 1.
#endif
#ifndef UPM_DENS /* if !=0, the UPM is calculated at such density */
#define UPM_DENS 0.
#endif
#ifndef EXP_MU
#define EXP_MU 1.
#endif

/* Variable for representation */
#ifndef FLAG_MSMCOMP
#define FLAG_MSMCOMP FALSE
#endif

/* used for BW represetation of MSM and other fractal sets */
#ifndef REVERSE
#define  REVERSE 1   /* if reverse=0 the background is black and 
		      * the mask is white; if reverse=1, the opposite. */
#endif
/* To represent Black/White masks; if REVERSE=0 the background is black
 * and the mask is white; if REVERSE=1, the opposite. */
#ifndef C0  // grey level for 0
#define C0 (char)0xff // (char)(REVERSE*255)
#endif
#ifndef CP  // grey level for 1 
#define CP (char)0x00 // (char)((1-REVERSE)*255)
#endif
#ifndef CM  // medium gray level, for indefinite cases
#define CM (char)0x7f       
#endif


/** Variables for WTMM analysis */

#ifndef WAVWTMM 
#define WAVWTMM WAVGAUSS /* Default choice for wavelet: gaussian */
#endif
/* Note: WAVGAUSS is defined in filter.h */
#ifndef ORDDERWTMM
#define ORDDERWTMM 2.
#endif

#ifndef DIRECT
#define DIRECT 1
#endif
#ifndef CANONICAL
#define CANONICAL (1-DIRECT) /* If selected, the calculation of WTMM exponents
			      * is done by canonical Legendre transform */
#endif
#ifndef FLAG_METHODWTMM
#define FLAG_METHODWTMM DIRECT
#endif

#ifndef FLAG_SUPWTMM
#define FLAG_SUPWTMM FALSE // Flag. If given, WTMM is evaluated over suprema
#endif
 
#ifndef MINMOMENT
#define MINMOMENT -5.
#endif
#ifndef MAXMOMENT
#define MAXMOMENT 5.
#endif
#ifndef DQ
#define DQ 0.2         
#endif

#ifndef SCMIN
#define SCMIN 2.
#endif
#ifndef SCMAX
#define SCMAX 100.
#endif

#ifndef SHIFTSPEC
#define SHIFTSPEC 0.
#endif

/** Don't modified the lines below !!! 
 ** See function prepara_multifractal for explaination **/

#ifndef FLAG_INVTRANS  // Enables implement explicit translational invariance
#define FLAG_INVTRANS ( /*if*/ \
		    ((TYPE_MF==TYPLOGLEVI) && /*and*/ \
		     (DEFFLAG_INVTRANS==TRUE) && /*and*/ \
		     (ALPHA==1.)) ? /*then it equals**/ \
		    FALSE : /*otherwise*/ \
		    DEFFLAG_INVTRANS /*default, all other cases*/ \
		     )
#endif

#ifndef ALPHA  // Exponent defining log-Levi
#define ALPHA (double)MAX( DEFALPHA, 0.01 ) 
#endif

#ifndef MU  // Average singularity (for log-Normal and log-Levi)
#define MU ( /*if*/ \
	     ((TYPE_MF==TYPLOGNORMAL) && /*and*/ \
	      (FLAG_INVTRANS==TRUE)) ? /*then it equals*/ \
	     (SIGMA*SIGMA/4.) : \
	     ( /*otherwise if*/ \
	       ((TYPE_MF==TYPLOGLEVI) && /*and*/ \
		(FLAG_INVTRANS==TRUE) && /*and*/ \
		(ALPHA!=1.)) ? /*then it equals*/ \
	       (pow(SIGMA/ALPHA,1.+1./(ALPHA-1.))*(ALPHA-1.)) : /*else*/ \
	       DEFMU /*default, all other cases*/ \
	       ) \
	     ) 
#endif
     
#ifndef H1  // Upper singularity bound in binomial model
#define H1  ( /*if*/ \
	      (TYPE_MFSIM==TYPBINOMIAL) ? /*then*/ \
	      ( /*if*/ \
		(FLAG_INVTRANS==TRUE) ? /*then*/ \
		( /*if*/ \
		  (DEFHINF>1.) ? /*then it equals*/ \
		  (-log(1-pow(0.5,DEFHINF))) : \
		  ( /*otherwise if*/ \
		    ((DEFH1>0.) && /*and*/ \
		     (DEFH1<1.)) ? /*then it equals*/ \
		    0.75 : /*else*/ \
		    DEFH1 \
		    ) \
		  ) : /*otherwise*/ \
		( /* if*/ \
		  (DEFH1<=DEFHINF) ? /*then it equals*/ \
		  (DEFHINF+DH) : /*else*/ \
		  DEFH1 \
		  ) \
	       ) : \
	      DEFH1  /*default, all other cases*/ \
	      )
#endif

#ifndef HINF  // Most singular exponent (Types 0, 3 and 4)
#define HINF  ( /*if*/ \
		(TYPE_MFSIM==TYPBINOMIAL) && /*and*/ \
		(FLAG_INVTRANS==TRUE) && /*then*/ \
		(DEFHINF<1.) ? /*then it equals*/ \
		(-log(1-pow(.5,H1))) : /*else*/ \
		DEFHINF /*default, all other cases*/ \
		)
#endif
     
#endif 	    /* !_MFRACTAL_H_ */

