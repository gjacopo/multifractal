/* ===================================
** mf_parse.h
** started on Thu Jan 18 18:06:02 2007 
** ===================================
*/

#ifndef   	_MF_PARSE_H_
#define   	_MF_PARSE_H_

#ifndef _PARSE_FRACTAL_PARAMETERS_
#define _PARSE_FRACTAL_PARAMETERS_
#endif

/** Definition of the structures for storing the different parameters
 ** of the fractal simulation and/or analysis 
 **/
#ifndef _PARFRAC_TYPE_
#define _PARFRAC_TYPE_
typedef struct parfrac {
  /** Variables related to multifractal simulation */
  int dim_space;
  int type_mf; /* Type of multifractal to be calculated:
		  * 0: Log-Poisson   2: Log-Levi   4: Monofractal
                  * 1: Log-Normal    3: Binomial    */
  int ndata; // Number of in/output data
  int lmax;  // Size of the generated series 
  int out_res; // Number of resolutions to be eliminated to smoothen
  int flag_invtrans; // Enables implement explicit translational invariance
  double hinf, // Most singular exponent (Types 0, 3 and 4)
    codinf;    // Most singular co-dimension (Types 0 and 3)
  double h0, // Lower singularity bound in binomial model
    h1;   // Upper singularity bound in binomial model
  double dh;  // Accepted conventional error
  double mu,  // Average singularity (for log-Normal and log-Levi)
    sigma,    // Singularity dispersion in log-Normal and log-Levi
    alpha,    // Exponent defining log-Levi
    dens;     // Density (monofractals only)
  double dispeta; // Log-Normal processing variables
  double tch; // Maximum effective dispersion (in sigmas)
  int nbox;  // Number of bins in the histogram
  int s0;

  /** Variables related to multifractal analysis **/
  int npoints;/* Number of scales in the multifractal evaluation
	       * -1 + number of points in the log-log linear regression 
	       * used to obtain the singularity exponents at every point */
  int flag_holder; // Analize of the measure (FALSE) or the own function
  int flag_savedouble; /* Flag for saving the data series in double format
			 * Default: DISABLE flag, the data are saved in double
			 * format */
  int flag_savefloat; // Ibid in float format
  int nsing;
  double msm_min, msm_max; /* boundary quantils used in the determination 
			    * of the most singular exponent */
  /* Specific to 2D analysis */
  int flag_1Danalysis // If set, the analysis is performed in 1D
  double thetau; // Direction for the 1D analysis
} ParFRAC;
#endif
  
  /** Variables related to MSM analysis */
#ifndef _PARMSM_TYPE_
#define _PARMSM_TYPE_
typedef struct parmsm {
  int mode_upm;
  int flag_msmcomp;
  double upm_thr, upm_dens; 
  double exp_mu;
  int nbox;  // Number of bins in the histogram
} ParMSM;
#endif

  /** Variables related to WTMM analysis */
#ifndef _PARWTMM_TYPE_
#define _PARWTMM_TYPE_
typedef struct parwtmm {
  int flag_methodwtmm;
  int flag_supwtmm;
  int flag_tauqwtmm;
  double shift;
  double scmin, scmax;
  double minmom, maxmom;
  double qstep;    
  int nqDef;
  double *qDefArray;
} ParWTMM;
#endif


#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */


#ifdef __cplusplus
}
#endif		/* __cplusplus */

#ifndef NARG
#define NARG 0
#endif

enum { ipar_=NARG, ipar_, ipar_, ipar_,
       ipar_, ipar_, ipar_, ipar_, 
       ipar_, 
       ipar_, ipar_, ipar_, frac_npar}

#ifdef FRAC_NPAR
#undef FRAC_NPAR
#endif
#define FRAC_NPAR frac_npar

#endif 	    /* !_MF_PARSE_H_ */

