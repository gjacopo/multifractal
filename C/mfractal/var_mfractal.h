/* ===================================
** var_fractal.h
** started on Tue Jan 23 15:14:18 2007 
** ===================================
*/

#ifndef   	_VAR_FRACTAL_H_
#define   	_VAR_FRACTAL_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

  const int nqDef = 65;
  const int qDefArray0={ -4.,-3.6,-3.2,-3.,-2.8,-2.6,-2.4,-2.2,
			 -2.,-1.8,-1.6,-1.4,-1.2,-1.1,
			 -1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,
			 0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,
			 0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,
			 1.,1.05,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
			 2.,2.3,2.6,3.,3.5,4.,5.,6.,7.,8. };

  /* ...and a default initialisation of it */  
  ParFRAC *p0_fil = {
    /** Variables related to WTMM analysis */
    DIRECT,               //  : method_wtmm
    WAVWTMM,             //  : wav_wtmm
    ORDDERWTMM,         //  : ord_der_wtmm
    MINSCALE, MAXSCALE,   //  : minscale, maxscale
    FALSE                 //  : flag_tauq_wtmm
    SHIFTSPEC             //  : shift
    NVOICES, SCTIME,      //  : nvoices, sctime
    SCRATIO               //  : scratio
    SCMIN, SCMAX,         //  : scmin, scmax
    MINMOMENT, MAXMOMENT, //  : minmom, maxmom
    DQ,                   //  : qstep
    65,
    qDefArray0            // qDefArray
  };

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !_VAR_FRACTAL_H_ */

#ifndef   	VAR_MSM_H_
#define   	VAR_MSM_H_

/***      Global variables               ***/

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

/*** Variables concerning the image ***/

/* If XMAX or YMAX are equal to zero it reads inrimage parameter */
int XMAX=0;  /* X size of processed window */
int YMAX=0;  /* Y size of processed window */
int X0=0;
int Y0=0; /* up left corner coordinates of the processed window */

int FOTO=0; /* If we are dealing with a photo or an inrimage */
int COLOR=0; /* In the first case, if it is a color image or not */

double BLOCK=1.; /* input resolution block size */
double BLOUT=1.; /* output resolution block size */


/*** Variables specific to the computation of manifolds ***/

int NBOX=400; /* number of boxes in the histograms */

int MULTIFRACTAL=0; /* If actual multifractal processing is performed
		     * or not */

int NPUNTOS=6;   /* -1 + number of points in the log-log linea regression 
		 * used to obtain the singularity exponents at every
		 * point */


#define BUEN 0.9  /* criterion on the linear regression coefficient */
int NSING=3;  /* number of manifolds to be represented */

double HINF=-0.5;
double DH=0.2; 	/* MSM parameters by default */

/* Graphic representation parameters */
#define  REVERSE 1   /* if reverse=0 the background is black and 
		      * the mask is white; if reverse=1, the opposite. */
char C0=(char)(REVERSE*255); /* grey level associated to 0 */
char CP=(char)((1-REVERSE)*255); /* grey level associated to 1  */
/* To represent Black/White masks; if reverse=0 the background is black
 * and the mask is white; if reverse=1, the opposite. */
char CM=0x7f;       /* medium gray level, for indefinite cases */


/*** Variables specific to the computation of sources ***/

double UPM_THR=1.;
double UPM_DENS=0.;
double EXP_MU=1.;


/*** Numerical parameters ***/

int CUT=0; /* to avoid recurrent reprocessing:
	    *  - 0 means full processing
	    *  - 1 the process computes and records the sources
	    *  - 2 the process starts with the recorded sources */

int OPTIMA=0; /* Option for (sligthly) optimizing the code */

int STDDEVCHROM=0; /* Computes a disparity measure between the original image 
		    * and the chromatically reduced one, which is the local   
		    * standard deviation of the difference of both signals. 
		    * It's an alternative solution to the computation of the
		    * sources. It provides a measure of the deviation between 
		    * both gradient fields. */ 
int GRAD=0; /* Provides the gradients that will be used for computing the MSM,
	     * the sources  and so on, instead of the spatial gradients of the 
	     * original image */
int DEVIATION=0; /* */
int DECONVOLUTION=0; /* Operates a deconvolution of an image: 
		      *  - the input image is supposed to be noisy,
		      *  - the entropy is computed over this image,
		      *  - the MSM is normally computed,
		      *  - but the reconstruction is performed with a 
		      *    reconstructing set composed with pixels of the MSM 
		      *    satisfying a given criterium.     */
int COMPLEMENT=0; /* Parameter for computing the reconstructed image from 
		   * the complementary set of the MSM 
		   * cas 0: cas normal: on calcule la MSM et on reconstruit 
		   * a partir de la MSM
		   * cas 1: on calcule l'ensemble complementaire de la MSM et 
		   * on effectue la reconstruction a partir de cet ensemble */
int UNITARY=0;  /* If the reduced image is assummed to be 
		 * the unitary image (1) or not (0) */

int TRUNCATE=0;


/*** Variables specific to entropy-based sources ***/

double ENTR_THR=0.75;
int ENTROPY=0; /* Computes the multiscale windowed entropy and compute the 
		* sources associated with this field */
int WINSIZE=10;  /* dimension de la fenetre analysante */
int PER=1; /* cyclicité des bords */
int GLOBAL=1; /* Normalisation globale */
int WEIGHT=-2; /* facteur de pondération multiéchelle */


/*** Graphic representation variables ***/

int ANGULO=0; /* Parameter for the representation of the orientation of 
	       * the different fields computed in the program 
	       * If 0, represent orientation of the vectorial field of sources 
	       * with a color image.
	       * If 0, represent values of the orientation of the sources */
int VIDEO=0; /* If the resulting images are presented or not as animated GIF's */
int VISUALIZE=1; /* flag to decide whether visualize the manifolds or
		    not */
int ORIGVISU=0; /* Save original normalized image */
int SAVE=0; /* Save the directionnal components of the source fields */

int VERBOSE=0;


#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !VAR_MSM_H_ */
