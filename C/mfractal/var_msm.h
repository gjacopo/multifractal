/********************************************************/
/*                                                      */
/*                  variables_msm.h                     */
/*           Version del 2 de Dekembriou, 2004          */
/*                                                      */
/*  Global variables used in the MSM analyzis scheme.   */
/*                                                      */
/********************************************************/

/*      Global variables               */

#ifndef   	VAR_MSM_H_
# define   	VAR_MSM_H_

/** Define the FLAG_MSM variable **/
#ifndef FLAG_MSM
#define FLAG_MSM
#endif

/** Variables associated to the file multifractal_generator.c */

//int NSERIES=1;     // Number of output series
//int D_space=1;     // Dimension of the space
//int LMAX=512;      // Size of the generated series
//int OUTRES=0;      // Number of resolutions to be elliminated to smoothen. 
//int TYPE_MF=0;        // Type of multifractal to be calculated.
                   // 0: Log-Poisson
                   // 1: Log-Normal
                   // 2: Log-Levi
                   // 3: Binomial
                   // 4: Monofractal


//int FLAG_INVTRANS=0;    Enables implement explicit translational invariance
//double H0=-1;      // Upper singularity bound in binomial model
//double DENS=0.1;    // Density (monofractals only)

//double TCH=3.;      // Maximum effective dispersion (in sigmas)

//int WAV_BASE=0;    
//int DER_WAV_BASE=2.;// Derivative order of the wavelet of choice

//int ANALYSIS=0;     
//int NBOX=128;       // Number of bins in the histogram



/* Parameters for the 1D and 2D wavelet analysis */

/* from multifractal_generator.c */
int NBOX=128;       // Number of bins in the histogram
const float DISP_ETA=0.5;     // Log-Normal processing variables
const double D0[3]= {1., .5,  1.}; // Minimum scale for each wavelet.
/**/

/* from multifractal.c */
/* minimum scales (SC0) and scale steps (QS).
 * They depend on the derivative order (first argument) and on
 * the wavelet, defined by WAV (second argument) 
 */
/* wavelet exponents 1D  (negative for gaussian) */
const double EXP_WL_1D[4] = { -1.  , 0.5  , 1.,  1.5   }; 
const double SC0_1D[3][4] =  
{ { 1.000, 1.000, 1.000, 1.000}, // Derivative order 0
  { 2.000, 2.000, 2.000, 2.000}, // Derivative order 1
  { 4.000, 4.000, 4.000, 4.000}  // Derivative order 2
}; 

/* wavelet exponents 2D  (negative for gaussian) */
const double EXP_WL_2D[4] = { -1.  , 1.,  1.5  , 2.  }; 
const double SC0_2D[3][4] =  
{ { 1.000, 0.500, 0.500, 1.0000}, // Derivative order 0
  { 2.000, 2.000, 2.000, 2.0000}, // Derivative order 1
  { 4.000, 4.000, 4.000, 4.0000}  // Derivative order 2
}; 
/**/

const double QS_1D[3][4] =   
{ { 1.750, 1.125, 1.750, 1.500}, // Derivative order 0
  { 3.000, 3.000, 3.000, 3.000}, // Derivative order 1
  { 3.000, 2.000, 2.000, 3.000}  // Derivative order 2
};
const double QS_2D[3][4] =   
{ { 1.750, 1.125, 1.250, 1.500}, // Derivative order 0
  { 3.000, 3.000, 3.000, 3.000}, // Derivative order 1
  { 3.000, 2.000, 2.000, 3.000}  // Derivative order 2
}; 

/* from Dh_evaluation.c */
int TYPE=0;        // Type of multifractal to be processed.
                   // 0: Log-Poisson
                   // 1: Log-Normal
                   // 2: Log-Levi
                   // 3: Binomial
int LEFF=512;      // Length/size of the series/images
int FROM_DH=0;     // Flag. If activated, the program retrieves the D(h)
                   // from previously recorded files, then calculates errors
int GEO_MAP=0;     // Flag. If activated, a map of qualities for series of the
                   // defined type and varying geometries is created
                   // (see parameters defining its geometry below)
int TYPE_MAP=0;    // Flag. If activated, a map of fixed geometry (1x16384)
                   // and varying parameters, according to the type selected,
                   // is generated (see the parameter choice below)
int NBOX=128;       // Initial number of bins in the histogram
const int Nh=20; // Number of points to be solved in the spectrum
//double HINF=-0.5;   // Most singular exponent (Types 0, 3 and 4)
//double CODINF=1.;   // Most singular co-dimension (Types 0 and 3)
//double H1=0.5;      // Upper singularity bound in binomial model
//double MU=0.5;      // Average singularity (for log-Normal and log-Levi)
//double SIGMA=1.;    // Singularity dispersion in log-Normal and log-Levi 
//double ALPHA=1.5;   // Exponent defining log-Levi
/*           Moment and WTMM methods    */
/* Define the method used for estimating the Legendre transform */
#define Nmom 65
double moms[Nmom]={-4.0,-3.6,-3.2,-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.3,2.6,3.0,3.5,4.0,5.0,6.0,7.0,8.0};
#define Ndist 10
int dist[Ndist]={4,5,6,8,10,12,15,18,22,30};
/*       Defining geometry for GEO_MAP    */
#define Nnums 5
#define Nleffs 4
int nums[Nnums]={1,10,100,1000,10000};
int leffs[Nleffs]={1024,4096,16384,65536};
/*          Defining parameters for TYPE_MAP    */
#define NLPs 10
const double hinfs[NLPs]=   {-.8, -.6, -.4, -.2, -.6, -.4, -.2, -.4, -.2, -.2};
const double codinfs[NLPs]={ 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.8, 0.6, 0.6, 0.4};
#define Nmeans 5
#define Nsigmas 4
const double means[Nmeans]= {-1.,-0.5,0,0.5,1.};
const double sigmas[Nsigmas]={0.33,0.5,1.,2.};
double DH;      // Accepted conventional error
/*       Internal global variables                  */
double shiftw;
double shiftg;


/**/

double *prob_levi; // To keep a copy of a Levi distribution of parameter alpha
double *prob_exp; // To keep experimental track of the distribution
double m_th[2]; // Extremes for h distribution
double acprob=0.;   // accumulated probability (for scarce prob. processes)

/** Variables associated to the file multifractal.c */

/*       Common to 1 and 2D analyses    */

//int NPOINTS=6;      
int WAV=0;          // wavelet of choice:
                    // 0: gaussian
                    // 1: Lorentzian, order 0.5
                    // 2: Lorentzian, order 1
                    // 3: Lorentzian, order 1.5

int ORDDER=0;      // Differentiation order
//int HOLDER=0;       // Analizes de measure (0) or the own function
double S0=1.;        // Initial scale for the wavelets

/*          Specific to 2D analysis                        */

int DIM1=0;         // Flag. If set, the analysis is performed in 1D
double THETAU=0.;     // Direction for the 1D analysis

/*  
minimum scales (SC0) and scale steps (QS). They depend on the derivative 
order (first argument) and on the wavelet, defined by WAV (second argument) 
*/

/* int FLAG_SAVEDOUBLE; /* Flag for saving the data series in double format
* Default: DISABLE flag, the data are saved in double format */

//int FLAG_SAVEFLOAT=0;



#endif 	    /* !VAR_MSM_H_ */
