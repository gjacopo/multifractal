/*  Program MF-analyzer.c. Version: October 11th, 2006 */

/*  Authors: Antonio Turiel and Oriol Pont         */

/*  Directorio de grabacion:  Z.OUTGOING           */

/*
    Provides: MF-analyzer. Requires: standard libraries + NetCDF library
    for full I/O support. You should install first NetCDF library on your
    computer. Alternatively, you can modify the read_data and write_data
    routines to get a non-NetCDF code.
*/

/*
    This program is supplied with no warranty, nor explicit either implicit.
    Neither the authors or their institutions are liable of responsibility
    for the consequences of the use of this code.
*/

/*
    You can freely use this code for your own calculations and modify it at
    will. The authors will be grateful to you if you recognise the source,
    authory and, eventually, credit some appropriate reference of them.
*/

/*
    Description: A huge, stand-alone piece of source C code, devoted to the
    obtention of singularity exponents and singularity spectra from given
    data. The program does not support all the possible functionalities;
    please search information on the Jfluid project for a more performing
    program.

    The program uses some simple, basic wavelets to perform singularity
    analysis; the wavelets are defined in the wavelet routine.

    The program accepts an arbitrarily long list of filenames, then it
    processes them together. Each data piece generates: an appropriate
    summary file on the appropriateness of the application of Microcanonical
    Multifractal Formalism on that piece of data; an output NetCDF file with
    the associated local singularity exponents at each position;
    an ASCII file with the associated histogram of singularities; and
    an ASCII file with the associated singularity spectrum. In addition,
    a standardized estimation of the singularity histogram and of the
    singularity spectrum associated to the whole datalist is generated.

    The program accepts some run-time parameters which modify the pace of
    execution; type MF-analyzer -h for a wider explanation.
*/

/*
   WARNING: This is a large text file, as all routines are build in this
   single file. This makes this file huge and hard to read, but easens
   compilation and compatibility for less experienced users. As we intend
   to disseminate science rather than perform academic programming, this
   is our choice. Suggestions on ways to improve future releases are of
   course welcome.
*/


/* The code starts here */


/*-------------------*/
/* HEADER DIRECTIVES */
/*-------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <netcdf.h>


/*-------------------------*/
/* SOME USEFUL DEFINITIONS */
/*-------------------------*/

#define MAX_NAME 90        // Maximum lenght for typical names
#define MAX_EXPLAIN 1024   // Maximum lenght for explanatory texts
#define MASKED ((char)255) // Value for masking exclosed points
#define VALID  ((char)0)   // Value to identify valid points
#define MAX_NDATA 10000    // Maximum number of data files possible
#define ASC_LEN 2000       // Maximum line length in ASCII files

#define XMAXVH 1536
#define YMAXVH 1024        // Dimensions of van Hateren files


/*---------------*/
/* BASIC STRUCTS */
/*---------------*/

/* Struct to contain data under analysis */
typedef struct
{
    char Name[MAX_NAME]; // Name of the variable
    int dimx;            // X size of the data matrix
    int dimy;            // Y size of the data matrix; 1 for 1-dimensional data
    double **Signal;     // Pointer to the array of values
    char **Mask;         // Mask array, in case data are masked
    double *p;           // Inner pointer to vectorized data
    char *pm;            // Inner pointer to vectorized mask

} Matrix;
#define _Matrix {"",0,0,NULL,NULL,NULL,NULL} // Matrix initalizer


/* Data reading function pointer prototype */
typedef int(*ReadMatrix)(char*, DataWindow, Matrix*);
        // char *filename, DataWindow select, Matrix* data

/* Struct to define the geometry of the processed file */
typedef struct
{
    ReadMatrix Reader; // Pointer to the reading routine to be used
  int dimx;
    int dimy;
    int dimv;
    int dimz; // Dimensions of the file

} DataStructure;
#define _DataStructure {NULL,0,0,0,0} // DataStructure initializer

/* Definition of the parse atom (for parse_argument routines) */
typedef struct
{
    char  argname[MAX_NAME];    // name of the argument
    char  valname[MAX_NAME];    // name of its value
    char  explain[MAX_EXPLAIN]; // explanation (as used by help)
    int   type;                 // type of variable: 0: flag, 1: int; 2: float, +3 if flagged, 7: string
    int   *flag;                // pointer to variable, if flag
    int   *var_i;               // pointer to variable, if int
    int   minvar_i;             // minimum value, if int
    int   maxvar_i;             // maximum value, if int
    float *var_f;               // pointer to variable, if float
    float minvar_f;             // minimum value, if float
    float maxvar_f;             // maximum value, if float
    char  *string;              // pointer to variable, if string

} Parse_atom;

/* Definition of the multifractal analysis summary list */
typedef struct
{
    char filename[MAX_NAME];
    int iv;
    int iz;
    double hmin;
    double hmax;
    double densgood; // density of good micro-canonical MF points
    double sc0;      // basic resolution scale

} MFList;
#define _MFList {"",0,0,1e30,-1e30,0.,1.} // MFList initializer


/*------------------------------------------*/
/*            RUN-TIME VARIABLES            */
/*------------------------------------------*/

// General
int VERBOSE=0;  // Flag. If enabled, information is printed to stdout
int SAVEDATA=0; // Flag. If enabled, data, singularities and h histogram are saved to disk, also individual spectra

// Data reading
int MODULUS=0;  // Flag. If set, the modulus of the vector data is taken
int LOGDATA=0;  // Flag. If set, data are logarithmized
int NDATA=1;    // Number of data files that the program expects to process
int AUTO_IN=0;  // File names are automatically read from stdin: ls *.nc | MF-analyzer.exe
                // (AUTO_EXT is assumed and NDATA is ignored)
int AUTO_EXT=0; // Extension is the default or command-line specified, not asked in runtime
char EXTENSION[MAX_NAME]="MF"; // Output extension for the generated files

int X0=0;       // Leftmost corner X-coordinate of the processed windows
int Y0=0;       // Bottonmost corner Y-coordinate of the processed windows
int V0=0;       // Lowest vectorial index of the processed windows
int Z0=0;       // Earliest time instant of the processed sequence
int XMAX=0;     // If non-zero, defines X size of the processed window
int YMAX=0;     // If non-zero, defines Y size of the processed window
int VMAX=0;     // If non-zero, defines V size of the processed window
int ZMAX=0;     // If non-zero, defines time extent of the processed sequence

int FORCEMASK=1;// Flag. Forces exponent matrices to have a mask (useful to exclude edge effects).
                // If unmasked data are processed, setting it to 0 can save a lot of memory.
                // For masked data, it has no effect.
int MEANMASK=0; // Flag. If gradients are processed, it sets all masked points to the mean gradient value,
                // to better estimate the multifractal measure of points near the mask and the edges. It may be
                // useful with unmasked data if FORCEMASK is not set. Avoid it when processing masked signals.

float BLOCK=1.; // Resolution change parameter; the number of points is
                // decreased when BLOCK > 1 and increased when BLOCK < 1

// Signal, gradient and singularity histograms
float SMIN=0.;  // Minimum signal value considered in histograms
float SMAX=0.;  // Maximum signal value considered in histograms
float GMIN=0.;  // Minimum gradient value considered in histograms
float GMAX=0.;  // Maximum gradient value considered in histograms
float HMIN=-1.; // Minimum singularity value considered in histograms
float HMAX= 2.; // Maximum singularity value considered in histograms
int NBOX=100;   // Number of histogram boxes

// Derivative parameters
int DER_MODE=0; // Mode use for calculating derivatives
                // 0: Half-pixel centered difference
                // 1: Shifted finite difference

// Analysis wavelet parameters
int PUNCTUAL=0;      // Flag. If activated, only one scale is taken into account
int WVINDEX=0;       // Index of the wavelet of choice
                     // 0: Gaussian
                     // 1: Lorentzian of order 0.5*d
                     // 2: Lorentzian of order 0.5 + 0.5*d
                     // 3: Lorentzian of order 1.0 + 0.5*d
int HOLDER=0;        // Flag. If set, the function itself is analyzed,
                     // instead of the measure
int WVPOINTS=7;      // Number of scales in the multifractal evaluation
float WAV_RANGE=10.; // Range of scales in wavelet analysis
int ORD_DER=0;       // Differentiation order on the wavelet
float GOODRHO=0.9;   // Goodness for singularity analysis regressions
int NON_TRANS=0;     // If set, data are not considered as translationally
                     // invariant. Three modes are available: 1 for whitened data
                     //                                       2 for data with 1/f^2 power spectrum
                     //                                      -1 the singularity spectrum is shifted to ensure
                     //                                         the scale invariance of the first order moment
                     // It has been found useful with the Benzi Model, but uneeded with real physical signals

/*------------------------------------------*/
/*     NON-EXTERNALLY ACCESSED PARAMETERS   */
/*------------------------------------------*/

const double SC0[3][4] =
{ { 1., 1., 1., 1.}, // Derivative order 0
  { 4., 4., 4., 4.}, // Derivative order 1
  { 8., 8., 8., 8.}  // Derivative order 2
};
//    Scale unit in the generate_wavelet routine

const double SC_W[3][4] =
{ {  16., 1.,  2.,  4.}, // Derivative order 0
  {  64., 4.,  8., 16.}, // Derivative order 1
  { 128., 8., 16., 32.}  // Derivative order 2
}; // Effective scales associated to the wavelets

const char originname[] =  "Program MF-analyzer v1.0";
const char factoryname[] = "Derived product by Institut de Ciencies del Mar, Barcelona (CSIC) & University of Barcelona";


/*------------------------------------------*/
/*          FUNCTION PROTOTYPES             */
/*------------------------------------------*/

/*          General purpose functions       */

int main(int argc, char *argv[]);
void parse_arguments(int argc, char *argv[]);

/*          Mathematical and other          */

double fMax( double a, double b);
double fMin( double a, double b);
int Max( int a, int b);
int Min( int a, int b);
int Mod( int a, int b);
double linear_fit( int N, double *x, double *y, double *a, double *b, double *corr);
int get_base( char *name, char splitter, char *base);

/*          Memory management routines      */

int fit_matrix( Matrix *data);
int initiate_matrix( Matrix a, Matrix *b);
int assign_matrix( Matrix a, Matrix b);
int change_resolution( double block, Matrix *data);

/*          Matricial routines              */

double matrix_mean( Matrix data);
double matrix_pnorm( int mode, double p, Matrix data);

/*          I/O routines                    */

DataStructure ping_file( char *filename);
int select_data_window( DataStructure InCourse, DataWindow *select);

// NetCDF I/O routines
DataStructure ping_netcdf( char *filename);
int simple_netcdf_reader( char *filename, DataWindow select, Matrix* data);
int nc_map_dv( int ncid, int ndims, int nvars, int *map_d, int *map_v);
int write_simple_netcdf( char *dataname, Matrix data);

// Unformatted van Hateren (read only) routines
DataStructure ping_vH( char *filename);
int read_vH( char *filename, DataWindow select, Matrix *data);

// ASCII I/O routines
DataStructure ping_ascii( char *filename);
int read_ascii( char *filename, DataWindow select, Matrix *data);
int column_count_ascii( char *filename, int *skip, int *Ncol, int *Nrow);
int line_column_count_ascii( FILE *chan, int *goahead);
int read_line_ascii(FILE *chan, int *Ncol, char *line);

/*          FFT routines                    */

int FFT( Matrix functionR, Matrix functionI, int sign);
int FFThorizontal( Matrix functionhR, Matrix functionhI, int sign);
void PFFThorizontal( int dima, int dimb, int dimc, int iy, double **uinR, double **uinI, double **uoutR, double **uoutI, int sign);
int FFTvertical( Matrix functionvR, Matrix functionvI, int sign);
void PFFTvertical( int dima, int dimb, int dimc, int ix, double **uinR, double **uinI, double **uoutR, double **uoutI, int sign);
int assign_convolution( Matrix f1, Matrix f2);

/*          Derivative routines             */

int gradient( Matrix gx, Matrix gy);
int gradient_naif( Matrix gx, Matrix gy);
int gradient_FFT( Matrix gx, Matrix gy);
int modgradient( Matrix signal);

int reconstruct( Matrix gx, Matrix gy);
int reconstruct_naif( Matrix gx, Matrix gy);
int reconstruct_FFT( Matrix gx, Matrix gy);

/*          Multifractal routines           */

double scale_linear( int dimx, int dimy);
double scale_wavelet( int dimx, int dimy);
int generate_wavelet( int dimx, int dimy, double sc, Matrix *wave,
		      int normalize_flag);
MFList calculate_multifractal( Matrix data,  Matrix *expon);
double non_translational_correction( Matrix data);
void Gaussian_derivative( double sc, Matrix data, Matrix *modg);
void wavelet_of_choice( char *output);
void standard_summary( char *filename, MFList summary);

void singularity_histogram( Matrix expon, double *h, double *histoh, double sc0);
void histogram( Matrix data, double *val, double *histo, double min, double max);
void histo_record( char *filename, double *val, double *histo, double min, double max);
void histoh_adapt( double sc0, double sc1, double *h, double *histoh);
void Dh_record( char *filename, double sc0, double *h, double *histoh);


/*------------------------------------------*/
/*          FUNCTION SPECIFICATIONS         */
/*------------------------------------------*/

/* General purpose functions */

int main(int argc, char *argv[])
{

/* INPUT DATA */

 Matrix data=_Matrix;
 char **name_in;

/* QUANTITIES TO BE CALCULATED */

 Matrix expon=_Matrix;
 MFList summary=_MFList;

 double *sig,*histosig;
 double *sig0,*histosig0;
 double sigmin,sigmax;

 double *grad,*histograd;
 double *grad0,*histograd0;
 double gradmin,gradmax;

 double *h,*histoh;
 double *h0,*histoh0;

/* AUXILIARY VARIABLES */

 DataStructure InCourse;
 DataWindow select;
 char anyname[MAX_NAME],exteff[MAX_NAME];
 Matrix aux=_Matrix;
 Matrix data_grad=_Matrix;
 double sc0=-1.;
 double expshift;
 int dimv,dimz;
 int iv0,iz0;
 int imod,imod0,dimmod;
 int in,iz,iv,ih;
 int ix,iy,ip;
 int n_data, ndatafiles, window_selected;
 int load_error, nomask;

/* Program */

/* On-line parameter reading */

 parse_arguments(argc,argv);

/* Filename reading */

 if (AUTO_IN)
 {
  /*
  In AUTO_IN mode, the number of data files is auto adjusted (NDATA is ignored)
  when stdin has EOF control, i.e., "ls *.nc | MF-analyzer.exe -a" or "MF-analyzer.exe -a < file_list".
  With keyboard input stdin, the control sequence "***" is used to specify that the last file has been entered.
  This mode is specially helpful when processing lots of files.
  The default extension is used, except for keyboard input mode (where it is asked)
  or if it has been specified in the command line.
  */
  printf("Reading the data filenames. Type *** to end the filename reading.\n");
  ndatafiles = MAX_NDATA;
 }
 else ndatafiles = NDATA; // the number of data files is the default (1) or specified in the command line

 name_in = (char **)calloc(ndatafiles, sizeof(char *));
 for (in=0; in<ndatafiles; in++)
 {
  name_in[in] = (char *)calloc(MAX_NAME, sizeof(char));
  printf("Type the name of the #%d data file, please\n", in+1);
  if (fscanf(stdin, "%s", name_in[in]) == EOF)
  {
   if (AUTO_EXT==0)
   {  // In case that it has not been specified, the default extension will be used
    AUTO_EXT=1;
    printf("Using the default extension: %s\n", EXTENSION);
   }
   break;
  }
  if ( strcmp(name_in[in], "***")==0 ) // The sentinel sequence specifying the last file
  {
   break;
  }
 }
 NDATA = in; // The last input is ignored, as it does not contain a file name
 if (AUTO_EXT==0) // If it has not been specified, ask the extension
 {
  printf("Type the output extension, please\n");
  scanf("%s", EXTENSION);
 }
 printf("Thank you.\n\n");

 if (VERBOSE)
 {
  printf("Multifractal analysis\n");
  printf("=====================\n");
  if (PUNCTUAL==0)
  {
   wavelet_of_choice(anyname);
   printf("Wavelet of choice: %s\n",anyname);
  }
  else
  {
   printf("Punctual method\n");
  }
 }

 if (VERBOSE && NDATA>1)
 {
  printf("Analyzing an ensemble of %d data file(s):\n", NDATA);
  for (in=0;in<NDATA;in++) printf(" #%d: %s\n", in+1, name_in[in]);
 }

 if (VERBOSE) printf("\n");

/* Memory initalization */

 sig       = (double *)calloc(NBOX,sizeof(double));
 histosig  = (double *)calloc(NBOX,sizeof(double));
 sig0      = (double *)calloc(NBOX,sizeof(double));
 histosig0 = (double *)calloc(NBOX,sizeof(double));
 sigmin=1e30; sigmax=-1e30;

 grad       = (double *)calloc(NBOX,sizeof(double));
 histograd  = (double *)calloc(NBOX,sizeof(double));
 grad0      = (double *)calloc(NBOX,sizeof(double));
 histograd0 = (double *)calloc(NBOX,sizeof(double));
 gradmin=1e30; gradmax=-1e30;

 h       = (double *)calloc(NBOX,sizeof(double));
 histoh  = (double *)calloc(NBOX,sizeof(double));
 h0      = (double *)calloc(NBOX,sizeof(double));
 histoh0 = (double *)calloc(NBOX,sizeof(double));

 expshift=0.; n_data=0;

  for (in=0;in<NDATA;in++)
  // Loops on the list: files are read and processed, then the histograms are accumlated
  {
   if (VERBOSE) printf("Processing file #%d:\n",in+1);

/* Testing the readibility of the file in course */

   InCourse=ping_file(name_in[in]);
   if (InCourse.Reader==NULL)
   {
    if (VERBOSE) printf("Unable to open file %s for reading, skipping it\n",name_in[in]);
    continue;
   }

/* Informing the user about the observed file geometry */

   if (VERBOSE)
   {
    printf("File %s open:\n",name_in[in]);
    if (InCourse.dimy==1)
     printf(" One-dimensional data record with %d elements\n",InCourse.dimx);
    else
     printf(" Two-dimensional data record with %dx%d elements\n",InCourse.dimx,InCourse.dimy);
    if(InCourse.dimv>1)
     printf(" Each element has %d dimensions\n",InCourse.dimv);
    if(InCourse.dimz>1)
     printf(" The file contains a sequence with %d images\n",InCourse.dimz);
   }

/* Re-adjusting data geometry to the desired processing window */

   window_selected = select_data_window(InCourse,&select);

   if (VERBOSE && window_selected)
   {
    printf("The specified hyper-cube has been taken into account\n");
    if (InCourse.dimx>1)
     printf("X range will go from %d to %d\n",select.ix0,select.ix0+select.dimx-1);
    if (InCourse.dimy>1)
     printf("Y range will go from %d to %d\n",select.iy0,select.iy0+select.dimy-1);
    if (InCourse.dimv>1)
     printf("V range will go from %d to %d\n",select.iv0,select.iv0+select.dimv-1);
    if (InCourse.dimz>1)
     printf("Z range will go from %d to %d\n",select.iz0,select.iz0+select.dimz-1);
   }

/*
   Although a more flexible programming is possible, we have decided
   here to independentize each component and time instant in the
   calculations; for that reason we modify the DataWindow select to
   just take a single component and time instant at once.
*/

   if (MODULUS)
   {
    iv0=0;
    dimv=1;
    imod0=select.iv0;
    dimmod=select.dimv;
   }
   else
   {
    iv0=select.iv0;
    dimv=select.dimv;
   }
   iz0=select.iz0;
   dimz=select.dimz;

   select.dimv=1;
   select.dimz=1;

/* the multicomponent/multisequence loop of readings is started here */

   for (iv=iv0;iv<iv0+dimv;iv++)
   {
   for (iz=iz0;iz<iz0+dimz;iz++)
   {
    select.iv0=iv;
    select.iz0=iz;
    get_base(name_in[in],'.',exteff);
    if (dimv>1) sprintf(exteff,"%s-V%02d",exteff,iv);
    if (dimz>1) sprintf(exteff,"%s-Z%02d",exteff,iz);

/*  Now, data are loaded  */

    load_error = 0;

    if (MODULUS)
    {
     create_matrix(select.dimx,select.dimy,1,&data); // Mask is assumed (as there might be masked components)
     for (imod=imod0;imod<dimmod;imod++)
     {
      if (VERBOSE) printf("Integrating %d component in modulus\n",imod);
      select.iv0=imod;
      if ( InCourse.Reader(name_in[in],select,&aux) )
      {
       load_error=-1;
       break;
      }
      for (ip=0;ip<aux.dimx*aux.dimy;ip++) data.p[ip] += aux.p[ip]*aux.p[ip];
      if (aux.pm!=NULL) for (ip=0;ip<aux.dimx*aux.dimy;ip++) if (aux.pm[ip]==MASKED) data.pm[ip]=MASKED;
      destroy_matrix(&aux);
     }
     for (ip=0;ip<data.dimx*data.dimy;ip++) data.p[ip]=sqrt(data.p[ip]);
    }
    else
    {
     load_error = InCourse.Reader(name_in[in],select,&data);
    }

    if (load_error)
    {
     printf("Warning! Read error found while reading %s\n", name_in[in]);
     continue; // Try with the next component/sequence (maybe not all the file is corrupted)
    }

    if ( VERBOSE && (dimv>1||dimz>1) ) printf("New data slice read\n");

/*  Resolution is changed, if required  */

    if (fabs(BLOCK-1.)>MIN_NUM_RESOLUTION)
    {
     change_resolution(BLOCK,&data);
     if (VERBOSE)
     {
      printf("Resolution of data changed by a factor %0.2f\n",BLOCK);
      printf("Data now contain %d x %d points\n",data.dimx,data.dimy);
     }
    }

/*  If data must be logarithmized, this is the moment for doing so  */

    if (LOGDATA)
    {
     if (VERBOSE) printf("The logarithm of data is taken\n");
     for (iy=0;iy<data.dimy;iy++)
     {
     for (ix=0;ix<data.dimx;ix++)
     {
      if (fabs(data.Signal[iy][ix]) > MIN_NUM_RESOLUTION)
       data.Signal[iy][ix] = log( fabs(data.Signal[iy][ix]) );
      else
       data.Signal[iy][ix] = log(MIN_NUM_RESOLUTION);
     }
     }
    }

/*  Normalized and gradient data are recorded. Their extrema are obtained.  */

    if (SAVEDATA)
    {
     if (SMIN>=SMAX) // so no boundaries (or incoherent) have been specified
     {
      for (ip=0;ip<data.dimx*data.dimy;ip++)
      {
       nomask=1; if (data.pm!=NULL) if (data.pm[ip]==MASKED) nomask=0;
       if (nomask)
       {
        sigmin = fMin(sigmin,data.p[ip]);
        sigmax = fMax(sigmax,data.p[ip]);
       }
      }
     }
     else // the specified values are used
     {
        sigmin = SMIN;
        sigmax = SMAX;
     }

     sprintf(anyname, "normalized_%s", exteff);
     write_simple_netcdf(anyname, data);

     create_matrix(data.dimx, data.dimy, (data.pm!=NULL || FORCEMASK)?1:0, &data_grad);
     assign_matrix(data, data_grad);
     modgradient(data_grad);

     if (GMIN>=GMAX)
     {
      for (ip=0;ip<data_grad.dimx*data_grad.dimy;ip++)
      {
       nomask=1; if (data_grad.pm!=NULL) if (data_grad.pm[ip]==MASKED) nomask=0;
       if (nomask)
       {
        gradmin = fMin(gradmin,data_grad.p[ip]);
        gradmax = fMax(gradmax,data_grad.p[ip]);
       }
      }
     }
     else
     {
        gradmin = GMIN;
        gradmax = GMAX;
     }

     sprintf(anyname, "gradient_%s", exteff);
     write_simple_netcdf(anyname, data_grad);

     destroy_matrix(&data_grad);
    }

/*  MULTIFRACTAL PROCESSING STARTS HERE  */

/*  First step: multifractal exponents are calculated and recorded  */

    summary=calculate_multifractal(data,&expon);
    if (SAVEDATA)
    {
     sprintf(anyname,"expon_%s",exteff);
     sprintf(expon.Name,"exponents");
     write_simple_netcdf(anyname,expon);
    }
    if (sc0<0.) sc0=summary.sc0; // The first time, resolution is kept as basic resolution

/*  Recording the associated MF output summary  */

    sprintf(summary.filename,"%s",name_in[in]);
    if (dimv>1) summary.iv=iv+1;
    if (dimz>1) summary.iz=iz+1;
    sprintf(anyname,"MF-summary_%s.txt",exteff);
    standard_summary(anyname,summary);

/*  Generating the partial singularity histogram and singularity spectrum  */

    singularity_histogram(expon,h,histoh,summary.sc0);

    if (SAVEDATA)
    {
     sprintf(anyname,"histoh_sng_%s.txt",exteff);
     histo_record(anyname,h,histoh,HMIN,HMAX);
    }
    sprintf(anyname,"Dh_sng_%s.txt",exteff);
    Dh_record(anyname,summary.sc0,h,histoh);

/*  Accumulating the partial histogram in the global one  */

/*  Histograms must be adapted to conform the same resolution scale  */

    histoh_adapt(sc0,summary.sc0,h,histoh);

    for (ih=0;ih<NBOX;ih++)
    {
     h0[ih]      += h[ih];
     histoh0[ih] += histoh[ih];
    }

/* Evaluating shift associated to non tranlational invariance */

    if (NON_TRANS>0)
    {
     expshift += non_translational_correction(data);
     n_data++;
    }

/*  Memory release before the end of the loop  */

    destroy_matrix(&expon);
    destroy_matrix(&data);
   }
   } // End of multicomponent/multisequence processing (iz, iv)
   if (VERBOSE) printf("\n");
  } // End of multisignal processing (in)

/*
   The process is now repeated to obtain the signal and gradient histograms.
   The code is the same as in the loop above, but more compact.
*/

 if (SAVEDATA)
 {
  if (VERBOSE) printf("Computing signal and gradient histograms\n");
  for (in=0;in<NDATA;in++)
  {
   InCourse=ping_file(name_in[in]); if (InCourse.Reader==NULL) continue;
   select_data_window(InCourse,&select);
   if (MODULUS)
   {
    iv0=0; dimv=1;
    imod0=select.iv0; dimmod=select.dimv;
   }
   else
   {
    iv0=select.iv0; dimv=select.dimv;
   }
   iz0=select.iz0; dimz=select.dimz;
   select.dimv=1; select.dimz=1;

   for (iv=iv0;iv<iv0+dimv;iv++)
   {
   for (iz=iz0;iz<iz0+dimz;iz++)
   {
    select.iv0=iv; select.iz0=iz;
    load_error = 0;
    if (MODULUS)
    {
     create_matrix(select.dimx,select.dimy,1,&data);
     for (imod=imod0;imod<dimmod;imod++)
     {
      select.iv0=imod;
      if ( InCourse.Reader(name_in[in],select,&aux) )
      {
       load_error=-1; break;
      }
      for (ip=0;ip<aux.dimx*aux.dimy;ip++) data.p[ip] += aux.p[ip]*aux.p[ip];
      if (aux.pm!=NULL) for (ip=0;ip<aux.dimx*aux.dimy;ip++) if (aux.pm[ip]==MASKED) data.pm[ip]=MASKED;
      destroy_matrix(&aux);
     }
     for (ip=0;ip<data.dimx*data.dimy;ip++) data.p[ip]=sqrt(data.p[ip]);
    }
    else load_error = InCourse.Reader(name_in[in],select,&data);
    if (load_error) continue;
    if (fabs(BLOCK-1.)>MIN_NUM_RESOLUTION) change_resolution(BLOCK,&data);
    if (LOGDATA) for (ip=0;ip<data.dimx*data.dimy;ip++)
    {
     if (fabs(data.p[ip])>MIN_NUM_RESOLUTION) data.p[ip]=log(fabs(data.p[ip]));
     else                                       data.p[ip]=log(MIN_NUM_RESOLUTION);
    }

/*  Signal and gradient histograms are made and accumulated  */

    create_matrix(data.dimx, data.dimy, (data.pm!=NULL || FORCEMASK)?1:0, &data_grad);
    assign_matrix(data, data_grad);
    modgradient(data_grad);

    histogram(data,sig,histosig,sigmin,sigmax); // min and max seen in the previous stage
    histogram(data_grad,grad,histograd,gradmin,gradmax);

    for (ih=0;ih<NBOX;ih++)
    {
     sig0[ih]      += sig[ih];
     histosig0[ih] += histosig[ih];

     grad0[ih]      += grad[ih];
     histograd0[ih] += histograd[ih];
    }

    destroy_matrix(&data_grad);
    destroy_matrix(&data);
   }
   }
  } // End of 2nd read
 } // End if (SAVEDATA)

 if (NON_TRANS > 0) // The global non-translational singularity shift is corrected
 {
  expshift /= (double)n_data;

  for (ih=0; ih<NBOX; ih++) h0[ih] += expshift*histoh0[ih];

  if (VERBOSE) printf("The global exponent shift due to non-translational effects is %0.2f\n", expshift);
 }

/* Recording the results associated to the join statistics */

 if (SAVEDATA)
 {
  sprintf(anyname,"histosig_%s.txt",EXTENSION);
  histo_record(anyname,sig0,histosig0,sigmin,sigmax);

  sprintf(anyname,"histograd_%s.txt",EXTENSION);
  histo_record(anyname,grad0,histograd0,gradmin,gradmax);

  sprintf(anyname,"histoh_%s.txt",EXTENSION);
  histo_record(anyname,h0,histoh0,HMIN,HMAX);
 }
 sprintf(anyname,"Dh_%s.txt",EXTENSION);
 Dh_record(anyname,sc0,h0,histoh0);

/* Memory release and end */

 free(sig);  free(histosig);
 free(sig0); free(histosig0);

 free(grad);  free(histograd);
 free(grad0); free(histograd0);

 free(h);  free(histoh);
 free(h0); free(histoh0);

 for (in=0;in<NDATA;in++) free(name_in[in]);
 free(name_in);

 exit(0);
}

void parse_arguments(int argc, char *argv[])
/*
   PURPOSE: To parse on-line arguments and modify the associated global
   variables, according to a pre-established scheme, define by an array
   of Parse_atoms (here named line).
*/
{
   Parse_atom *line;
   const int MaxNarg=100;
   int Narg,in;
   int lar,flagv;

/* Defining the accepted run-time arguments */

   line=(Parse_atom *) calloc(MaxNarg,sizeof(Parse_atom));
   in=0;


//    Argument SAVEDATA. Type 0: flag

    sprintf(line[in].argname,"%s","-full");
    sprintf(line[in].explain," %s : %s\n",line[in].argname,
     "Flag. If enabled, the program shows outputs many intermediate files.\nDefault: DISABLED");
    line[in].type=0;
    line[in].flag=&SAVEDATA;
    in++;

//   Argument NDATA. Type 1: Integer

    strcpy(line[in].argname,"-N");
    strcpy(line[in].valname,"#files");
    sprintf(line[in].explain," %s : %s %d\n",line[in].argname,
     "Number of input data files. Default: ",NDATA);
    line[in].type=1;
    line[in].var_i=&NDATA;
    line[in].minvar_i=1;
    line[in].maxvar_i=MAX_NDATA;
    in++;

//   Argument AUTO_IN. Type 0: Flag

    strcpy(line[in].argname,"-a");
    sprintf(line[in].explain," %s : %s\n",line[in].argname,
     "Flag. Enables automated filename reading from a stdin file. Default: DISABLED");
    line[in].type=0;
    line[in].flag=&AUTO_IN;
    in++;

//   Argument EXTENSION. Type 7: String

    strcpy(line[in].argname,"-ext");
    strcpy(line[in].valname,"extension");
    sprintf(line[in].explain," %s : %s %s\n",line[in].argname,
     "Specify an extension for the generated files. Default: ",EXTENSION);
    line[in].type=7;
    line[in].flag=&AUTO_EXT;
    line[in].string=EXTENSION;
    in++;

//   Argument MODULUS. Type 0: Flag

    strcpy(line[in].argname,"-mod");
    sprintf(line[in].explain," %s : %s \n",line[in].argname,
     "Flag. If enabled, the given V components are combined to form the modulus. Default: DISABLED");
    line[in].type=0;
    line[in].flag=&MODULUS;
    in++;

//   Argument LOGDATA. Type 0: Flag

    strcpy(line[in].argname,"-log");
    sprintf(line[in].explain," %s : %s \n",line[in].argname,
     "Flag. Enables to process log data. Default: DISABLED");
    line[in].type=0;
    line[in].flag=&LOGDATA;
    in++;

//   Argument X0. Type 1: Integer

    strcpy(line[in].argname,"-ix");
    strcpy(line[in].valname,"coord");
    sprintf(line[in].explain," %s : %s %d\n",line[in].argname,
     "Leftmost X coordinate of the processing window. Default:",X0);
    line[in].type=1;
    line[in].var_i=&X0;
    line[in].minvar_i=0;
    line[in].maxvar_i=10000;
    in++;

//   Argument Y0. Type 1: Integer

    strcpy(line[in].argname,"-iy");
    strcpy(line[in].valname,"coord");
    sprintf(line[in].explain," %s : %s %d\n",line[in].argname,
     "Bottommost Y coordinate of the processing window. Default:",Y0);
    line[in].type=1;
    line[in].var_i=&Y0;
    line[in].minvar_i=0;
    line[in].maxvar_i=10000;
    in++;

//   Argument V0. Type 1: Integer

    strcpy(line[in].argname,"-iv");
    strcpy(line[in].valname,"component");
    sprintf(line[in].explain," %s : %s %d\n",line[in].argname,
     "Selected vector component of the processing window.\nFor multi-component data only. Default:",V0);
    line[in].type=1;
    line[in].var_i=&V0;
    line[in].minvar_i=0;
    line[in].maxvar_i=10000;
    in++;

//   Argument Z0. Type 1: Integer

    strcpy(line[in].argname,"-iz");
    strcpy(line[in].valname,"instant");
    sprintf(line[in].explain," %s : %s %d\n",line[in].argname,
     "Selected time instant in a sequence. For sequence data only.\n Default:",Z0);
    line[in].type=1;
    line[in].var_i=&Z0;
    line[in].minvar_i=0;
    line[in].maxvar_i=10000;
    in++;

//   Argument XMAX. Type 1: Integer

    strcpy(line[in].argname,"-sizex");
    strcpy(line[in].valname,"extent");
    sprintf(line[in].explain," %s : %s %d\n",line[in].argname,
     "Size in X direction of the processing window. Default:",XMAX);
    line[in].type=1;
    line[in].var_i=&XMAX;
    line[in].minvar_i=1;
    line[in].maxvar_i=10000;
    in++;

//   Argument YMAX. Type 1: Integer

    strcpy(line[in].argname,"-sizey");
    strcpy(line[in].valname,"extent");
    sprintf(line[in].explain," %s : %s %d\n",line[in].argname,
     "Size in Y direction of the processing window. Default:",YMAX);
    line[in].type=1;
    line[in].var_i=&YMAX;
    line[in].minvar_i=1;
    line[in].maxvar_i=10000;
    in++;

//   Argument VMAX. Type 1: Integer

    strcpy(line[in].argname,"-sizev");
    strcpy(line[in].valname,"#components");
    sprintf(line[in].explain," %s : %s %d\n",line[in].argname,
     "Number of processed vector components. For multi-component data only.\n Default:",VMAX);
    line[in].type=1;
    line[in].var_i=&VMAX;
    line[in].minvar_i=1;
    line[in].maxvar_i=10000;
    in++;

//   Argument ZMAX. Type 1: Integer

    strcpy(line[in].argname,"-sizez");
    strcpy(line[in].valname,"interval");
    sprintf(line[in].explain," %s : %s %d\n",line[in].argname,
     "Number of time instants processed. For sequence data only.\n Default:",ZMAX);
    line[in].type=1;
    line[in].var_i=&ZMAX;
    line[in].minvar_i=1;
    line[in].maxvar_i=10000;
    in++;

//   Argument BLOCK. Type 2: Float

    strcpy(line[in].argname,"-Bin");
    strcpy(line[in].valname,"factor");
    sprintf(line[in].explain," %s : %s  %0.2f\n",line[in].argname,
     "Change in resolution of processed data. The final number of data points\nis reduced (if given factor > 1) or increased (if factor < 1).\n Default:",BLOCK);
    line[in].type=2;
    line[in].var_f=&BLOCK;
    line[in].minvar_f=0.0001;
    line[in].maxvar_f=10000.;
    in++;

//   Argument NBOX. Type 1: Integer

    strcpy(line[in].argname,"-Nbox");
    strcpy(line[in].valname,"#boxes");
    sprintf(line[in].explain,"%s\n %s : %s %d\n",
     "\nHISTOGRAM VARIABLES\n===================",line[in].argname,
     "Histogram variable. Number of histogram boxes. Default:",NBOX);
    line[in].type=1;
    line[in].var_i=&NBOX;
    line[in].minvar_i=2;
    line[in].maxvar_i=100000000;
    in++;

//   Argument SMIN. Type 2: Float

    strcpy(line[in].argname,"-smin");
    strcpy(line[in].valname,"smin");
    sprintf(line[in].explain," %s : %s  %0.2f\n",line[in].argname,
     "Minimum signal value in the signal histogram. Default:",SMIN);
    line[in].type=2;
    line[in].var_f=&SMIN;
    line[in].minvar_f=-1e30;
    line[in].maxvar_f=1e30;
    in++;

//   Argument SMAX. Type 2: Float

    strcpy(line[in].argname,"-smax");
    strcpy(line[in].valname,"smax");
    sprintf(line[in].explain," %s : %s  %0.2f\n",line[in].argname,
     "Maximum signal value in the signal histogram. If min and max coincide, both are autoadjusted. Default:",SMAX);
    line[in].type=2;
    line[in].var_f=&SMAX;
    line[in].minvar_f=-1e30;
    line[in].maxvar_f=1e30;
    in++;

//   Argument GMIN. Type 2: Float

    strcpy(line[in].argname,"-gmin");
    strcpy(line[in].valname,"gmin");
    sprintf(line[in].explain," %s : %s  %0.2f\n",line[in].argname,
     "Minimum gradient value in the gradient histogram. Default:",GMIN);
    line[in].type=2;
    line[in].var_f=&GMIN;
    line[in].minvar_f=-1e30;
    line[in].maxvar_f=1e30;
    in++;

//   Argument GMAX. Type 2: Float

    strcpy(line[in].argname,"-gmax");
    strcpy(line[in].valname,"gmax");
    sprintf(line[in].explain," %s : %s  %0.2f\n",line[in].argname,
     "Maximum gradient value in the gradient histogram. If min and max coincide, both are autoadjusted. Default:",GMAX);
    line[in].type=2;
    line[in].var_f=&GMAX;
    line[in].minvar_f=-1e30;
    line[in].maxvar_f=1e30;
    in++;

    //   Argument HMIN. Type 2: Float
    //   "-hmin" : HMIN

    //   Argument HMAX. Type 2: Float
    //   "-hmax" : HMAX

    //   Argument DER_MODE. Type 1: Integer
    //   "-der_mode" : DER_MODE
    
    //   Argument WVINDEX. Type 1: Integer
    //   "-wav" : WVINDEX 

    //   Argument ORD_DER. Type 1: Integer
    //   "-der" : ORD_DER
    
    //   Argument HOLDER. Type 0: Flag
    //   "-hold" : HOLDER

//   Argument WVPOINTS. Type 1: Int
    strcpy(line[in].argname,"-wvpoints");
    strcpy(line[in].valname,"#scale_points");
    sprintf(line[in].explain," %s : %s %d\n",line[in].argname,
     "Number of scale points to be used in wavelet analysis. Default:",WVPOINTS);
    line[in].type=1;
    line[in].var_i=&WVPOINTS;
    line[in].minvar_i=3;
    line[in].maxvar_i=1000;
    in++;

//   Argument WAV_RANGE. Type 2: Float
//   "-range" : WAV_RANGE

//   Argument GOODRHO. Type 2: Float
    strcpy(line[in].argname,"-goodness");
    strcpy(line[in].valname,"regres_thresh");
    sprintf(line[in].explain," %s : %s %f\n",line[in].argname,
     "Minimum absolute value of the regression coefficient, in order to\nconsider the local scaling at point as correct in wavelet analysis.\n Default:",GOODRHO);
    line[in].type=2;
    line[in].var_f=&GOODRHO;
    line[in].minvar_f=0.;
    line[in].maxvar_f=1.;
    in++;

//   Argument NON_TRANS. Type 1: Int
    strcpy(line[in].argname,"-NonTrans");
    strcpy(line[in].valname,"mode");
    sprintf(line[in].explain," %s : %s %d\n",line[in].argname,
     "Non translational corrections:\n   0: Disabled\n   1: For whitened data\n   2: For f^{-2} data\n   -1: Shift the singularities to impose the conservation of the first order moment\nDefault:",NON_TRANS);
    line[in].type=1;
    line[in].var_i=&NON_TRANS;
    line[in].minvar_i=-1;
    line[in].maxvar_i=2;
    in++;

//   Argument PUNCTUAL. Type 0: Flag
    strcpy(line[in].argname,"-punct");
    sprintf(line[in].explain," %s : %s \n",line[in].argname,
     "Flag. If enabled, singularities are estimated punctually.\n Default: DISABLED");
    line[in].type=0;
    line[in].flag=&PUNCTUAL;
    in++;

/*      Closing argument definition    */

    Narg=in;

/*     Parsing the input to check for them        */

/*                   Parsing loop                         */

    for (lar=1,flagv=1; (lar<argc)&&(flagv==1); lar++) {
      flagv=0;
      for (in=0; (in<Narg)&&(flagv==0); in++)
	if ((!strcmp(argv[lar],"-h"))||(!strcmp(argv[lar],"-help"))) flagv=2;
	else if(!strcmp(argv[lar],line[in].argname)) flagv=1;
      
      if (flagv==1)  {
	in--;

	if (line[in].type==7)  { // String type
	  lar++;
	  sscanf(argv[lar], "%s", line[in].string);
	  *(line[in].flag)=1; // string types are assumed to be flagged

	} else if (line[in].type%3==2) {   // Float type
	  lar++;
	  sscanf(argv[lar], "%f", line[in].var_f);
	  if ((*(line[in].var_f)<line[in].minvar_f) ||
	      (*(line[in].var_f)>line[in].maxvar_f)) {
	    printf("\nValue out of range (%f - %f) for argument %s\n\n",
		   line[in].minvar_f,line[in].maxvar_f,line[in].argname);
	    flagv=0;
	  }
	  if (line[in].type==5) *(line[in].flag)=1;

	} else if (line[in].type%3==1)   { // Integer type
	  lar++;
	  sscanf(argv[lar], "%d", line[in].var_i);
	  if ((*(line[in].var_i)<line[in].minvar_i) ||
	      (*(line[in].var_i)>line[in].maxvar_i))	    {
	    printf("\nValue out of range (%d - %d) for argument %s\n\n",
		   line[in].minvar_i,line[in].maxvar_i,line[in].argname);
	    flagv=0;
	  }
	  if (line[in].type==4) *(line[in].flag)=1;

	}   else  {   // Flag type
	  *(line[in].flag)=1;
	  if (line[in].type==3) 
	    *(line[in].var_i)=1; // Flagged flag type (3) uses var_i to store the parent flag
	}
      }
    }
    
    if (flagv!=1) {
      printf("Usage: %s",argv[0]);
      for (in=0; in<Narg; in++)
	if(line[in].type%3==0) printf(" [%s]",line[in].argname);
	else printf(" [%s %s]",line[in].argname,line[in].valname);
      printf("\n");
    }
    
    if (flagv==2) for(in=0; in<Narg; in++) printf("%s",line[in].explain);
    
    /*                Termination              */
    
    if (flagv!=1) exit(-1); // parsing throws error, we stop the program
}


/*       double fMax( double a, double b)
* unchanged 
*/

/*       double fMin( double a, double b)
* unchanged 
*/

/*       int Max( int a, int b)
* unchanged 
*/

/*       int Min( int a, int b)
* unchanged 
*/

/*       int Mod( int a, int b)
* unchanged 
*/

/*       double linear_fit( int N, double *x, double *y, double *a, double *b,
 *                          double *corr );
 * replaced by :
 *       double fit( double *x, double *y, int N, double *a, double *b,
 *                   double *corr );
*/

/*       int get_base( char *name, char splitter, char *base );
 * unchanged 
 */

/*      Memory management routines        */

/***************************************************************************/
int fit_matrix( Matrix *data) {
  /***************************************************************************/
  /* PURPOSE: To give inner coherency to the *data matrix struct  */
  int iy;
  
  /* Early consistency checks  */
  if ( (data->p==NULL) || (data->dimx<1) || (data->dimy<1) ) return(-1);
  
  /* Prepares data->Signal to fit the structure */
  if (data->Signal!=NULL) free(data->Signal);
   data->Signal=(double **) calloc(data->dimy,sizeof(double *));
   
   /* Does the same with the mask if required */
   if (data->pm!=NULL)   {
     if (data->Mask!=NULL) free(data->Mask);
     data->Mask=(char **) calloc(data->dimy,sizeof(char *));
   }
   
   /* The array pointers are now directed to the starts of the corresponding
    * horizontal rows */
   for (iy=0;iy<data->dimy;iy++) {
     data->Signal[iy]=&(data->p[iy*data->dimx]);
     if (data->pm!=NULL) data->Mask[iy]=&(data->pm[iy*data->dimx]);
   }
 
   return(0);  /* Exiting */
}

/***************************************************************************/
int create_matrix( int dimx, int dimy, int simask, Matrix *data) {
  /* replaced by 
   *          Signal*  create_signal ( int dimx, int dimy, int *flag ) */
  /***************************************************************************/
  /* PURPOSE: To allocate space for a new matrix */
  int ix;
  
  /* Early checks */
  if ( (dimx<1) || (dimy<1) ) return(-1);
  
  /* Cleaning any previous structure on data */
  destroy_matrix(data);
  
  /* Allocating memory */
  sprintf(data->Name,"dummy");
  data->dimx=dimx;
  data->dimy=dimy;
  data->p=(double *) calloc(dimx*dimy,sizeof(double));
  if (simask)   {
    data->pm=(char *) calloc(dimx*dimy,sizeof(char));
    for (ix=0; ix<dimx*dimy; ix++) data->pm[ix] = VALID;
  }
  
  return(fit_matrix(data)); /* Exiting by refitting */
}

/*       int destroy_matrix( Matrix *data)
 * replaced by:
 *       int free_signal( Signal* sig ) {
 */


/***************************************************************************/
 int initiate_matrix( Matrix a, Matrix *b) {
   /* replaced by:
    *       int init_signal(Signal* sig, Signal def);  */
   /***************************************************************************/
   /* PURPOSE: Both to allocate memory on b and to copy the content of a */   
   int status=0;
   
   if((status=create_matrix(a.dimx,a.dimy,(a.pm!=NULL)?1:0,b))) 
     return status;
   else
     return (assign_matrix(a,*b));
 }

/***************************************************************************/
int assign_matrix( Matrix a, Matrix b) {
  /* replaced by:
   *       int copy_signal( Signal *sin, Signal* sdest, char **m ); */
  /***************************************************************************/
  /* PURPOSE: Copy the values of matrix a to matrix b on the common
   * overlapping area */
  int dimx,dimy;
  int ix,iy;
  
  dimx=Min(a.dimx,b.dimx);
  dimy=Min(a.dimy,b.dimy);
  
  /* Early check */  
  if ( (dimx<1) || (dimy<1) ) return(-1);
  
   /* Copying the values */
   for (iy=0;iy<dimy;iy++)  
     for (ix=0;ix<dimx;ix++)   {
       b.Signal[iy][ix]=a.Signal[iy][ix];
       if ((b.Mask!=NULL)&&(a.Mask!=NULL)) b.Mask[iy][ix]=a.Mask[iy][ix];
     }
   
   return(0);   /* Exiting */
}

int change_resolution( double block, Matrix *data)
/*
   PURPOSE: To change the resolution of the matrix data according to the
   parameter block.
*/
{
   Matrix aux=_Matrix;
   double bsize;
   double buff;
   char dmask;
   int beffx,beffy;
   int xsize,ysize;
   int ix,iy;
   int bx,by;

/* Early checks */

   if ( (data->dimx<1) || (data->dimy<1) || (data->p==NULL) ) return(-1);
   if (fabs(block-1.)<MIN_NUM_RESOLUTION) return(0); // For tiny res. changes do nothing!

/* Copying data and deallocating memory for the new shaped matrix */

   initiate_matrix(*data,&aux);
   destroy_matrix(data);

   if (block<=1.) {// Increased number of points
     
     beffx=(int) (1./block);
     beffy=(aux.dimy>1)?beffx:1;
     
     create_matrix(beffx*aux.dimx,beffy*aux.dimy,(aux.pm!=NULL)?1:0,data);
     for (iy=0;iy<aux.dimy;iy++)    
       for(ix=0;ix<aux.dimx;ix++)    {
	 
	 /* Expanding the data points */
	 for (by=0;by<beffy;by++)
	   for (bx=0;bx<beffx;bx++)	    
	     data->Signal[iy*beffy+by][ix*beffx+bx]=aux.Signal[iy][ix];
	 
	 if (aux.pm!=NULL)   // Processing the mask	 
	   for (by=0;by<beffy;by++)	     
	     for (bx=0;bx<beffx;bx++)	       
	       data->Mask[iy*beffy+by][ix*beffx+bx]=aux.Mask[iy][ix];
       }
     
   } else {// decreased number of points
     
     beffx=(int) block;
     beffy=(aux.dimy>1)?beffx:1;
     bsize=(double)(beffx*beffy);
     xsize=aux.dimx/beffx;
     ysize=aux.dimy/beffy;
     create_matrix(xsize,ysize,(aux.pm!=NULL)?1:0,data);
     
     for (iy=0;iy<ysize;iy++) {
       for (ix=0;ix<xsize;ix++)  {
	 buff=0.;
	 for (by=0;by<beffy;by++) 
	   for (bx=0;bx<beffx;bx++)  
	     buff+=aux.Signal[beffy*iy+by][beffx*ix+bx];
	 data->Signal[iy][ix]=buff/bsize;
	 if (aux.Mask!=NULL) {  // Masked points are dealt as absorbing
	   
	   dmask=VALID;
	   for (by=0;(by<beffy)&&(dmask==VALID);by++)  
	     for (bx=0;(bx<beffx)&&(dmask==VALID);bx++) 
	       dmask=aux.Mask[beffy*iy+by][beffx*ix+bx];
	   data->Mask[iy][ix]=dmask;
	 }
       }
     }
   }
   
   /* Memory release and end */
   
   destroy_matrix(&aux);
    return(0);
}

/* Matricial routines */

double matrix_mean( Matrix data)
/*
   PURPOSE: Computes the mean value
*/
{
   double mean,Nev;
   int nomask;
   int ix,iy;

   mean=Nev=0.;
   for (iy=0;iy<data.dimy;iy++) {
     for (ix=0;ix<data.dimx;ix++)   {
       nomask=1;
       if (data.Mask!=NULL)	 {
         if (data.Mask[iy][ix]==MASKED) nomask=0;
       }
       if (nomask)	 {
         Nev+=1.;
         mean+=data.Signal[iy][ix];
       }
     }
   }
   
   if (Nev>0.5) mean/=Nev;
   else         mean=1.;

   return(mean);
}

double matrix_pnorm( int mode, double p, Matrix data) {
  /* PURPOSE: Computes the p-norm. If mode=1, normalizes data */

  double norm,Nev;
  int nomask;
  int ix,iy;
  
  if (p<=0.) return(-1.);
  
  norm=0.;
  Nev=0.;
  for (iy=0;iy<data.dimy;iy++) 
    for (ix=0;ix<data.dimx;ix++) {
      nomask=1;
      if (data.Mask!=NULL)	
	if(data.Mask[iy][ix]==MASKED) nomask=0;
      
      if (nomask)	 {
	Nev+=1.;
	norm+=pow(fabs(data.Signal[iy][ix]),p);
      }
    }
  
   if (Nev>0.5) norm/=Nev;
   else         norm=1.;
   norm=pow(norm,1./p);
   
   if (mode) 
     for (iy=0;iy<data.dimy;iy++) 
       for (ix=0;ix<data.dimx;ix++)	 {
	 nomask=1;
	 if (data.Mask!=NULL)	 
	   if (data.Mask[iy][ix]==MASKED) nomask=0;
	 if (nomask) data.Signal[iy][ix]/=norm;
       }
   
   return(norm);
}

/*              I/O routines                */

DataStructure ping_file( char *filename)
/*
   PURPOSE:  To determine the geometry of the data, for later use in defining
   the appropriate DataWindow which will be effectively read. In addition, the
   routine that should be used to read the data is also recorded inside
   the output DataStructure.
*/
{
   DataStructure output=_DataStructure;

   output=ping_netcdf(filename); // First try NetCDF format
   if (output.Reader==NULL)      // if data do not accomodate NetCDF...
    output=ping_vH(filename);    // let's try for van Hatteren type
   if (output.Reader==NULL)      // not yet found...
    output=ping_ascii(filename); // looking for ASCII multicolum 1D

   return(output);
}

int select_data_window( DataStructure InCourse, DataWindow *select)
/*
   PURPOSE: Define the data window compatible with data geometry
*/
{
   int status=0;

   select->ix0=Min(X0,InCourse.dimx-1); // ix0 does not go beyond limit
   if (XMAX>0) // If a X size was defined
    select->dimx=Min(XMAX,InCourse.dimx-select->ix0);
   else
    select->dimx=InCourse.dimx-select->ix0;

   select->iy0=Min(Y0,InCourse.dimy-1); // ix0 does not go beyond limit
   if (YMAX>0) // If a Y size was defined
    select->dimy=Min(YMAX,InCourse.dimy-select->iy0);
   else
    select->dimy=InCourse.dimy-select->iy0;

   select->iv0=Min(V0,InCourse.dimv-1); // iv0 does not go beyond limit
   if (VMAX>0) // If a V size was defined
    select->dimv=Min(VMAX,InCourse.dimv-select->iv0);
   else
    select->dimv=InCourse.dimv-select->iv0;

   select->iz0=Min(Z0,InCourse.dimz-1); // iz0 does not go beyond limit
   if (ZMAX>0) // If a Z size was defined
    select->dimz=Min(ZMAX,InCourse.dimz-select->iz0);
   else
    select->dimz=InCourse.dimz-select->iz0;

   status = X0 + XMAX + Y0 + YMAX + V0 + VMAX + Z0 + ZMAX; // If a window has been specified, it is different from 0

   return(status);
}

//    NetCDF I/O routines
DataStructure ping_netcdf( char *filename);
int simple_netcdf_reader( char *filename, DataWindow select, Matrix* data);
int nc_map_dv( int ncid, int ndims, int nvars, int *map_d, int *map_v);
int write_simple_netcdf( char *dataname, Matrix data);




// Unformatted van Hateren (read only) routines
DataStructure ping_vH( char *filename)
/*
   PURPOSE: To determine the structure of .imc van Hateren files
*/
{
   DataStructure output=_DataStructure;

   FILE *canal;
   long int length;

   canal=fopen(filename,"rb");
   if (canal==NULL) return(output);
   fseek(canal,0,SEEK_SET);
   length=-ftell(canal);
   fseek(canal,0,SEEK_END);
   length+=ftell(canal);
   fclose(canal);
   if (length==3145728) // type van Hateren
   {
      output.Reader=&read_vH;
      output.dimx=XMAXVH;
      output.dimy=YMAXVH;
      output.dimv=output.dimz=1;
   }

/* Exiting the output */

   return(output);
}

int read_vH( char *filename, DataWindow select, Matrix *data)
/*
   PURPOSE: Simple reading interface to van Hateren .imc files
*/
{
   FILE* canal;

   int littleindian=0; // van Hatteren images are big endian
   int value;
   int dat1,dat2,ix,iy;

/* Allocating memory previous to reading */

   create_matrix(select.dimx,select.dimy,0,data); // has no mask

/* Picking up the selected data window */

   canal=fopen(filename,"rb");
   for (iy=0;iy<YMAXVH;iy++)
   {
   for (ix=0;ix<XMAXVH;ix++)
   {
      dat1=(int) getc(canal);
      if (dat1<0) dat1=dat1+256;

      dat2=(int) getc(canal);
      if (dat2<0) dat2=dat2+256;

      if (littleindian) value=dat1+256*dat2;
      else              value=dat2+256*dat1;

      if ((ix>=select.ix0)&&(ix<select.ix0+select.dimx)&&(iy>=select.iy0)&&(iy<select.iy0+select.dimy))
      data->Signal[data->dimy-1-(iy-select.iy0)][ix-select.ix0]=(double)value;
//                         |
//                         |
//                         |
//   van Hateren images follow the image processing convention of placing
//   the origin in the top left corner; we shift the image to accomodate
//   our convention and so all data are represented with the origin at the
//   bottom left corner.

   }
   }
   fclose(canal);

   return(0);
}

// ASCII I/O routines
DataStructure ping_ascii( char *filename)
/*
   PURPOSE: To determine the geometry of a multicolumn ASCII file
*/
{
   DataStructure output=_DataStructure;
   int status=0;
   int skip;
   int Ncol,Nrow;

   status=column_count_ascii(filename,&skip,&Ncol,&Nrow);
   if ( (status==-1) || (Ncol<1) ) return(output);

/* Columns are taken as dimv dimension */

   output.Reader=&read_ascii;
   output.dimx=Nrow;
   output.dimy=1;
   output.dimv=Ncol;
   output.dimz=1;

/* Exiting the output */

   return(output);
}

int read_ascii( char *filename, DataWindow select, Matrix *data)
/*
   PURPOSE: To define a simple reading routine for multicolum ASCII files
*/
{
   FILE *chan;
   char buffer[ASC_LEN];

   float output;

   int goahead;
   int status;
   int skip;
   int Ncol,Nrow;
   int it,ic;

   status=column_count_ascii(filename,&skip,&Ncol,&Nrow);
   if ( (status==-1) || (Ncol<1) )
   {
      return(-1);
   }

/* Now we read data */

   create_matrix(select.dimx,select.dimy,0,data);

   chan=fopen(filename,"rt");

/* Skipping the first lines */

   for (it=0,goahead=1;(it<skip)&&goahead;it++) line_column_count_ascii(chan,&goahead);

/* Reading the others */

   for (it=0;(it<select.ix0+select.dimx)&&goahead;it++)
   {
      for(ic=0;(ic<Ncol)&&goahead;ic++)
      {
         if (fscanf(chan,"%s",buffer)==EOF) goahead=0;
         if ((it>=select.ix0)&&(ic==select.iv0))
         {
            sscanf(buffer,"%f",&output);
            data->p[it]=(double)output;
         }
      }
   }
   fclose(chan);

   return(status);
}

int column_count_ascii( char *filename, int *skip, int *Ncol, int *Nrow)
{
   FILE *chan;
   const int Nrowmin=5;
   int goahead;
   int Ncol0;
   int ic;

   if ( (chan=fopen(filename,"rt")) == NULL ) return(-1);

/* First, determine how many lines to skip */

   *skip=-Nrowmin;
   Ncol0=line_column_count_ascii(chan,&goahead);
   for (ic=0;(ic<Nrowmin)&&goahead;ic++,(*skip)++)
   {
    *Ncol=line_column_count_ascii(chan,&goahead);
    if (*Ncol!=Ncol0) ic--;
    Ncol0=*Ncol;
   }

/* Then, determine total file lenght */

   rewind(chan);
   *Nrow=0;
   if (goahead==0) return(-1);

   for (ic=0;ic<*skip;ic++)
    line_column_count_ascii(chan,&goahead);

   for (*Nrow=0;goahead;(*Nrow)++)
    Ncol0=line_column_count_ascii(chan,&goahead);
   if (Ncol0==0) (*Nrow)--;

/* Closing file and finishing */

   fclose(chan);

   return(0);
}

int line_column_count_ascii( FILE *chan, int *goahead)
{
 char buffer[ASC_LEN];
 int Ncol;

 *goahead=read_line_ascii(chan,&Ncol,buffer);

 if (*goahead==0) return(0);
 else             return(Ncol);
}

int read_line_ascii(FILE *chan, int *Ncol, char *line) {
  int goahead,stop;
  int left=0;
  int point,len;
  
  *Ncol = 0;
  for (len=0,stop=0,goahead=1;stop<1;point=fgetc(chan)) {
    
    if (point==EOF)  {
      stop=1;
      goahead=0;
    }
    
    else if (point==0x0a)      {
      line[len]=(char)0;
      stop=1;
      
    }    else {
      line[len]=(char)point;
      len++;
      if (len>=ASC_LEN) return(0);
      if (point<0x21)	{
	if (left)	  {
	  left=0;
	  (*Ncol)++;
	}
      } else left=1;
    }
 }
  
  if (left) (*Ncol)++;
  fseek(chan,-sizeof(char),SEEK_CUR);
  
  return(goahead);
}

/*       int FFT( Matrix functionR, Matrix functionI, int sign)
 * replaced by:
 *       int FFT_signal( Signal* funcR, Signal* funcI, int sign ) {
 */

/*       int FFThorizontal( Matrix functionR, Matrix functionI, int sign)
 * replaced by:
 *       int FFThorizontal_signal( Signal* funcR, Signal* funcI, int sign ) {
 */

/*       void PFFThorizontal( int dima, int dimb, int dimc, int iy, 
 * 	         	      double **uinR, double **uinI, 
 * 		              double **uoutR, double **uoutI, int sign )
 * unchanged 
 */

/*       int FFTvertical( Matrix functionR, Matrix functionI, int sign)
 * replaced by:
 *       int FFTvertical_signal( Signal* funcR, Signal* funcI, int sign );
 */

/*       void PFFTvertical( int dima, int dimb, int dimc, int ix, 
 * 	         	    double **uinR, double **uinI, 
 * 	         	    double **uoutR, double **uoutI, int sign)l
 * unchanged
 */


int assign_convolution( Matrix f1, Matrix f2)
/*
   PURPOSE: The convolution of f1 with f2 is calculated, then stored in f2.
   f1 is conserved and f2 is lost.
   This routine uses FFT.
*/
{
   Matrix Rf1=_Matrix;
   Matrix If1=_Matrix;
   Matrix If2=_Matrix;

   double buffR,buffI;
   double norma;
   int status=0;
   int ix,iy;

   if ( (f1.dimx!=f2.dimx) || (f1.dimy!=f2.dimy) ) return(-1); // incompatible matrices
   norma=sqrt((double)(f1.dimx*f1.dimy));

   create_matrix(f1.dimx,f1.dimy,0,&Rf1);
   create_matrix(f1.dimx,f1.dimy,0,&If1);
   create_matrix(f1.dimx,f1.dimy,0,&If2);

   assign_matrix(f1,Rf1);

   FFT(Rf1,If1,-1);
   FFT(f2,If2,-1);

   for (iy=0;iy<Rf1.dimy;iy++)
   {
   for (ix=0;ix<Rf1.dimx;ix++)
   {
    buffR = Rf1.Signal[iy][ix] * f2.Signal[iy][ix]  - If1.Signal[iy][ix] * If2.Signal[iy][ix];
    buffI = Rf1.Signal[iy][ix] * If2.Signal[iy][ix] + If1.Signal[iy][ix] * f2.Signal[iy][ix];
    f2.Signal[iy][ix]  = norma*buffR;
    If2.Signal[iy][ix] = norma*buffI;
   }
   }

   FFT(f2,If2,+1);

   destroy_matrix(&Rf1);
   destroy_matrix(&If1);
   destroy_matrix(&If2);

   return(status);
}

/*            Derivative routines         */

int gradient( Matrix gx, Matrix gy)
/*
   PURPOSE: Common enter point for the derivative routines
*/
{
   int status;

   switch (DER_MODE)
   {
    case 0:
     status=gradient_FFT(gx,gy);
     break;
    case 1:
     status=gradient_naif(gx,gy);
     break;
    default:
     printf("Unrecognized derivative mode\n");
     return(-1);
   }

   return(status);
}

int gradient_naif( Matrix gx, Matrix gy)
/*
   PURPOSE: Gradient estimate by finite differences
*/
{
   Matrix gxI=_Matrix,gyI=_Matrix;

   double dxR,dxI,dyR,dyI;
   double aux,x,y;
   int ix,iy;

   if ((gx.dimx<1)||(gx.dimy<1)||(gx.dimx!=gy.dimx)||(gx.dimy!=gy.dimy)) return(-1);

   create_matrix(gx.dimx,gx.dimy,0,&gxI);
   create_matrix(gy.dimx,gy.dimy,0,&gyI);

   FFT(gx,gxI,-1);
   for (iy=0;iy<gx.dimy;iy++)
   {
    y=((double)iy)/((double)gx.dimy); //for dimy=1, y=0
    dyR=cos(2.*PI*y)-1.;
    dyI=sin(2.*PI*y); // and so dyR and dyI are zero too.
    for (ix=0;ix<gx.dimx;ix++)
    {
     x=((double)ix)/((double)gx.dimx);
     dxR=cos(2.*PI*x)-1.;
     dxI=sin(2.*PI*x);
     gy.Signal[iy][ix]=dyR*gx.Signal[iy][ix]-dyI*gxI.Signal[iy][ix];
     gyI.Signal[iy][ix]=dyR*gxI.Signal[iy][ix]+dyI*gx.Signal[iy][ix];
     aux=dxR*gx.Signal[iy][ix]-dxI*gxI.Signal[iy][ix];
     gxI.Signal[iy][ix]=dxR*gxI.Signal[iy][ix]+dxI*gx.Signal[iy][ix];
     gx.Signal[iy][ix]=aux;
    }
   }
   FFT(gx,gxI,1);
   FFT(gy,gyI,1);

/* Memory release and end */

   destroy_matrix(&gxI);
   destroy_matrix(&gyI);

   return(0);
}

int gradient_FFT( Matrix gx, Matrix gy)
/*
   PURPOSE: A FFT surrogate to produce a finite-difference derivative
   estimate, of finite displacement 1 pixel but centered around the basis
   point. We call this oddity a "half-pixel centered finite difference".
   It has many interesting properties, although a certain tendence to create
   some aliasing effects. However, it fits incredibly well the behaviour of
   many power spectra.
*/
{
   Matrix gxI=_Matrix,gyI=_Matrix;
   double aux,x,y;
   int ix,iy;

   if ((gx.dimx<1)||(gx.dimy<1)||(gx.dimx!=gy.dimx)||(gx.dimy!=gy.dimy)) return(-1);

   create_matrix(gx.dimx,gx.dimy,0,&gxI);
   create_matrix(gy.dimx,gy.dimy,0,&gyI);

   FFT(gx,gxI,-1);
   for (iy=0;iy<gx.dimy;iy++)
   {
    y=((double)iy)/((double)gx.dimy);
    if ((iy>0)&&(iy>gx.dimy/2)) y-=1.;
    y=2.*sin(PI*y); //y=0 for dimy=1
   for(ix=0;ix<gx.dimx;ix++)
   {
    x=((double)ix)/((double)gx.dimx);
    if ((ix>0)&&(ix>gx.dimx/2)) x-=1.;
    x=2.*sin(PI*x); //x=0 for dimx=1

    gy.Signal[iy][ix]=-y*gxI.Signal[iy][ix];
    gyI.Signal[iy][ix]=y*gx.Signal[iy][ix];
    aux=-x*gxI.Signal[iy][ix];
    gxI.Signal[iy][ix]=x*gx.Signal[iy][ix];
    gx.Signal[iy][ix]=aux;
   }
   }
   FFT(gx,gxI,1);
   FFT(gy,gyI,1);

/* Memory release and end */

   destroy_matrix(&gxI);
   destroy_matrix(&gyI);

   return(0);
}

int modgradient( Matrix data)
/*
   PURPOSE: A simple routine to obtain the modulus of the gradient.
   High gradients associated to masked points are filtered out.
*/
{
   Matrix gx=_Matrix,gy=_Matrix;
   Matrix filtered=_Matrix;
   double mean_gradient;
   int ix,iy;

   if ( (data.dimx<1) || (data.dimy<1) ) return(-1);

   initiate_matrix(data,&gx);
   create_matrix(data.dimx,data.dimy,0,&gy);
   gradient(gx,gy);

   for (iy=0;iy<data.dimy;iy++)
   {
   for (ix=0;ix<data.dimx;ix++)
   {
      data.Signal[iy][ix]=sqrt(gx.Signal[iy][ix]*gx.Signal[iy][ix]+gy.Signal[iy][ix]*gy.Signal[iy][ix]);
   }
   }

   destroy_matrix(&gx);
   destroy_matrix(&gy);

/* Filtering the neighborhood of the mask and the edges */

   create_matrix(data.dimx,data.dimy,1,&filtered);
   assign_matrix(data, filtered);

   // Exclude the edges:

      for (iy=0,ix=0;iy<data.dimy;iy++)           filtered.Mask[iy][ix] = MASKED;
      for (iy=0,ix=data.dimx-1;iy<data.dimy;iy++) filtered.Mask[iy][ix] = MASKED;
   if (data.dimy>1)
   {
      for (iy=0,ix=0;ix<data.dimx;ix++)           filtered.Mask[iy][ix] = MASKED;
      for (iy=data.dimy-1,ix=0;ix<data.dimx;ix++) filtered.Mask[iy][ix] = MASKED;
   }

   // Exclude the mask's neighbors:
   if (data.Mask!=NULL)
   {
      for (iy=1;iy<data.dimy-1;iy++)
      {
      for (ix=1;ix<data.dimx-1;ix++)
      {
         if (data.Mask[iy][ix]==MASKED)
         {
            filtered.Mask[iy-1][ix] = MASKED;
            filtered.Mask[iy+1][ix] = MASKED;
            filtered.Mask[iy][ix-1] = MASKED;
            filtered.Mask[iy][ix+1] = MASKED;
         }
      }
      }
   }

   if (MEANMASK)
   {
      // The mean value of unmasked points is assigned to the masked points
      mean_gradient = matrix_mean(filtered);
      for (ix=0; ix<filtered.dimx*filtered.dimy; ix++) if (filtered.pm[ix] == MASKED) filtered.p[ix] = mean_gradient;
   }

   assign_matrix(filtered, data);

   destroy_matrix(&filtered);

   return(0);
}

/*       int reconstruct( Matrix gx, Matrix gy){
 * replaced by:
 *       int reconstruct_signal( Signal* gx, Signal* gy){
 */

/*       int reconstruct_naif( Matrix gx, Matrix gy)
 * replaced by:
 *       int reconstruct_signal_naif( Signal* gx, Signal* gy )
 */

/*       int reconstruct_FFT( Matrix gx, Matrix gy);
 * replaced by:
 *       int reconstruct_signal_FFT( Signal* gx, Signal* gy );
 */

/***************************************************************************/
int exponent_filter( double expon, Matrix data)
/***************************************************************************/
{
   Matrix auxR=_Matrix,auxI=_Matrix;
   double x,y,f;
   int ix,iy;

   if ( (data.dimx<1)||(data.dimy<1) ) return(-1);

   initiate_matrix(data,&auxR);
   create_matrix(data.dimx,data.dimy,0,&auxI);

   FFT(auxR,auxI,-1);

   for (iy=0;iy<data.dimy;iy++)
   {
      y=((double)iy)/((double)data.dimy);
      if (iy>data.dimy/2) y-=1.;
      y=2.*sin(PI*y);
   for (ix=0;ix<data.dimx;ix++)
   {
      x=((double)ix)/((double)data.dimx);
      if (ix>data.dimx/2) x-=1.;
      x=2.*sin(PI*x);
      f=sqrt(x*x+y*y);
      if (f>MIN_NUM_RESOLUTION) f=pow(f,expon);
      else                        f=0.;
      auxR.Signal[iy][ix]*=f;
      auxI.Signal[iy][ix]*=f;
   }
   }

   FFT(auxR,auxI,1);
   assign_matrix(auxR,data);

   destroy_matrix(&auxR);
   destroy_matrix(&auxI);

   return(0);
}

/* Multifractal processing routines */

/***************************************************************************/
double scale_linear( int dimx, int dimy)
/***************************************************************************/
/*
   PURPOSE: to provide the resolution scale
*/
{
    return( fMin( 1./((double)dimx), 1./((double)dimy) ) );
}


/***************************************************************************/
double scale_wavelet_bak( int dimx, int dimy)
/***************************************************************************/
/*
   PURPOSE: to provide a wavelet-dependent resolution scale
*/
{
   Matrix wave=_Matrix;
   double sc_w;
   double norma=0.;
   int ix,iy;

   if (PUNCTUAL)
   {
      return(scale_linear(dimx,dimy));
   }
   else
   {
      generate_wavelet(dimx, dimy, 1., &wave, 0);
      if (dimy>1)
      {
         for (iy=0; iy<dimy; iy++)
         {
         for (ix=0; ix<dimx; ix++)
         {
            norma += fabs(wave.Signal[iy][ix]);
         }
         }
      }
      else
      {
         for (ix=0; ix<dimx; ix++) norma += wave.Signal[0][ix]*wave.Signal[0][ix];
      }
      destroy_matrix(&wave);

      sc_w = sqrt(norma/((double)(dimx*dimy)));
      if (dimy>1) sc_w /= sqrt(PI);

      return(sc_w);
   }
}

/***************************************************************************/
double scale_wavelet( int dimx, int dimy)
/***************************************************************************/
/*
   PURPOSE: to provide a wavelet-dependent resolution scale
*/
{
   double sc_l;

   sc_l = scale_linear(dimx,dimy);
   if (PUNCTUAL) return(sc_l);
   else          return(sc_l*SC_W[ORD_DER][WVINDEX]);
}


/***************************************************************************/
int generate_wavelet( int dimx, int dimy, double sc, Matrix *wave, int normalize_flag)
/***************************************************************************/
/*
   PURPOSE: To generate a wavelet of the appropriate type, order and
   degree of derivation to analyze the data. Notice that the wavelet
   itself will not be differenciated; the derivatives are applied on
   the data directly, to enhance spatial detection. However, the changes
   in the normalization factor due to the derivatives must be taken into
   account in the definition of the wavelet.
*/
{
   double D_space;
   double expon;
   double x,y,base;
   double norm;
   int ix,iy;

/* Early checks */

   if ( (dimx<1)||(dimy<1)||(sc<MIN_NUM_RESOLUTION) ) return(-1);

/* Initialization */

   if (dimy==1) D_space=1.; // one-dimensional signal
   else         D_space=2.; // two-dimensional signal

   create_matrix(dimx, dimy, 0, wave);
   expon = 0.5*(D_space+WVINDEX-1);

/* Defining wavelet without normalization */

   norm=0.;
   for (iy=0;iy<dimy;iy++)
   {
      y = (double)iy;
      if (iy>dimy/2) y -= (double)dimy;
      y /= sc*SC0[ORD_DER][WVINDEX];
      for (ix=0;ix<dimx;ix++)
      {
         x = (double)ix;
         if (ix>dimx/2) x -= (double)dimx;
         x /= sc*SC0[ORD_DER][WVINDEX];

         if (WVINDEX>0) base = pow(1.+x*x+y*y,-expon);
         else           base = exp(-0.5*(x*x+y*y));
         wave->Signal[iy][ix] = base;
         norm += fabs(base);
      }
   }
   norm *= pow(sc, -ORD_DER); // include the shift due to derivatives

/* Normalize the wavelet */

   if (normalize_flag)
   {
      for(iy=0;iy<dimy;iy++)
      {
      for(ix=0;ix<dimx;ix++)
      {
         wave->Signal[iy][ix] /= norm;
      }
      }
   }

/* Exiting */

   return(0);
}

/***************************************************************************/
MFList calculate_multifractal( Matrix data,  Matrix *expon)
/***************************************************************************/
/*
    PURPOSE: To obtain the singularity exponents associated to a given signal.
    In addition to the matrix expon, data are outputted in a summary record
*/
{
   MFList output=_MFList;
   Matrix *wave=NULL;
   Matrix dmasa=_Matrix;

   double *xx,*yy;
   double mean;
   double wp;
   double Nev;
   double sc;
   double a,b,corr;
   double Nbuen;
   double qs;

   int nomask;
   int ix,iy,ip;

/* Initialization */

   if (VERBOSE) printf("Multifractal calculations:\n");
   if (VERBOSE) printf("+ Allocating memory\n");

   wave = (Matrix *)calloc(WVPOINTS,sizeof(Matrix));
   for (ip=0;ip<WVPOINTS;ip++)
   {
      wave[ip].dimx=wave[ip].dimy=0;
   }

   xx = (double *)calloc(WVPOINTS,sizeof(double));
   yy = (double *)calloc(WVPOINTS,sizeof(double));

   create_matrix(data.dimx, data.dimy, (data.pm!=NULL || FORCEMASK)?1:0, &dmasa);
   create_matrix(data.dimx, data.dimy, (data.pm!=NULL || FORCEMASK)?1:0, expon);

/* We proceed with the standard calculations */

   if (VERBOSE) printf("+ Constructing the standard analysis mass\n");

   assign_matrix(data,dmasa);

   if (HOLDER==0) modgradient(dmasa);
   if (ORD_DER)   exponent_filter((double)ORD_DER,dmasa);

   output.sc0 = scale_wavelet(data.dimx, data.dimy);

   if (PUNCTUAL) // Singularities are computed punctually
   {
      if (VERBOSE) printf("+ Punctual projection\n");

/* Obtaining the mean value */

      mean=matrix_mean(dmasa);

/* Normalizing directly by mean and logarithmically by scale */

      for (iy=0; iy<data.dimy; iy++)
      {
      for (ix=0; ix<data.dimx; ix++)
      {
         nomask=1;
         if (dmasa.Mask!=NULL)
         {
            if (dmasa.Mask[iy][ix]==MASKED) nomask=0;
         }
         if (nomask)
         {
            if (fabs(dmasa.Signal[iy][ix]/mean) > MIN_NUM_RESOLUTION)
               expon->Signal[iy][ix] = log(fabs(dmasa.Signal[iy][ix])/mean) / log(output.sc0);
            else
               expon->Signal[iy][ix] = log(MIN_NUM_RESOLUTION)            / log(output.sc0);
            output.hmin = fMin(output.hmin, expon->Signal[iy][ix]);
            output.hmax = fMax(output.hmax, expon->Signal[iy][ix]);
         }
         else
         {
            expon->Mask[iy][ix]=MASKED;
            expon->Signal[iy][ix]=(double)HMAX; // conventional value
             // (even so these points are excluded from the analysis)
         }
      }
      }
      output.densgood=1.; // No regression is done, so the qualitiy of the MF fit cannot be estimated in punctual mode
   }
   else // Singularitites are computed at multiple scales
   {
      if (VERBOSE) printf("+ Wavelet projections:\n");
      qs = pow(WAV_RANGE, 1./((double)(WVPOINTS-1)));

      /* Calculating the wavelet projections */

      for (ip=0,sc=1.; ip<WVPOINTS; ip++,sc*=qs)
      {
         if (VERBOSE) printf("  - Projection %d of %d .",ip+1,WVPOINTS);

         xx[ip]=1./log(sc*output.sc0);
         if (VERBOSE) printf(".");

         generate_wavelet(dmasa.dimx,dmasa.dimy,sc,&(wave[ip]),1);
         if (VERBOSE) printf(".");

         assign_convolution(dmasa,wave[ip]);
         if (VERBOSE) printf(". complete\n");
      }

      /* Exponents are evaluated by log-log regression at any valid point */

      if (VERBOSE) printf("+ Computing exponents\n");

      Nev=Nbuen=0.;
      for (iy=0; iy<data.dimy; iy++)
      {
      for (ix=0; ix<data.dimx; ix++)
      {
         nomask=1;
         if (dmasa.Mask!=NULL)
         {
            if (dmasa.Mask[iy][ix]==MASKED) nomask=0;
         }
         if (nomask)
         {
            Nev+=1.;
            for (ip=0;ip<WVPOINTS;ip++)
            {
               wp = fabs(wave[ip].Signal[iy][ix]);
               if (wp>MIN_NUM_RESOLUTION) yy[ip]=log(wp);
               else yy[ip]=log(MIN_NUM_RESOLUTION);
               yy[ip]*=xx[ip];
            }
            linear_fit(WVPOINTS,xx,yy,&a,&b,&corr); // fits yy = a * xx + b
            if (fabs(corr)>GOODRHO) Nbuen+=1.;
            expon->Signal[iy][ix] = b;
            output.hmin = fMin(output.hmin,expon->Signal[iy][ix]);
            output.hmax = fMax(output.hmax,expon->Signal[iy][ix]);
         }
         else
         {
            expon->Mask[iy][ix]=MASKED;
            expon->Signal[iy][ix]=(double)HMAX; // conventional value
             // (even so these points are excluded from the analysis)
         }
      }
      }
      if (Nev>0) output.densgood = Nbuen/Nev;
      else       output.densgood = 0.;
   }

   if (VERBOSE)
   {
      printf("Singularities: minimum: %f and maximum: %f\n", output.hmin,output.hmax);
      if (PUNCTUAL==0) printf("Percentage of good regression points: %f %%\n", 100.*output.densgood);
   }

/* Memory release and final outputting */

   for (ip=0;ip<WVPOINTS;ip++) destroy_matrix(&(wave[ip]));
   free(wave);
   destroy_matrix(&dmasa);
   free(xx);
   free(yy);

   return(output);
}

/***************************************************************************/
double non_translational_correction( Matrix data)
/***************************************************************************/
/*
   PURPOSE: To evaluate the dependende with the scale of the estimated mean
   derivative, to obtain the scaling exponent associated to a lack of
   translational invariance in the system.
*/
{
   Matrix aux=_Matrix;
   double *xx,*yy;
   double sc,qs;
   double expshift;
   double b,corr;

   int ip;

/* Initialization */

   qs = pow(WAV_RANGE, 1./((double)(WVPOINTS-1)));

   xx = (double *)calloc(WVPOINTS,sizeof(double));
   yy = (double *)calloc(WVPOINTS,sizeof(double));

/* Calculating the wavelet projections */

   if (VERBOSE) printf("Obtaining the dependence of the derivative with the scale\n");

   for (ip=0,sc=1.; ip<WVPOINTS; ip++,sc*=qs)
   {
      if (VERBOSE) printf("+ Projection %d of %d...",ip+1,WVPOINTS);
      xx[ip]=log(sc);
      Gaussian_derivative(sc,data,&aux);
      yy[ip]=log(matrix_mean(aux));

      destroy_matrix(&aux);
      if (VERBOSE) printf("ok\n");
   }
   linear_fit(WVPOINTS,xx,yy,&expshift,&b,&corr);
   if (NON_TRANS==1) expshift-=1.;
   else if (NON_TRANS==2) expshift-=2.;
   else printf("Warning! Bugged NonTrans mode!\n");
   if (VERBOSE) printf("Exponent correction: %0.2f (regression coeff.): %0.2f\n", expshift,corr);

/* Memory release and end */

   free(yy);
   free(xx);
   return(expshift);
}

/***************************************************************************/
void Gaussian_derivative( double sc, Matrix data, Matrix *modg)
/***************************************************************************/
/*
   PURPOSE: Furnishes a simple implementation of the Gaussian derivative
*/
{
   Matrix wave=_Matrix;
   double base;
   double x,y;
   int dimx,dimy;
   int ix,iy;

/* Initialization */

   dimx=data.dimx;
   dimy=data.dimy;
   create_matrix(dimx,dimy,0,&wave);
   create_matrix(dimx,dimy,0,modg);

/* Creating the Gaussian versions of finite increments */

   if (NON_TRANS==1)
   {
      for (iy=0;iy<dimy;iy++)
      {
         y = (double)iy;
         if (iy>dimy/2) y -= (double)dimy;
         y /= sc;
      for (ix=0;ix<dimx;ix++)
      {
         x = (double)ix;
         if (ix>dimx/2) x -= (double)dimx;
         x /= sc;
         base=exp(-0.5*(x*x+y*y));
         wave.Signal[iy][ix]=x*base;
         modg->Signal[iy][ix]=y*base;
      }
      }
      assign_convolution(data,wave);
      assign_convolution(data,*modg);
   }
   else if (NON_TRANS==2)
   {
      for (iy=0;iy<dimy;iy++)
      {
         y = (double)iy;
         if (iy>dimy/2) y -= (double)dimy;
         y /= sc;
      for (ix=0;ix<dimx;ix++)
      {
         x = (double)ix;
         if (ix>dimx/2) x -= (double)dimx;
         x /= sc;
         base=exp(-0.5*(x*x+y*y));
         wave.Signal[iy][ix]=(x*x+y*y-1.)*base;
         modg->Signal[iy][ix]=data.Signal[iy][ix];
      }
      }
      assign_convolution(wave,*modg);
   }
   else printf("Warning! Bugged NonTrans mode!\n");

/* Obtaining the modgradient */

   for (iy=0;iy<dimy;iy++)
   {
   for (ix=0;ix<dimx;ix++)
   {
      modg->Signal[iy][ix]=sqrt(wave.Signal[iy][ix]*wave.Signal[iy][ix]+modg->Signal[iy][ix]*modg->Signal[iy][ix]);
   }
   }

/* Memory release and end */

   destroy_matrix(&wave);
}


/***************************************************************************/
void wavelet_of_choice( char *output) {
  /***************************************************************************/
  /* PURPOSE: To generate a character chain referring the used wavelet */
  
  switch (WVINDEX) {
  case 0:
    sprintf(output,"Gaussian wavelet");
    break;
  case 1:
    sprintf(output,"Lorentzian wavelet of exponent 0.5*D_space");
    break;
  default:
     sprintf(output,"Lorentzian wavelet of exponent 0.5*(D_space+%d)",WVINDEX-1);
     break;
  }
  if (ORD_DER) sprintf(output,"%s in %d-th derivative\n",output,ORD_DER);
}

/***************************************************************************/
void standard_summary( char *filename, MFList summary)
/***************************************************************************/
/*
   PURPOSE: To ouput a text file with the summary of MF calculations
*/
{
   FILE *chan;
   char anyname[MAX_NAME];

   wavelet_of_choice(anyname);
   chan=fopen(filename,"wt");
   fprintf(chan, "Data filename: %s\n", summary.filename);
   if (summary.iv>0) fprintf(chan, "Data component: %d\n", summary.iv-1);
   if (summary.iz>0) fprintf(chan, "Data slice: %d\n", summary.iz-1);
   fprintf(chan, "Processing wavelet: %s\n", anyname);
   fprintf(chan, "Experimental singularity range: (%0.2f,%0.2f)\n", summary.hmin, summary.hmax);
   if (PUNCTUAL==0)
      fprintf(chan, "Percertage of points for which muCF is well verified: %0.2f%%\n", 100.*summary.densgood);
   else
      fprintf(chan, "The degree of muCF validity cannot be estimated by punctual analysis.\n");
   fclose(chan);
}

/***************************************************************************/
void singularity_histogram( Matrix expon, double *h, double *histoh, double sc0)
/***************************************************************************/
/*
   PURPOSE: Generates the histogram associated to singularity exponents
*/
{
   double auxdist, hshift;
   int ih, ihmode;

   histogram(expon, h, histoh, HMIN, HMAX);

   /* We numerically ensure the theoretical translational invariance.
    * A singularity shift is performed, if required.
    * Explanation: 
    * We define "translational invariance" as the scale independence of the first
    * order moment, i.e., the expected value, of the normalized multifractal
    * measure: <T_r> ~ r^\tau_1 with \tau_1 = 0, which is a very reasonable
    * requirement for any physical signal. This is explicitly assumed in the
    * punctual mode, when we have divided the T_r(x) by its mean to ensure that
    * the prefactor a(x) in the multifractal scaling, T_r(x) = a(x) r^h(x) + o(x),
    * has small influence, also implicitly assumed in the multiscale regression
    * mode, when wavelet projections are normalized. On the other hand, the \tau_p
    * functional is the Legendre-transform of the reduced singularity spectrum D(h),
    * i.e., \tau_p = {ph - D(h)} so \tau_1 = {h - D(h)} = 0 implies that the
    * spectrum is tangent to the f(h)=h line. However, numerical deviations can be
    * expected, since numerically computed means are not always translational
    * invariant. If this is the case, we adapt the estimated spectrum to satisfy 
    * this property.                    */
   if (NON_TRANS==-1) {
     
    ihmode = 0;
    for (ih=0; ih<NBOX; ih++) if (histoh[ih] > histoh[ihmode]) ihmode = ih;
    
    hshift = h[ihmode]/histoh[ihmode]; // initial guess
    for (ih=ihmode; ih>=0; ih--)
    {
     auxdist = (histoh[ih] > 0) ? (h[ih]/histoh[ih]) - (-log(histoh[ih]/histoh[ihmode])/log(sc0)) : auxdist;
     hshift = fMin(hshift, auxdist);
    }

    for (ih=0; ih<NBOX; ih++) h[ih] -= hshift*histoh[ih];

    if (VERBOSE) printf("To satisfy translational invariance, a singularity shift of %0.2f is required\n", -hshift);
   }
}

/*       void histogram( Matrix data, double *val, double *histo, double min,
 *                       double max );
 * replaced by:
 *       int compute_histo( dimx, dimy, s, mms, hist, mask );
 * with
 *       double *mms, **s;
 *       Histo *hist;
 *       char **mask 
 *       int dimx=data.dimx, dimy=data.dimy;
 *       hist->val=val, hist->h=histo; 
 *       mms[0]=min, mms[1]=max, mms[2]=max-min;
 *       s=data.data, mask=data.mask;
 */

/***************************************************************************/
void histo_record( char *filename, double *val, double *histo, double min, double max)
/***************************************************************************/
/* PURPOSE: Recording the histogram of values of val */
{
   FILE *chan;
   double d;
   double eff0, eff1;
   int ih;

   /* Initalization */
   d=(max-min)/((double)NBOX);
   
   /* Structure of the file: 1st column: uniform representative of val
    *                        2nd column: average representative of val
    *                        3rd column: number of cases
   */

   chan=fopen(filename,"wt");
   for (ih=0;ih<NBOX;ih++)  {
    eff0=min+(0.5+(double)ih)*d;
    if (histo[ih]>0.5) eff1 = val[ih]/histo[ih]; // there is at least one case
    else               eff1 = eff0;
    fprintf(chan,"%f  %f  %f\n",eff0,eff1,histo[ih]);
   }
   fclose(chan);
}

/***************************************************************************/
void histoh_adapt( double sc0, double sc1, double *h, double *histoh ) {
  /* replaced by:
   *       int scale_histo( Histo *h0, double sc0, Histo *h1,
   *       double sc1 );
   */
  /***************************************************************************/
  /* Adapt singularity histograms coming from data with different resolutions
   * sc1 to a common basic resolution sc0. */
  double hmode;
  double Nnew;
  int ih;
  
  /* Obtain the histogram mode */
  hmode = 0.;
  for (ih=0;ih<NBOX;ih++) hmode=fMax(hmode,histoh[ih]);
  
  /* We now modify the histogram taking the following into account:
   * histoh(h)= hmode * r^{d-D(h)}, and now r must change from sc1 (that
   * of data) to sc0 (the new basic resolution). The average h associated
   * to the bin box must also be changed accordingly. */
  
  for( ih=0; ih<NBOX; ih++ ) 
    if(histoh[ih] > 0.5) { // there is at least one case
      
      Nnew = hmode * exp(log(histoh[ih]/hmode) * log(sc0) / log(sc1));
      /* Notice that, even if in order to propose this rule for the scale
       * transformation we have implicitly introduced the existence of the
       * singularity spectrum, the transformation rule is more general and
       * would be valid for any quantity exhibiting a general power-law scaling
       * in the histogram. In addition, notice that if sc0=sc1, Nnew=histo[ih].
       * So that, this routine does not induce any change in the data if they are
       * completely homogeneous in size. Notice also that, as the dependence in
       * sc0, sc1 is logarithmic, the induced change is very small if data do
       * not differ in size by a considerable amount.   */
      
      /* The values are now updated */
      h[ih] *= Nnew / histoh[ih];
      histoh[ih] = Nnew;
    }
  
  return OK;
} // end of histoh_adapt

/***************************************************************************/
void Dh_record( char *filename, double sc0, double *h, double *histoh) {
/***************************************************************************/
   FILE *chan;
   double *histo_r,*h_r;
   double *Dh,*errDh;
   double prob,dp,Nev;
   double maxprob;
   double cump;
   const double min_ev=100;
   const double Ks=3; // 3 sigmas corresponds to a 99.7% confidence level
   int Nh;
   int ip,ir;

/* Inicialization */

   histo_r = (double *)calloc(NBOX,sizeof(double));
   h_r     = (double *)calloc(NBOX,sizeof(double));

/* Filtering histogram to avoid low-probability distorsions */

   cump=0.;
   ir=0;
   for (ip=0;ip<NBOX;ip++)
   {
    cump+=histoh[ip];
    histo_r[ir]+=histoh[ip]*histoh[ip];
    h_r[ir]+=h[ip];
    if(cump>min_ev)
    {
     histo_r[ir]/=cump;
     h_r[ir]/=cump;
     cump=0.;
     ir++;
    }
   }

   if (cump>0.)
   {
    histo_r[ir]/=cump;
    h_r[ir]/=cump;
    ir++;
   }
   Nh=ir;

/* The number of exponents is now known; we reservate memory accordingly */

   Dh    = (double *)calloc(Nh,sizeof(double));
   errDh = (double *)calloc(Nh,sizeof(double));

/* Producing the error bars */

// First, we compute the total number of events

   Nev=0.;
   for (ip=0;ip<Nh;ip++) Nev+=histo_r[ip];

/*
   The confidence range is taken as +- Ks sigmas in the distribution
   of probability boxes, which is a renormalized binomial. In this routine
   we directly propagate the confidence range to the D(h) calculation
*/

   for (ip=0;ip<Nh;ip++)
   {
    prob=histo_r[ip]/Nev; // probability of that interval
    if (prob>MIN_NUM_RESOLUTION) dp=Ks*sqrt((1.-prob)/(prob*Nev)); // Ks sigmas in the prob. distribution
    else dp=0.;
/*
   The error bar is constructed by propagation; as the propagated interval is
   asymetric, we take the bar as the maximum of the two distances to the
   central value. This is given by the lower bound in the logarithm
*/
    if (dp>=1.) errDh[ip]=1.;
    else errDh[ip]=log(1.-dp)/log(sc0);
   }

/* Finding and normalizing by the mode */

   maxprob=histo_r[0];
   for (ip=1;ip<Nh;ip++) maxprob=fMax(maxprob,histo_r[ip]);
   if (maxprob>MIN_NUM_RESOLUTION) for(ip=0;ip<Nh;ip++) histo_r[ip]/=maxprob;

/* Evaluating experimental reduced D(h) (red{D(h)} = D(h)-D_space */

   for (ip=0;ip<Nh;ip++)
   {
    if (histo_r[ip]>MIN_NUM_RESOLUTION) Dh[ip]=-log(histo_r[ip])/log(sc0);
    else Dh[ip]=-log(MIN_NUM_RESOLUTION)/log(sc0);
   }

/* Recording the results in an appropriate text file */

/* First column: h; second: reduced D(h); third: associated error bar */

   chan=fopen(filename,"wt");
   for (ip=0;ip<Nh;ip++) fprintf(chan,"%0.3f  %0.3f  %0.3f\n",h_r[ip],Dh[ip],errDh[ip]);
   fclose(chan);

/* Memory release and exiting */

   free(Dh);
   free(errDh);
   free(histo_r);
   free(h_r);
}
