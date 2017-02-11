/********************************************************/
/*                                                      */
/*    parse_args.c. - Version tou 2 Dekembriou 2004     */
/*                                                      */
/* Note that in the following the main modifications    */
/* with regard to the original parse functions are      */
/* designed between MODIF and END MODIF anchors         */
/*                                                      */
/********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Personnal libraries */
#include <utils.h>
#include <utl_alloc.h>
#include <mf_parse.h>


#ifndef M_PI
#define M_PI   3.1415926535997932
#endif

/* Global variables passed to the functions */


/* parsing all other parameters of the analysis
 * See file const-msm.h for description  */
extern int WAV_BASE;
extern int DER_WAV_BASE;
extern int ANALYSIS;

extern double DISP_ETA;

extern double D0[3]; 

extern double *prob_levi; 
extern double *prob_exp; 
extern double m_th[2]; 
extern double acprob;  




/***************************************************************************/
int parse_fractal1D( int in0, int siflag, int *deflag, 
		     ParseArg *atom ) {
  /***************************************************************************/
  
  int in=in0;

  return(in);
}

/***************************************************************************/
int parse_fractal2D( int in0, int siflag, int *deflag, 
		     ParseArg *atom ) {
  /***************************************************************************/
  int in=in0;

  //   Argument DIM1. Type 0: Flag
  strcpy(atom[in]argname,"-dim1");
  sprintf(atom[in].explain,"%s\n %s : %s %d\n",
	  "\nMULTIFRACTAL VARIABLES\n======================",atom[in].argname,
	  "Flag. If enabled, analysis is carried out in 1D."
	  "\n  Default:",DIM1);
  if(siflag)    {
    atom[in].type=3;
    atom[in].var.i=&(p_frac->dim1);
    atom[in].flag=deflag;
  }  else {
    atom[in].type=0;
    atom[in].flag=&(p_frac->dim1);
  }
  atom[in].num=1; in++;
	
   //   Argument THETAU. Type 2: Double
  strcpy(atom[in].argname,"-theta");
  strcpy(atom[in].valname,"angle_1D");
  sprintf(atom[in].explain," %s : %s %g\n",atom[in].argname,
	  "Multifractal variable. Direction for 1D analysis. Default:",THETAU);
  if(siflag)    {
    atom[in].type=5;
    atom[in].flag=deflag;
  }  else atom[in].type=2;
  /*  MODIF  */ atom[in].num=1; /*  END MODIF  */
  atom[in].var.d=&THETAU;
  atom[in].extval[0].d=0;
  atom[in].extval[1].d=M_PI/2.;
  in++;

  //   Argument HOLDER . Type 0: Flag
  strcpy(atom[in].argname,"-hold");
  sprintf(atom[in].explain," %s : %s\n",atom[in].argname,
	  "Flag. It enables function processing (in opposition to measure"
	  "\n  analysis). Default:",FLAG_HOLDER);
  if(siflag)    {
    atom[in].type=3;
    atom[in].flag=deflag;
    atom[in].var.d=&(p_frac->flag_holder);
  }  else {
    atom[in].type=0;
    atom[in].flag=&(p_frac->flag_holder);
  }
  atom[in].num=1; in++;



  return(in);
}

/***************************************************************************/
int filparser_init_wtmm( int in, int siflag, int *deflag, 
			 ParseArg *atom )  {
  /***************************************************************************/
  
  if(SelectPar[ipar_method_wtmm]+SelectPar[ipar_wav_wtmm]+
     SelectPar[ipar_der_wtmm]+SelectPar[ipar_sup_wtmm]+
     SelectPar[ipar_scrange]+SelectPar[ipar_momq]+
     SelectPar[ipar_shift]+SelectPar[ipar_tauq] >= TRUE)
    sprintf(atom[in].explain,"\nWTMM PARAMETERS\n===============\n");
  
  //    Argument WTMMMethod. Type 0: Flag
  IF(SelectPar[ipar_method_wtmm]) {
    sprintf(atom[in].argname,"%s","-canon"); 
    sprintf(atom[in].explain," %s : %s %d\n", atom[in].argname,
	    "Flag. If enabled, the method used in WTMM scheme to approximate the"
	    "\n  Legendre Transform is the canonical approach."
	    " Default:",FLAG_METHODWTMM); 
    atom[in].type=0;
    atom[in].flag= &(p_wtmm->flag_methodwtmm); // WTMMMethod;
    atom[in].num=1;    in++;
  }
  
  //    Argument FLAG_SUPWTMMM. Type 0: Flag
  IF(SelectPar[ipar_sup_wtmm]) {
    sprintf(atom[in].argname,"%s","-sup_wtmm"); 
    sprintf(atom[in].explain," %s : %s %d\n", atom[in].argname,
	    "Flag. If enabled,  Default:",FLAG_SUPWTMM); 
    atom[in].type=0;
    atom[in].flag= &(p_wtmm->flag_supwtmm); 
    atom[in].num=1;    in++;
  }

  
    //   Argument SCmin & SCmax. Type 2: Double
    strcpy(atom[in].argname,"-regsc");
    //  DEFFLAG(atom[in].argname);
    strcpy(atom[in].valname,"sc_min sc_max");
    sprintf(atom[in].explain," %s : %s [%g,%g]\n",atom[in].argname,
	    "Range of scales used to perform the estimation (through regression) of the"
	    "\n  multifractal exponents in WTMM scheme. Default:",SCMIN, SCMAX);
    if(siflag)    {
      atom[in].type=5;
      atom[in].flag=deflag;
    }  else atom[in].type=2;
    atom[in].var.d= &(p_wtmm->scmin); 
    atom[in+1].var.d= &(p_wtmm->scmax); 
    atom[in].extval[0].d=atom[in+1].extval[0].d=0.1;
    atom[in].extval[1].d=atom[in+1].extval[1].d=1000.;
    atom[in].num=2;    in+=2;
  
  IF(SelectPar[ipar_momq]) {
    //   Arguments MinMoment & MaxMoment. Type 2: Double
    strcpy(atom[in].argname,"-q");
    DEFFLAG(atom[in].argname);
    strcpy(atom[in].valname,"min_moment max_moment");
    sprintf(atom[in].explain," %s : %s [%g,%g]\n",atom[in].argname,
	    "Range of moments used to perform the estimation of the multifractal exponents"
	    "\n  in WTMM scheme. Default:", MINMOMENT, MAXMOMENT);
    if(siflag)    {
      atom[in].type=5;
      atom[in].flag=deflag;
    }  else atom[in].type=2;
    atom[in].var.d= &(p_wtmm->minmom); 
    atom[in+1].var.d= &(p_wtmm->maxmom); 
    atom[in].extval[0].d=atom[in+1].extval[0].d=-100.;
    atom[in].extval[1].d=atom[in+1].extval[1].d=100.;
    atom[in].num=2;    in+=2;
    
    //   Argument QStep. Type 2: Double
    strcpy(atom[in].argname,"-dq");
    DEFFLAG(atom[in].argname);
    strcpy(atom[in].valname,"moment_step");
    sprintf(atom[in].explain," %s : %s %g\n",atom[in].argname,
	    "Moment step of the WTMM analysis. If moment_step is choosed as 0, then a "
	    "\n  list of (irregularly spaced) moments will be used instead of the range"
	    "\n  of moments [min_moment, max_moment] (see variable qDefArray)."
	    "\n  Default:",DQ);
    if(siflag)    {
      atom[in].type=5;
      atom[in].flag=deflag;
    }    else atom[in].type=2;
    atom[in].var.d= &(p_wtmm->qstep); 
    atom[in].extval[0].d=0.,    atom[in].extval[1].d=10.;
    atom[in].num=1;    in++;
  }

  //   Argument ShiftSpectrum. Type 2: Double
  IF(SelectPar[ipar_shift]) {
    strcpy(atom[in].argname,"-cdh");
    strcpy(atom[in].valname,"shift_spectrum");
    sprintf(atom[in].explain," %s : %s %g\n",atom[in].argname,
	    "Constant to add to the values of the WTMM spectrum when storing it."
	    "\n  Default: ",SHIFTSPEC);
    if(siflag)    {
      atom[in].type=5;
      atom[in].flag=deflag;
    }  else atom[in].type=2;
    atom[in].var.d= &(p_wtmm->shift); // ShiftSpectrum;
    atom[in].extval[0].d=-2.,   atom[in].extval[1].d=2.;
    atom[in].num=1;    in++;
  }
  
  //    Argument flagTauq. Type 0: flag  
  IF(SelectPar[ipar_tauq]) {
    sprintf(atom[in].argname,"%s","-tauq"); 
    sprintf(atom[in].explain," %s : %s %d\n",
	    atom[in].argname,
	    "Flag. If enabled, the multifractal exponents computed by the WTMM scheme"
	    "\n  are saved. Default:",FLAG_TAUQWTMM); 
    atom[in].type=0;
    atom[in].flag= &(p_wtmm->flag_tauqwtmm); 
    atom[in].num=1;    in++;
  }
  
  return(in);
}

/***************************************************************************/
void parse_fractal(int argc, char *argv[]){
  /***************************************************************************/
  char **olarg,**olval,**olexp;
  double **ptrvar,**ptrval;
  int **ptrvar_i,**ptrval_i;
  int **ptrflag;
  int *type, *olnumb;
  int lar;
  int flagv;
  int arglen;
  
  /* MODIF */
  int i, cur;
  int Narg0=50;  // Initialization value; it should be greater than 
                 // (but not necessarily equal to) the number or arguments
  /* END MODIF */  

  int in,Narg;

  /*         Defining the expected arguments            */
  //use alloc_parser	


  in=0;
	


//    Argument H1. Type 2: float
  sprintf(atom[in].argname,"%s","-h1");
  sprintf(atom[in].valname,"%s","max_sing");
  sprintf(atom[in].explain," %s : %s %0.2f\n", atom[in].argname,
	  "Maximum singularity in binomial MFs. Default:",H1); 
  atom[in].type=2;
  atom[in].var.d=&(p_frac->h1);
  atom[in].extval[0].d=-1.,  atom[in].extval[1].d=2.;
  atom[in].num=1; in++;
  
//   Argument HMIN. Type 2: Float
    strcpy(atom[in].argname,"-hmin");
    strcpy(atom[in].valname,"hmin");
    sprintf(atom[in].explain," %s : %s  %0.2f\n",atom[in].argname,
     "Multifractal variable. Minimum allowed exponent. Default:",HMIN);
    atom[in].type=2;
    atom[in].var.d=&(p_frac->hmin);
    atom[in].extval[0].d=-1000.,    atom[in].extval[1].d=1000.;
    atom[in].num=1; in++;

//   Argument HMAX. Type 2: Float
    strcpy(atom[in].argname,"-hmax");
    strcpy(atom[in].valname,"hmax");
    sprintf(atom[in].explain," %s : %s  %0.2f\n",atom[in].argname,
	    "Multifractal variable. Maximum allowed exponent. No autoadjusting. Default:",HMAX);
    atom[in].type=2;
    atom[in].var.d=&(p_frac->hmax);
    atom[in].extval[0].d=-1000.,    atom[in].extval[1].d=1000.;
    atom[in].num=1; in++;

//    Argument FROM_DH. Type 0: flag
    sprintf(atom[in].argname,"%s","-fromDh");
    sprintf(atom[in].explain," %s : %s\n",
	    atom[in].argname,
	    "Flag. If enabled, the program takes previously computed D(h) files"
	    "\n  and estimates the error from them.\n Default: DISABLED"); 
    atom[in].type=0;
    atom[in].flag=&(p_frac->flag_fromdh);    
    atom[in].num=1; in++;
    
    //    Argument GEO_MAP. Type 0: flag
    sprintf(atom[in].argname,"%s","-geomap");
    sprintf(atom[in].explain," %s : %s\n",
	    atom[in].argname,
	    "Flag. If enabled, the program tries to generate quality maps"
	    "\n  for each method, changing geometry but keeping the given MF type."
	    "\n  Default:",FLAG_GEOMAP); 
    atom[in].type=0;
    atom[in].flag=&(p_frac->flag_geomap);
    atom[in].num=1; in++;

//    Argument TYPE_MAP. Type 0: flag
    sprintf(atom[in].argname,"%s","-typemap");
    sprintf(atom[in].explain," %s : %s %d\n",
	    atom[in].argname,
	    "Flag. If enabled, the program tries to generate quality maps"
	    "\n  for each method, for fixed geometry (1x16384) and changing parameters"
	    "\n  in the given MF type.\n Default:",FLAG_TYPEMAP); 
    atom[in].type=0;
    atom[in].flag=&(p_frac->flag_typemap);
    atom[in].num=1; in++;
    
  //    Argument NSERIES. Type 1: integer
  sprintf(atom[in].argname,"%s","-N");
  sprintf(atom[in].valname,"%s","#series");
  sprintf(atom[in].explain," %s : %s %d\n",atom[in].argname,
	  "Number of series to be processed. Default:",NSERIES); 
  atom[in].type=1;
  atom[in].var.i=&(p_frac->nseries);
  atom[in].extval[0].i=1,  atom[in].extval[1].i=10000;
  atom[in].num=1; in++;


//    Argument LEFF. Type 1: integer

  sprintf(atom[in].argname,"%s","-dim");
  sprintf(atom[in].valname,"%s","length");
  sprintf(atom[in].explain," %s : %s %d\n",atom[in].argname,
	  "Size of series to be processed. Default:",LEFF); 
  atom[in].type=1;
  atom[in].var.i=&(p_frac->leff);
  atom[in].extval[0].i=256,  atom[in].extval[1].i=65536;
  atom[in].num=1; in++;

   //    Argument FLAG_SAVEFLOAT. Type 0: flag
  sprintf(atom[in].argname,"%s","-float");
  sprintf(atom[in].explain," %s : %s %d\n",atom[in].argname,
	  "Flag. If enabled, the data series are saved in float format instead"
	  "\n  of double. Default:",FLAG_SAVEFLOAT); 
  atom[in].type=0;
  atom[in].flag=&(p_frac->flag_savefloat);
  atom[in].num=1; in++;
  
  //    Argument LMAX. Type 1: integer
  sprintf(atom[in].argname,"%s","-dim");
  sprintf(atom[in].valname,"%s","size");
  sprintf(atom[in].explain," %s : %s %d\n",atom[in].argname,
	  "Linear size for generated series. It will be rounded to"
	  "\n  the least power of 2 greater than this value. Default:",LMAX); 
  atom[in].type=1;
  atom[in].var.i=&(p_frac->lmax);
  atom[in].extval[0].i=2,  atom[in].extval[1].i=1000000;
  atom[in].num=1; in++;

  //    Argument OUTRES. Type 1: integer
  sprintf(atom[in].argname,"%s","-outres");
  sprintf(atom[in].valname,"%s","#res_levels");
  sprintf(atom[in].explain," %s : %s %d\n",atom[in].argname,
	  "Number of dyadic resolutions to be smoothened. Default:",
	  OUTRES); 
  atom[in].type=1;
  atom[in].var.i=&(p_frac->outres);
  atom[in].extval[0].i=0,  atom[in].extval[1].i=10;
  atom[in].num=1; in++;

  //    Argument FLAG_INVTRANS. Type 0: flag
  sprintf(atom[in].argname,"%s","-invtrans");
  sprintf(atom[in].explain," %s : %s %d\n",atom[in].argname,
	  "Flag. If enabled, the program produces translational invariant multifractals."
	  "\n  Default:",FLAG_INVTRANS); 
  atom[in].type=0;
  atom[in].flag=&(p_frac->flag_invtrans);
  atom[in].num=1;  in++;

  //    Argument DIM_SPACE. Type 1: integer
  sprintf(atom[in].argname,"%s","-d_space");
  sprintf(atom[in].valname,"%s","dimension");
  sprintf(atom[in].explain," %s : %s %d\n", atom[in].argname,
	  "Dimension of the embedding space. Default:",
	  DSPACE); 
  atom[in].type=1;
  atom[in].var.i=&(p_frac->dim_space);
  atom[in].extval[0].i=DIM1D,  atom[in].extval[1].i=DIM2D;
  atom[in].num=1; in++;

  //    Argument TYPE_MF. Type 1: integer
  sprintf(atom[in].argname,"%s","-type");
  sprintf(atom[in].valname,"%s","mult_type");
  sprintf(atom[in].explain," %s : %s %d\n", atom[in].argname,
	  "Type of multifractal to be generated.\n   "
	  "0: Log-Poisson\n   "
	  "1: Log-Normal\n   "
	  "2: Log-Levi\n   "
	  "3: Binomial (bi-fractal) \n   "
	  "4: Monofractal\n  Default:",TYPE_MF); 
  atom[in].type=1;
  atom[in].var.i=&(p_frac->type_mfsim);
  atom[in].extval[0].i=TYPLOGPOISSON,  atom[in].extval[1].i=TYPMONOFRACTAL;
  atom[in].num=1; in++;

  //    Argument HINF. Type 2: double
  sprintf(atom[in].argname,"%s","-hinf");
  sprintf(atom[in].valname,"%s","min_sing");
  sprintf(atom[in].explain," %s : %s %g\n", atom[in].argname,
	  "Most singular exponent. Valid for log-Poisson, mono and binomials. "
	  "Default:",HINF); 
  atom[in].type=2;
  atom[in].var.d=&(p_frac->hinf);
  atom[in].extval[0].d=-1.,  atom[in].extval[1].d=0.;
  atom[in].num=1; in++;

  //    Argument CODINF. Type 2: double
  sprintf(atom[in].argname,"%s","-Codinf");
  sprintf(atom[in].valname,"%s","min_sing_cod");
  sprintf(atom[in].explain," %s : %s %g\n", atom[in].argname,
	  "Most singular codimension. Valid for log-Poisson, mono and binomials. "
	  "Default:",CODINF); 
  atom[in].type=2;
  atom[in].var.d=&(p_frac->codinf);
  atom[in].extval[0].d=0.,  atom[in].extval[1].d=DIM2D;
  atom[in].num=1; in++;

  //    Argument H1. Type 2: double
  sprintf(atom[in].argname,"%s","-h1");
  sprintf(atom[in].valname,"%s","max_sing");
  sprintf(atom[in].explain," %s : %s %g\n",atom[in].argname,
	  "Maximum singularity in binomial MFs. Default:",H1); 
  atom[in].type=2;
  atom[in].var.d=&(p_frac->h1);
  atom[in].extval[0].d=-1.,  atom[in].extval[1].d=2.;
  atom[in].num=1;    in++;

  //    Argument MU. Type 2: double
  sprintf(atom[in].argname,"%s","-mu");
  sprintf(atom[in].valname,"%s","sing_av");
  sprintf(atom[in].explain," %s : %s %g\n",atom[in].argname,
	  "Singularity mean. Valid for log-Normal and log-Levi. Default:",MU); 
  atom[in].type=2;
  atom[in].var.d=&(p_frac->mu);
  atom[in].extval[0].d=-5.,  atom[in].extval[1].d=5.;
  atom[in].num=1; in++;

  //    Argument SIGMA. Type 2: double
  sprintf(atom[in].argname,"%s","-sigma");
  sprintf(atom[in].valname,"%s","disp_sing");
  sprintf(atom[in].explain," %s : %s %g\n",atom[in].argname,
	  "Singularity dispersion. Valid for log-Normal and log-Levi."
	  "\n  Default:",SIGMA); 
  atom[in].type=2;
  atom[in].var.d=&(p_frac->sigma);
  atom[in].extval[0].d=0.,  atom[in].extval[1].d=5.;
  atom[in].num=1; in++;

  //    Argument TCH. Type 2: double
  sprintf(atom[in].argname,"%s","-max_disp");
  sprintf(atom[in].valname,"%s","#sigmas");
  sprintf(atom[in].explain," %s : %s %g\n",atom[in].argname,
	  "Singularity range, expressed in sigmas. Valid for log-Normal and log-Levi."
	  "\n  Default:",TCH); 
  atom[in].type=2;
  atom[in].var.d=&(p_frac->tch);
  atom[in].extval[0].d=1.,  atom[in].extval[1].d=20.;
  atom[in].num=1; in++;

  //    Argument ALPHA. Type 2: double
  sprintf(atom[in].argname,"%s","-alpha");
  sprintf(atom[in].valname,"%s","exponent");
  sprintf(atom[in].explain," %s : %s %g\n",atom[in].argname,
	  "Valid for log-Levi only: log-Levi exponent. Default:",
	  ALPHA); 
  atom[in].type=2;
  atom[in].var.d=&(p_frac->alpha);
  atom[in].extval[0].d=0.,  atom[in].extval[1].d=2.;
  atom[in].num=1; in++;

  //    Argument DENS. Type 2: double
  sprintf(atom[in].argname,"%s","-dens");
  sprintf(atom[in].valname,"%s","value");
  sprintf(atom[in].explain," %s : %s %g\n",atom[in].argname,
	  "Monofractal density; valid for monofractals only Default:",DENS); 
  atom[in].type=2;
  atom[in].var.d=&(p_frac->dens);
  atom[in].extval[0].d=0.,  atom[in].extval[1].d=1.;
  atom[in].num=1; in++;

  //    Memory usage parameters (included in <FFT1D.c>)
  in=parsing_memory(in,0,0,olarg,olval,olexp,ptrvar,ptrval, 
		    ptrvar_i,ptrval_i,ptrflag,type, 
		    /*  MODIF */olnumb/*  END MODIF */);

  //   Derivative parameters (included in <derivacion_1D.c>)  
  in=parsing_derivacion(in,0,0,olarg,olval,olexp,ptrvar,ptrval,
			ptrvar_i,ptrval_i,ptrflag,type, 
			/*  MODIF */olnumb/*  END MODIF */);
  
  //   multifractal parameters (included in <multimf_1D.c>)  
  in=parsing_multimf_1D(in,0,0,olarg,olval,olexp,ptrvar,ptrval,
			     ptrvar_i,ptrval_i,ptrflag,type, 
			     /*  MODIF */olnumb/*  END MODIF */);

  /*  MODIF */  
  // Wavelet Transform Modulus Maxima method parameters
#ifdef FLAG_WTMM
  in=parsing_wtmm( in, 0, 0, olarg, olval, olexp, ptrvar, ptrval, 
		   ptrvar_i, ptrval_i, ptrflag, type, olnumb);
#endif
  /* END MODIF */  
  
  //    Argument WAV_BASE. Type 1: integer
  sprintf(atom[in].argname,"%s","-wv_basis");
  sprintf(atom[in].valname,"%s","choice");
  sprintf(atom[in].explain,"%s\n %s : %s %d\n",
	  "\nPARAMETERS DEFINING THE REPRESENTATION BASIS\n"
	  "===========================================",atom[in].argname,
	  "Wavelet of choice for the basis."
	  "\n    0: Gaussian wavelet"
	  "\n    1: Lorentzian wavelet"
	  "\n    2: Diagonal Haar\n  Default:",WAVBASE); 
  atom[in].type=1;
  /*  MODIF  */ /*  END MODIF  */
  atom[in].var.i=&(p_frac->wavbase);
  atom[in].extval[0].i=0,  atom[in].extval[1].i=2;
  atom[in].num=1;  in++;

  //    Argument DER_WAV_BASE. Type 1: integer
  sprintf(atom[in].argname,"%s","-der_wv_basis");
  sprintf(atom[in].valname,"%s","order");
  sprintf(atom[in].explain," %s : %s %d\n",atom[in].argname,
	  "Order of the derivatives in the wavelet. Default:",DERWAVBASE); 
  atom[in].type=1;
  atom[in].var.i=&(p_frac->derwavbase);
  atom[in].extval[0].i=0,  atom[in].extval[1].i=2;
  atom[in].num=1;  in++;
  
  //    Argument ANALYSIS. Type 0: Flag  
  sprintf(atom[in].argname,"%s","-analysis");
  sprintf(atom[in].explain,"%s\n %s : %s\n",
	  "\nPARAMETER LAUNCHING THE ANALYSIS\n"
	  "========================================",atom[in].argname,
	  "If enabled, the program analyzes and obtains the singularity spectra"
	  "\n  of the data in course of generation thanks to the different methods: Histogram,"
	  "\n  Singularity Analysis"
	  ".  Default: DISABLED"); 
  atom[in].type=0;
  /*  MODIF  */ atom[in].num=1; /*  END MODIF  */
  atom[in].flag=&ANALYSIS;
  in++;
  /* END MODIF */  

  //    Argument NBOX. Type 1: integer
  sprintf(atom[in].argname,"%s","-Nbin");
  sprintf(atom[in].valname,"%s","#bins");
  sprintf(atom[in].explain,"%s\n %s : %s %d\n",
	  "\nPARAMETER CONCERNING HISTOGRAM ANALYSIS\n"
	  "========================================",atom[in].argname,
	  "Number of histogram bins. Default:",NBOX); 
  atom[in].type=1;
  atom[in].var.i=&(p_frac->nbox);
  atom[in].extval[0].i=2,  atom[in].extval[1].i=16384;
  atom[in].num=1;  in++;

  /*             End of parameters definition      */

  /* initialisa parser */
  Narg=in;
  // use init_parser


  
  /*            Freeing memory before terminating        */
  // use free_parser
	
  /*                Termination              */
  if(flagv!=1) exit(-1);

}


/***************************************************************************/
int fracparse_update_wtmm() {
  /***************************************************************************/

  if(p_wtmm->shift == CRAZY) p_wtmm->shift=SHIFTSPEC;
  
  if (p_wtmm->minmom >= p_wtmm->maxmom) {
    p_wtmm->minmom = MINMOMENT;
    p_wtmm->maxmom = MAXMOMENT;
  }
  
  if(p_wtmm->qstep == CRAZY) p_wtmm->qstep = DQ;
  
  if(p_wtmm->scmin >= p_wtmm->scmax) {
    p_wtmm->scmin = SCMIN; 
    p_wtmm->scmax = SCMAX;
  }
  
  IFVERBOSE {
    if(QStep != 0.)
      WarningVV("Range of moments and moment step used to perform WTMM estimation:"
		"\n      [min_moment=%g, max_moment=%g] - moment_step=%g",
		p_wtmm->minmom, p_wtmm->maxmom, p_wtmm->qstep);
    else {
      int i;
      printf("\n List of moments used to perform WTMM estimation:\n [");
      for( i=0; i<nqDef; i++ ) printf(" %g",qDefArray[i]);
      printf(" ]");
    }
    WarningVV(" Range of (log)scales used to run regression for WTMM estimation:"
	      "\n     [log(sc_min)=%g, log(sc_max)=%g]", 
	      LOG(p_wtmm->scmin), LOG(p_wtmm->scmax));
  }
  
  return OK;
} // end of filparse_update_wtmm


#ifdef DEBUG
int display_variables() {

  WarningV("NSERIES=%d",NSERIES);
  WarningV("D_space=%d",D_space);
  WarningV("LMAX=%d",LMAX);
  WarningV("OUTRES=%d",OUTRES);
  WarningV("TYPE=%d",TYPE);
  WarningV("FLAG_INVTRANS=%d",FLAG_INVTRANS);

  WarningV("HINF=%f",HINF);
  WarningV("CODINF=%f",CODINF);
  WarningV("H1=%f",H1);
  WarningV("H0=%f",H0);
  WarningV("MU=%f",MU);
  WarningV("SIGMA=%f",SIGMA);
  WarningV("ALPHA=%f",ALPHA);
  WarningV("DENS=%f",DENS);
  WarningV("TCH=%f",TCH);

  WarningV("WAV_BASE=%d",WAV_BASE);
  WarningV("DER_WAV_BASE=%d",DER_WAV_BASE);
  WarningV("ANALYSIS=%d",ANALYSIS);
  WarningV("NBOX=%d",NBOX);

  WarningV("DISP_ETA=%f",DISP_ETA);
  WarningV("DH=%f",DH);
  WarningV("D0[0]=%f",D0[0]);
  WarningV("D0[1]=%f",D0[1]);
  WarningV("D0[2]=%f",D0[2]);

  WarningV("m_th[0]=%f",m_th[0]);
  WarningV("m_th[1]=%f",m_th[1]);
  WarningV("acprob=%f",acprob);

  WarningV("NPOINTS=%d",NPOINTS);
  WarningV("WAV=%d",WAV);
  WarningV("ORDDER=%d",ORDDER);
  WarningV("HOLDER=%d",HOLDER);
  WarningV("DIM1=%d",DIM1);

  WarningV("S0=%f",S0);
  WarningV("THETAU=%f",THETAU);

  return OK;
}

#endif
