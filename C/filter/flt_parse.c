/* ===================================
** flt_parse.c
** started on Thu Jan 18 11:01:48 2007 
** ===================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <utils.h>
#include <utl_alloc.h>
#include <utl_char.h>
#include <utl_operator.h>

#include <filter.h>
#include <flt_parse.h>


extern ParFILT *p_fil;
extern ParWAV *p_wav;
extern char SelectPar[];


/***************************************************************************/
ParFILT* filparse_alloc() {
  /***************************************************************************/
  ParFILT *par;
  
  if( (par = (ParFILT*)malloc(sizeof(ParFILT))) == NULL )
    return NULL;
  else
    return par;
} // end of filparse_alloc


/***************************************************************************/
int filparse_free( ParFILT *par ) {
  /***************************************************************************/
  Free( par );
  return OK;
} // end of filparse_free


/***************************************************************************/
ParFILT* filparse_create( ) {
/***************************************************************************/
  ParFILT *par;
  
  if( (par=filparse_alloc()) == NULL) return NULL;
  filparse_default(par);
  
  return par;
} // end of filparse_init


/***************************************************************************/
int filparse_default( ParFILT *par ) {
/***************************************************************************/

  if(SelectPar[ipar_mem] == TRUE) /* We have to select the flag '-mem' because
				   * by default FLAG_MEMORY is not set */
    par->flag_memory = FALSE; 
  else         /* We have to select the flag '-nomem' because by default
		* FLAG_MEMORY is set */
    par->flag_memory = TRUE;
  //  par->flag_nomem =  1-par->flag_memory;

  return OK;
} // end of filparse_default_fft


/***************************************************************************/
int filparse_init( int in0, int siflag, int *deflag, 
		   char **olarg, char **olval, char **olexp,
		   double **ptrvar, double **ptrval, 
		   int **ptrvar_i, int **ptrval_i,
		   int **ptrflag, int *type, int*olnumb ) {
  /***************************************************************************/
  int in=in0;

  /** DERIVATIVE PARAMETERS **/
  if(SelectPar[ipar_dermode]+SelectPar[ipar_freq] >= TRUE) 
    sprintf(olexp[in],"\nDERIVATIVE PARAMETERS\n=================\n");
  
  //    Argument mode_deriva. Type 1: Integer  
  IF(SelectPar[ipar_dermode]) {
    sprintf(olarg[in],"-dermode");
    strcpy(olval[in],"mode");
    sprintf(olexp[in]," %s : %s %d\n",olarg[in],
	    "Method for calculating derivatives"
	    "\n    0: Half-pixel FFT kernel derivative"
	    "\n    1: Direct 1-pixel increment"
	    "\n    2: Derivative interpolation\n  Default:",MODE_DERIVA);
    if(siflag)    {
      type[in]=4;
      ptrflag[in]=deflag;
    }  else type[in]=1;
    ptrvar_i[in]=&(p_fil->mode_deriva); 
    ptrval_i[in][0]=0;
    ptrval_i[in][1]=2;
    olnumb[in]=1; in++;
  }
  
  //    Argument mode_freq. Type 0: Flag  
  IF(SelectPar[ipar_freq]) {
    sprintf(olarg[in],"-freq");
    sprintf(olexp[in]," %s : %s\n", olarg[in],
	    "Flag. If enabled, the definition of the frequency vector is"
	    "\n        f_x=sin(PI*x) f_y=sin(PI*y),"
	    "\n  otherwise:"
	    "\n        f_x=x f_y=y,"
	    "\n  with (x,y) the coordinates in Fourier space. Default: DISABLED");
    if(siflag)    {
      type[in]=3;
      ptrflag[in]=deflag;
    }  else  type[in]=0;
    ptrflag[in]=&(p_fil->mode_freq); 
    olnumb[in]=1; in++;
  }
  
  /** FLAG_MEMORY PARAMETERS **/
  if(SelectPar[ipar_mem] >= TRUE)
    sprintf(olexp[in],"\nFLAG_MEMORY PARAMETERS\n=================\n");
  
  //    Argument FLAG_MEMORY. Type 0: Flag
  IF(SelectPar[ipar_mem]) {
    sprintf(olarg[in],"-memory");
    sprintf(olexp[in]," %s : %s\n", olarg[in],
	    "Flag. If enabled, the program makes a better use of memory (at the"
	    "\n  cost of longer processing times). Default: DISABLED");
    if(siflag)    {
      type[in]=3;
      ptrflag[in]=deflag;
    }  else type[in]=0;
    ptrvar_i[in]=&(p_fil->flag_memory); //&FLAG_MEMORY;
    olnumb[in]=1; in++;
  }
  
  //   Argument NFLAG_MEMORY. Type 0: Flag
  /*   IF(SelectPar[ipar_nmem]) { */
  /*     strcpy(olarg[in],"-nomem"); */
  /*     sprintf(olexp[in]," %s : %s\n",olarg[in], */
  /* 	    "Flag for no optimization of the memory (ex: images are locally copied)."); */
  /*     type[in]=0; */
  /*     ptrflag[in]=&(p_fil->flag_nomem);  */
  /*     olnumb[in]=1;   in++; */
  /*   } */
  
return in;
}


/***************************************************************************/
int wavparse_init( int in0, int siflag, int *deflag, 
		   char **olarg, char **olval, char **olexp,
		   double **ptrvar, double **ptrval, 
		   int **ptrvar_i, int **ptrval_i,
		   int **ptrflag, int *type, int*olnumb ) {
  /***************************************************************************/
  int in=in0;

  /** WAVELET VARIABLES **/ 
  if(SelectPar[ipar_wav]+SelectPar[ipar_der]+
     SelectPar[ipar_scrange]+SelectPar[ipar_theta] >= TRUE)
    sprintf(olexp[in],"\nWAVELET VARIABLES\n===============\n");
  
  //   Argument WAV. Type 1: Integer
  IF(SelectPar[ipar_wav]) {
    strcpy(olarg[in],"-wav");
    strcpy(olval[in],"wav_index");
    sprintf(olexp[in]," %s : %s %d\n", olarg[in],
	    "Wavelet to be used, either:"
	    "\n    -2   : Morlet"
	    "\n    -1   : Gaussian"
	    "\n     0   : Haar"
	    "\n   [1-3] : Lorentzian at exponent <wav_index>/2"
	    /* "\n   1: Lorentzian at exponent 0.5"
	       "\n   2: Lorentzian"
	       "\n   3: Lorentzian at 1.5.\n*/  
	    "\n Default:",WAV);
    if(siflag)    {
      type[in]=4;
      ptrflag[in]=deflag;
    }  else type[in]=1;
    ptrvar_i[in]=&(p_wav->wav);
    ptrval_i[in][0]=-2;
    ptrval_i[in][1]=4;
    olnumb[in]=1;    in++;
  }
  
  //   Argument ord_der. Type 1: Integer
  IF(SelectPar[ipar_der]) {
    strcpy(olarg[in],"-der");
    strcpy(olval[in],"deriv_order");
    sprintf(olexp[in]," %s : %s %d\n",olarg[in],
	    "Derivative order, used only when the wavelet is gaussian or lorentzian."
	    "\n  Default:", ORDDER);
    if(siflag)    {
      type[in]=4;
      ptrflag[in]=deflag;
    }  else type[in]=1;
    ptrvar_i[in]= &(p_wav->ord_der); 
    ptrval_i[in][0]=0;
    ptrval_i[in][1]=5;
    olnumb[in]=1;    in++;
  }
  
  IF(SelectPar[ipar_scrange]) {
    //   Argument WAVRANGE . Type 2: Float
    strcpy(olarg[in],"-scrng");
    strcpy(olval[in],"scale_range");
    sprintf(olexp[in]," %s : %s %f\n",olarg[in],
	    "Range of scales in wavelet analysis. Default:",WAVRANGE);
    if(siflag)	{
      type[in]=5;
      ptrflag[in]=deflag;
    }	else type[in]=2;
    ptrvar[in]=&(p_wav->wav_range);
    ptrval[in][0]=1.00001;
    ptrval[in][1]=1000.;
    in++;
	
    //   Argument MinScale. Type 2: Double
    strcpy(olarg[in],"-sc0");
    DEFFLAG(olarg[in]);
    strcpy(olval[in],"scale_0");
    sprintf(olexp[in]," %s : %s %g\n",olarg[in],
	    "Initial scale of analysis for wavelets. Default:",MINSCALE); //S0
    if(siflag)    {
      type[in]=5;
      ptrflag[in]=deflag;
    }  else type[in]=2;
    ptrvar[in]= &(p_wav->minscale); //&(p_wav->sc0)
    ptrval[in][0]=.33; 
    ptrval[in][1]=100.;
    olnumb[in]=1;    in++;
  
    //   Argument MaxScale. Type 2: Double
    strcpy(olarg[in],"-scmax");
    DEFFLAG(olarg[in]);
    strcpy(olval[in],"scale_max");
    sprintf(olexp[in]," %s : %s\n",olarg[in],
	    "Maximal scale of analysis. By default, if it is not given,  it is computed"
	    "\n automatically according to the lenght of the input signal (see options"
	    "\n  '-tsc' and '-rsc' below).");
    if(siflag) {
    type[in]=5;
    ptrflag[in]=deflag;
    }  else type[in]=2;
    ptrvar[in]= &(p_wav->maxscale); 
    ptrval[in][0]=1; 
    ptrval[in][1]=1000.;
    olnumb[in]=1;  in++;
    
    //   Argument NVoices. Type 1: Integer
    strcpy(olarg[in],"-nvoi");
    strcpy(olval[in],"no_voices");
    sprintf(olexp[in]," %s : %s %d\n",olarg[in],
	    "Number of voices per octave. The (multiplicative) scale step in WTMM scheme"
	    "\n  will be 2^<1/no_voices>. Default: ", NVOICES);
    if(siflag)   {
      type[in]=4;
      ptrflag[in]=deflag;
    }  else type[in]=1;
    ptrvar_i[in]=&(p_wav->nvoices);
    ptrval_i[in][0]=1;  
    ptrval_i[in][1]=100.;
    olnumb[in]=1;    in++;
    
    //   Argument ScTime. Type 1: Integer
    strcpy(olarg[in],"-tsc");
    strcpy(olval[in],"time_scale");
    sprintf(olexp[in]," %s : %s %g\n",olarg[in],
	    "Factor of convolution range in WTMM scheme, i.e.:\n"
	    "        [wavelet box] = #{-time_scale*nsc,...,time_scale*nsc}."
	    "\n  Default: ", SCTIME);
    if(siflag)    {
      type[in]=4;
      ptrflag[in]=deflag;
    }    else type[in]=1;
    ptrvar_i[in]= &(p_wav->sctime); 
    ptrval_i[in][0]=1;  
    ptrval_i[in][1]=100.;
    olnumb[in]=1;    in++;
  
    //   Argument ScRatio. Type 2: Double
    strcpy(olarg[in],"-rsc");
    strcpy(olval[in],"ratio");
    sprintf(olexp[in]," %s : %s %g\n",olarg[in],
	    "Ratio between maximum scale and signal length. Default: ",
	    SCRATIO);
    if(siflag)      {
      type[in]=5;
      ptrflag[in]=deflag;
    } else type[in]=2;
    ptrvar[in]= &(p_wav->scratio); 
    ptrval[in][0]=.001;
    ptrval[in][1]=10.;
    olnumb[in]=1;    in++;
  }

  //   Argument THETAU. Type 2: Double
  IF(SelectPar[ipar_theta]) {
    strcpy(olarg[in],"-theta");
    strcpy(olval[in],"angle_1D");
    sprintf(olexp[in]," %s : %s %g\n",olarg[in],
	    "Direction for 1D analysis. Default:",THETA0);
    if(siflag)    {
      type[in]=5;
      ptrflag[in]=deflag;
    }  else type[in]=2;
    ptrvar[in]=&(p_wav->thetau);
    ptrval[in][0]=0;
    ptrval[in][1]=M_PI/2.;
    olnumb[in]=1;  in++;
  }
  
  return in;
}

  


/***************************************************************************/
int filparse_update(  ) {
/***************************************************************************/

  return OK;
} // end of filparse_update


/***************************************************************************/
int filparse_update_fft(  ) {
/***************************************************************************/

//  if(p_fil->flag_nomem == TRUE) p_fil->flag_memory = FALSE;    

  return OK;
} // end of filparse_update_fft


/***************************************************************************/
int filparse_update_wav() {
  /***************************************************************************
   * Put the default parse of WTMM analysis scheme if necessary
   ***************************************************************************/
  
  int flag; // useless temporary flag
  
  if(p_wav->wav == CRAZY) p_wav->wav = WAVWTMM;
  if(p_wav->ordder == CRAZY) p_wav->ordder = ORDDERWTMM;
  
  if(p_wav->sctime == CRAZY) p_wav->sctime = SCTIME;

  if(p_wav->scratio == CRAZY) p_wav->scratio = RATIO;
  
  if(p_wav->minscale == CRAZY) p_wav->minscale = MINSCALE;
  
  if (p_wav->maxscale == CRAZY || p_wav->minscale>=p_wav->maxscale) {
    p_wav->maxscale = p_wav->scratio*(((float)dimx)-1.0)/(float)p_wav->sctime;
    flag=YES;
  } else     flag=NO;
  
  if(p_wav->nvoices == CRAZY)    p_wav->nvoices = NVOICES;
  p_wav->scstep = pow(2.,(1./(double)p_wav->nvoices));
  
  IFVERBOSE {
    WarningVV(" Wavelet to be used in WTMM scheme:   wav_wtmm=%d => %s",
	      p_wav->wav,wav_name[p_wav->wav+1]);
    if(p_wav->wav>0) 
      WarningV(" Derivative order of the gaussian wavelet:"
	       "      deriv_order=%d ",p_wav->ordder);
    WarningV(" Number of voices: no_voices=%d", p_wav->nvoices);
    WarningVV(" Range of scales of analysis in WTMM scheme:"
	      "\n      [scale_0=%g, scale_max=%g]", p_wav->minscale, p_wav->maxscale);
    IF(flag)  
      WarningVV(" Factor of convolution range in WTMM scheme:"
		"      time_scale=%d" 
		"\n Ratio between maximum scale and signal length:"
		"      ratio=%g", p_wav->sctime, p_wav->scratio);
  }

  return OK;
}
