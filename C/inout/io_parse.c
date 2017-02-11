#include <stdio.h>
#include <math.h> 
#include <string.h>

#include <utils.h>
#include <utl_alloc.h>
#include <utl_char.h>
#include <utl_parse.h>

#include <inout.h>
#include <io_parse.h>


extern ParIO *p_io;
extern char SelectPar[];

/***************************************************************************/
ParIO* ioparse_alloc() {
  /***************************************************************************/
  ParIO *par;
  
  if( (par = (ParIO*)malloc(sizeof(ParIO))) == NULL ) 
    return NULL;
  else
  return par;
} // end of ioparse_alloc


/***************************************************************************/
int ioparse_free( ParIO *par ) {
  /***************************************************************************/
  Free( par );
  return OK;
} // end of ioparse_free


/***************************************************************************/
ParIO* ioparse_create( ) {
/***************************************************************************/
  ParIO *par;
  
  if( (par=ioparse_alloc()) == NULL )      return NULL;
  ioparse_default(par);

  return par;
}


/***************************************************************************/
int ioparse_default( ParIO *par ) {
/***************************************************************************/
  /* default values */
  par->dim_space = DSPACE;
  par->flag_inr = FALSE;
  par->flag_foto = FALSE;
  par->flag_color = FALSE;
  par->flag_video = FALSE;
  par->flag_visu = FALSE;
  par->flag_window = par->flag_res = FALSE;
  par->block_in = par->block_out = DEFZOOM;

  return OK;
}


/***************************************************************************/
int ioparse_init( int in, int siflag, int *deflag, 
		    char **olarg, char **olval, char **olexp,
		    double **ptrvar, double **ptrval, 
		    int **iptrvar, int **iptrval, char **cptrvar,
		    int **ptrflag, int *type, int*olnumb ) {
  /***************************************************************************/
  
  //    Argument D_space. Type 1: integer
  IF(SelectPar[ipar_dspace]) {
    sprintf(olarg[in],"%s","-d_space");
    sprintf(olval[in],"%s","dimension");
    sprintf(olexp[in]," %s : %s %d\n", olarg[in],
	    "Dimension of the embedding space. Default:",
	    DSPACE); 
    type[in]=1;
    iptrvar[in]=&(p_io->dim_space); 
    iptrval[in][0]=1;
    iptrval[in][1]=2;
    olnumb[in]=1; in++;
  }

  //    Argument INPUT. Type -1: char list
  IF(SelectPar[ipar_in]) {
    sprintf(olarg[in],"%s","-i");
    sprintf(olval[in],"%s","input");
    sprintf(olexp[in]," %s : %s\n", olarg[in],
	    "The name of the input image to be processed."); 
    type[in]=-1;
    cptrvar[in]=&(p_io->in[0]); 
    olnumb[in]=1;  in++;
  }
  
  //    Argument out. Type -1: char list
  IF(SelectPar[ipar_out]) {
    sprintf(olarg[in],"%s","-o");
    sprintf(olval[in],"%s","output");
    sprintf(olexp[in]," %s : %s\n", olarg[in],
	    "The name of the output generic name (def: root of the input name)."); 
    type[in]=-1;
    cptrvar[in]=&(p_io->out[0]); 
    olnumb[in]=1;    in++;
  }
  
  //    Argument out. Type -1: char list
  IF(SelectPar[ipar_ext]) {
    sprintf(olarg[in],"%s","-e");
    sprintf(olval[in],"%s","extension");
    sprintf(olexp[in]," %s : %s\n", olarg[in],
	    "The name of the extension to the output name."); 
    type[in]=-1;
    cptrvar[in]=&(p_io->ext[0]); 
    olnumb[in]=1;    in++;
  }

  //   Argument x0 and y0. Type 1: Integer
  IF(SelectPar[ipar_dim0]) {
    strcpy(olarg[in],"-x0");
    strcpy(olval[in],"x0 y0");
    sprintf(olexp[in]," %s : %s\n",olarg[in],
	    "Extract window: parsing horizontal and vertical coordinates of the left"
	    "\n   up corner.");
    type[in]=1;
    iptrvar[in]=&(p_io->x0); 
    iptrvar[in+1]=&(p_io->y0); 
    iptrval[in][0]=iptrval[in+1][0]=0;
    iptrval[in][1]=iptrval[in+1][1]= INT_MAX;
    olnumb[in]=2;    in+=2;
  }
  
  //   Argument x and y. Type 1: Integer
  IF(SelectPar[ipar_dim]) {
    strcpy(olarg[in],"-x");
    strcpy(olval[in],"x y");
    sprintf(olexp[in]," %s : %s\n",olarg[in],
	    "Extract window: parsing horizontal and vertical sizes.");
    type[in]=1;
    iptrvar[in]=&(p_io->x); 
    iptrvar[in+1]=&(p_io->y); 
    iptrval[in][0]=iptrval[in+1][0]=1;
    iptrval[in][1]=iptrval[in+1][1]= INT_MAX;
    olnumb[in]=2;    in+=2;
  }
  
  //   Argument BIN and BOUT. Type 2: Double
  IF(SelectPar[ipar_res]) {
    strcpy(olarg[in],"-zoom");
    strcpy(olval[in],"block_in block_out");
    sprintf(olexp[in]," %s : %s\n",olarg[in],
	    "Input and output magnifying resolutions (def: no zoom).");
    type[in]=2;
    ptrvar[in]=&(p_io->block_in); 
    ptrvar[in+1]=&(p_io->block_out); 
    ptrval[in][0]=ptrval[in+1][0]=1.;
    ptrval[in][1]=ptrval[in+1][1]=10.;
    olnumb[in]=2;    in+=2;
  }

  //   Argument foto. Type 0: Flag
  IF(SelectPar[ipar_fot]) {
    strcpy(olarg[in],"-pic");
    sprintf(olexp[in]," %s : %s\n",olarg[in],
	    "Flag for picture images (GIF, PPM, PGM).");
    type[in]=0;
    ptrflag[in]=&(p_io->flag_foto); 
    olnumb[in]=1;     in++;
  }
  
  //   Argument INR. Type 0: Flag
  IF(SelectPar[ipar_inr]) {
    strcpy(olarg[in],"-inr");
    sprintf(olexp[in]," %s : %s\n",olarg[in],
	    "Flag for images in inrimage format.");
    type[in]=0;
    ptrflag[in]=&(p_io->flag_inr); 
    olnumb[in]=1;     in++;
  }
  
  //   Argument color. Type 0: Flag
  IF(SelectPar[ipar_col]) {
    strcpy(olarg[in],"-color");
    sprintf(olexp[in]," %s : %s\n",olarg[in],
	    "Flag for color images.");
    type[in]=0;
    ptrflag[in]=&(p_io->flag_color); 
    olnumb[in]=1;   in++;
  }
  
  //   Argument video. Type 0: Flag
  IF(SelectPar[ipar_vid]) {
    strcpy(olarg[in],"-video");
    sprintf(olexp[in]," %s : %s\n",olarg[in],
	    "Flag for video storage. ");
    type[in]=0;
    ptrflag[in]=&(p_io->flag_video); 
    olnumb[in]=1;  in++;
  }
  
  //   Argument visu. Type 0: Flag
  IF(SelectPar[ipar_vis]) {
    strcpy(olarg[in],"-nvis");
    sprintf(olexp[in]," %s : %s\n",olarg[in],
	    "Flag for visualizing normalized original images in GIF format (be" 
	    "\n   carefull that it may increase dramatically the computation time).");
    type[in]=0;
    ptrflag[in]=&(p_io->flag_visu); 
    olnumb[in]=1;   in++;
  }

  return in;
} // end of ioparser_init


/***************************************************************************/
int ioparse_check( ) {
  /***************************************************************************/
  char *str;
  
  IF(SelectPar[ipar_in]) 
    if(!strlen(p_io->in)) {Error("No input given (option '-i')");}
    else str = strrchr( p_io->in, '.' );

  IF(SelectPar[ipar_inr] || SelectPar[ipar_fot]) 
    if(p_io->flag_foto==TRUE && p_io->flag_inr==TRUE)
      Error("Incompatible option '-inr' and '-pic'");
  
  IF(SelectPar[ipar_fot]) 
    if(p_io->flag_foto == FALSE &&
       (!strcmp(str,".gif") || !strcmp(str,".GIF") ||
	!strcmp(str,".pgm") || !strcmp(str,".PGM") ||
	!strcmp(str,".ppm") || !strcmp(str,".PPM")) ) {
      Warning("* Warning: flag '-pic' forced because of input image extension *");
    }
  
  IF(SelectPar[ipar_inr]) 
    if(p_io->flag_inr == FALSE &&
       (!strcmp(str,".inr") || !strcmp(str,".INR")) ) {
      Warning("* Warning: flag '-inr' forced because of input"
	      " image extension *");
    }
  
  IF(SelectPar[ipar_out]) 
    if(strlen(p_io->out) == 0 ) Error("No output given (option '-o')");
  
  return OK;
} // end of ioparse_check


/***************************************************************************/
int ioparse_update( ) {
  /***************************************************************************/
  char *str= strrchr( p_io->in, '.' );
  
  if(!strlen(p_io->ext)) default_name(p_io->in,p_io->ext,FALSE);
    
  if(p_io->flag_foto == FALSE &&
     (!strcmp(str,".gif") || !strcmp(str,".GIF") ||
      !strcmp(str,".pgm") || !strcmp(str,".PGM") ||
      !strcmp(str,".ppm") || !strcmp(str,".PPM")) ) {
    p_io->flag_foto = TRUE;
    p_io->flag_inr = FALSE;
  }
  
  if( p_io->flag_inr == FALSE &&
      (!strcmp(str,".inr") || !strcmp(str,".INR")) ) {
    p_io->flag_inr = TRUE;
    p_io->flag_foto = FALSE;
  } else if(p_io->flag_inr == TRUE) p_io->flag_foto = -TRUE;
  
  return OK;
} // end of ioparse_update



/***************************************************************************/
int ioparse_display( ParIO *par ) {
/***************************************************************************/

  Warning("\n == display IO option ==");
  WarningV("in=%s",par->in);
  WarningV("out=%s",par->out);
  WarningV("ext=%s",par->ext);
  WarningVV("block_in=%f, block_out=%f",par->block_in,par->block_out);
  WarningVV("x0=%d, y0=%d",par->x0,par->y0);
  WarningVV("x=%d, y=%d",par->x,par->y);
  WarningV("pic=%d",par->flag_foto);
  WarningV("color=%d",par->flag_color);
  WarningV("visu=%d",par->flag_visu);
  WarningV("video=%d",par->flag_video);

  return OK;
} // end of ioparse_display


#ifdef Debug
/***************************************************************************/
int ioparse_test (int argc, char *argv[]) {
/***************************************************************************/
  Option* p_io;
  p_io = ioparse_init();
  TrackError( parse_parses(argc,argv,/*FINISH*/), "Error parsing options" );
  ioparse_display( p_io );
  ioparse_free(p_io);
  
  return OK;
} // end of ioparse_test
#endif
