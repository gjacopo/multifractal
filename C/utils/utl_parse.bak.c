#include <stdio.h>
#include <math.h> 
#include <string.h>

#include <utils.h>
#include <utl_alloc.h>
#include <utl_char.h>
#include <utl_parse.h>

extern char SelectPar[];
extern int Verbose;


/***************************************************************************/
int parser_alloc( int Narg, char ***olarg, char***olval,char***olexp, 
		  double ***ptrvar, double***ptrval,
		  int ***iptrvar, int***iptrval, int ***ptrflag,
		  char ***cptrvar, int **type, int**olnumb ) {
  /***************************************************************************/ 

  /* Name of the variable */
  TrackNullAlloc( ((*olarg)=cmatrix2D(Narg,MAXCHARLENGTH)) );
  
  /* Explanation and usasge */
  TrackNullAlloc( ((*olexp)=cmatrix2D(Narg,MAXTEXTLENGTH)) ); 
  
  /* Type of variable */
  TrackNullAlloc( ((*type)=(int*)calloc(Narg,sizeof(int))) ); 

  /* Pointer to the DOUBLE variables */
  TrackNullAlloc( ((*ptrvar)=(double**)calloc(Narg,sizeof(double *))) ); 
  TrackNullAlloc( ((*ptrval)=matrix2D(Narg,2)) ); // Allowed min and max
  
  /* Pointer to the INT variables */
  TrackNullAlloc( ((*iptrvar)=(int**)calloc(Narg,sizeof(int *))) ); 
  TrackNullAlloc( ((*iptrval)=imatrix2D(Narg,2)) ); // Allowed min and max
  
  /*  Pointer to the CHAR variables */
  TrackNullAlloc( ((*cptrvar)=(char**)calloc(Narg,sizeof(char*))) );

  /* Pointer to the flags */
  TrackNullAlloc( ((*ptrflag)=(int**)calloc(Narg,sizeof(int *))) ); 
    
  /* Number of arguments to follow each variable */ 
  TrackNullAlloc( ((*olnumb)=(int*)calloc(Narg,sizeof(int))) );
  
  /* Name of the argument associated to each variable */
  TrackNullAlloc( ((*olval)=cmatrix2D(Narg,MAXCHARLENGTH)) );    
  
  return OK;
} // end of alloc_parser


/***************************************************************************/
int parser_init( char **olarg, char **olval, char **olexp,
		 double **ptrvar, double **ptrval, 
		 int **iptrvar, int **iptrval, char **cptrvar,
		 int **ptrflag, int *type, int*olnumb ) {
  /***************************************************************************/
  int in=0;
  
  /* These flags are always displayed */
  
  //   Argument flag VER. Type 0: Flag
  strcpy(olarg[in],"-ver");
  sprintf(olexp[in]," %s : %s%d\n",olarg[in],
	  "Flag for verbose mode. Default: ", FLAG_VERBOSE);
  type[in]=0;
  ptrflag[in] = &Verbose;
  olnumb[in]=1;   in++;
 
  /* These others also don't need pointers */
  //   Argument flag USAGE. Type 0: Flag
  strcpy(olarg[in],"-usage");
  sprintf(olexp[in]," %s : %s\n",olarg[in],
	  "Flag for displaying progam usage.");
  in++;
  
  //   Argument flag HELP. Type 0: Flag
  strcpy(olarg[in],"-help");
  sprintf(olexp[in]," %s : %s\n",olarg[in],
	  "Flag for help.");
  in++;
  
  return in;
} // end of parser_init


/***************************************************************************/
int parser_free( int Narg, char **olarg, char**olval,char**olexp, 
		  double **ptrvar, double**ptrval,
		  int **iptrvar, int**iptrval, int **ptrflag,
		  char **cptrvar, int *type, int*olnumb ) {
  /***************************************************************************/ 

  if(olarg)   free_cmatrix2D(olarg,Narg);
  if(olexp)   free_cmatrix2D(olexp,Narg); 
  if(type)    free(type);
  if(ptrvar)  free(ptrvar);
  if(ptrval)  free_matrix2D(ptrval,Narg); 
  if(iptrvar) free(iptrvar);
  if(iptrval) free_imatrix2D(iptrval,Narg);
  if(cptrvar) free(cptrvar);
  if(ptrflag) free(ptrflag);
  if(olnumb)  free(olnumb);
  if(olval)   free_cmatrix2D(olval,Narg);    
  
  return OK;
} // end of free_parser


/***************************************************************************/
int parser_load( int Narg, int argc, char *argv[],
		 char **olarg, char **olval, char **olexp,
		 double **ptrvar, double **ptrval, 
		 int **iptrvar, int **iptrval,
		 char **cptrvar,
		 int **ptrflag, int *type, int*olnumb ) {
  /***************************************************************************/

  int flagv;
  int lar, arglen;
  int in, i, cur;
  char usage[MAXTEXTLENGTH];
  int l;

  /*                   Parsing loop                         */
  for( lar=1,flagv=1; (lar<argc)&&(flagv==1); lar++ ) {

    flagv = 0;
    for( in=0; (in<Narg)&&(flagv==0); in++ )
      if( !strcmp(argv[lar],"-h") || !strcmp(argv[lar],"-help") ) 
	flagv = 2;
      else if( !strcmp(argv[lar],"-usage") ) 
	flagv = 3;
      else if(!strcmp(argv[lar],olarg[in])) {
	flagv = 1; /* break; 
		    * the test (flagv==0) ensures us we leave the loop 'for( in=0;...'
		    * when both strings argv[lar] and olarg[in] match */
      }

    /* Error: none of the known flags match with the entered one:
     * leave the loop 'for( lar=1,flagv=1; ...' */
    if(in>=Narg && flagv==0) break;

    in--; /* the right index is in-1 */
    
    if(flagv == 1)
      for ( i=0; i<olnumb[in]; i++ ) { 
	/* check all parameters following the flag 
	 * in : there should be olnumb[in] in total */
	
	if (type[in] == -1) { /* the parameter is of type CHAR list */
	  lar++;
	  if((argv[lar]==NULL) || (argv[lar][0]=='-'))
	    ErrorV("Missing argument for option %s",argv[lar-i-1]);
	  strcpy(cptrvar[in+i],argv[lar]);

	} else if(type[in]%3 == 1)  { /* the parameter is of type INT */
	  lar++;
	  sscanf(argv[lar],"%d",iptrvar[in+i]);
	  if((*iptrvar[in+i]<iptrval[in+i][0]) ||
	     (*iptrvar[in+i]>iptrval[in+i][1])) {
	    WarningVVV("Value out of range (%d - %d) for argument %s\n",
		       iptrval[in+i][0],iptrval[in+i][1],olarg[in]);
	    flagv=0;
	    /* old flagv=-1; */ 
	    break;
	  }
	  if(type[in]==4) *ptrflag[in+i]=YES;
	  
	}  else if(type[in]%3 == 2)  { /* the parameter is of type DOUBLE */
	  lar++;
	  sscanf(argv[lar],"%lf",ptrvar[in+i]);
	  if((*ptrvar[in+i]<ptrval[in+i][0]) ||
	     (*ptrvar[in+i]>ptrval[in+i][1]))	{
	    WarningVVV("Value out of range (%f - %f) for argument %s\n",
		       ptrval[in+i][0],ptrval[in+i][1],olarg[in]);
	    flagv=0;
	    /* old flagv=-1; */ 
	    break;
	  }
	  if(type[in]==5) *ptrflag[in+i]=YES;
	  
	} else /* if(type[in]%3 == 0) : the parameter is a FLAG */
	  /* if(argv[lar][0]=='-') *ptrflag[in+i]=SIGNMINUS*TRUE;
	     else if(argv[lar][0]=='+') *ptrflag[in+i]=SIGNPLUS*TRUE; */
	  *ptrflag[in+i] = TRUE; 	/* otherwise FALSE */
	  if(type[in]==3) *iptrvar[in+i]=1;
      }
      
    
    /*            End of parsing loop             */
  }	

  if(flagv == 0)  ErrorV("Unknown option %s",argv[lar]);

  if(flagv == -1) 
    ErrorV("Parameters for option %s out of the range of possible values",
	   argv[lar-i-2]);

  if(/*argc==1 ||*/flagv != 1) {
    WarningV("Usage : %s \n======= \n",argv[0]);
    l=0;
    for( in=0; in<Narg; in++ ) {
      if(l > MAXLINELENGTH) {Warning(""); l=0;} // jump a line
      if(type[in]%3 == 0) {l += fprintf(stderr," [%s]",olarg[in]);}
      else l += fprintf(stderr," [%s <%s>]",olarg[in],olval[in]);
      for( i=0,cur=in; i<olnumb[cur]-1; i++ ) 
	in++; // skip some indices 
    }
    Warning("");
  }
  
  if(flagv == 2) {
    Warning("Input parameters :\n==================\n");
    for(in=0;in<Narg;in++) WarningV("%s",olexp[in]);
  }

  if(flagv != 1) ExitOK;

  return OK;
} // end of init_parser


/***************************************************************************/
int parse_param( int argc, char *argv[], int (*parsers[])(), int nparse ) {
/***************************************************************************/

  char **olarg,**olval,**olexp;
  double **ptrvar,**ptrval;
  int **iptrvar,**iptrval;
  int **ptrflag;
  int *type, *olnumb;
  int narg;
  char **cptrvar;
  //  int nparse=NELEMS(parsers);
  
  parser_alloc( NARGMAX, &olarg, &olval, &olexp, &ptrvar, &ptrval,
		&iptrvar, &iptrval, &ptrflag, &cptrvar, &type, &olnumb );
  
  narg = parser_init( olarg, olval, olexp,
		      ptrvar, ptrval, iptrvar, iptrval,
		      cptrvar,ptrflag, type, olnumb );
  
  while(nparse-- >= 0)
    narg += (parsers[nparse])( narg, olarg, olval, olexp,
			       ptrvar, ptrval, iptrvar, iptrval,
			       cptrvar,ptrflag, type, olnumb );
  
  parser_load( narg, argc, argv, olarg, olval, olexp,
	       ptrvar, ptrval, iptrvar, iptrval,
	       cptrvar,ptrflag, type, olnumb );
  
  parser_free( narg, olarg, olval,olexp, ptrvar, ptrval,
	       iptrvar, iptrval, ptrflag, cptrvar, type, olnumb );
  
  return OK;
} // end of parse_param


/***************************************************************************/
int parser_select( const int choice[], int nopt ) {
  /***************************************************************************/
  int i;
  for( i=0; i<nopt; i++ ) SelectPar[i] = choice[i];
  for( i=nopt; i<NARGMAX;i++ ) SelectPar[i] = FALSE;
  return OK;
} // end of parser_select

