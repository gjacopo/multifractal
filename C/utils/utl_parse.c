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
ParseATOM* parser_alloc( int Narg ) {
  /***************************************************************************/ 
  ParseATOM* atom;
  
  TrackNullAlloc( atom=(ParseATOM*)calloc(Narg,sizeof(ParseATOM)) );
  return OK;
} // end of alloc_parser


/***************************************************************************/
int parser_init( ParseATOM* atom ) {
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
int parser_free( ParseATOM* atom ) {
  /***************************************************************************/ 
  Free(atom);
  return OK;
} // end of free_parser


/***************************************************************************/
int parser_load( int Narg, int argc, char *argv[], ParseATOM* atom ) {
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
      else if(!strcmp(argv[lar],atom[in].argname)) {
	flagv = 1; /* break; 
		    * the test (flagv==0) ensures us we leave the loop 'for( in=0;...'
		    * when both strings argv[lar] and atom[in].argname match */
      }
    
    /* Error: none of the known flags match with the entered one:
     * leave the loop 'for( lar=1,flagv=1; ...' */
    if(in>=Narg && flagv==0) break;
    
    in--; /* the right index is in-1 */
    
    if(flagv == 1)
      for ( i=0; i<atom[in].num; i++ ) { 
	/* check all parameters following the flag 
	 * in : there should be olnumb[in] in total */
	
	if (atom[in].type == IS_STRING) { /* the parameter is of type CHAR list */
	  lar++;
	  if((argv[lar]==NULL) || (argv[lar][0]=='-'))
	    ErrorV("Missing argument for option %s",argv[lar-i-1]);
	  strcpy(atom[in+i].cvar,argv[lar]);
	  *(atom[in+i].flag) = TRUE; // string types are assumed to be flagged
	  
	} else if(atom[in].type%3 == IS_INT)  { /* the parameter is of type INT */
	  lar++;
	  sscanf(argv[lar],"%d",atom[in+i].var.i);
	  if((*(atom[in+i].var.i)<atom[in+i].minval.i) ||
	     (*(atom[in+i].var.i)>atom[in+i].maxval.i)) {
	    WarningVVV("Value out of range (%d - %d) for argument %s\n",
		       atom[in+i].minval.i,atom[in+i].maxval.i,atom[in].argname);
	    flagv=0;
	    /* old flagv=-1; */ 
	    break;
	  }
	  if(atom[in].type == 4) *(atom[in+i].flag)=YES;
	  
	}  else if(atom[in].type%3 == IS_DOUBLE)  { /* the parameter is of type DOUBLE */
	  lar++;
	  sscanf(argv[lar],"%lf",atom[in+i].var.d);
	  if((*(atom[in+i].var.d)<atom[in+i].minval.d) ||
	     (*(atom[in+i].var.d)>atom[in+i].minval.d))	{
	    WarningVVV("Value out of range (%f - %f) for argument %s\n",
		       atom[in+i].minval.d,atom[in+i].maxval.d,atom[in].argname);
	    flagv=0;
	    /* old flagv=-1; */ 
	    break;
	  }
	  if(atom[in].type == 5) *(atom[in+i].flag)=YES;
	  
	} else /* if(type[in]%3 == IS_FLAG) : the parameter is a FLAG */
	  /* if(argv[lar][0]=='-') *ptrflag[in+i]=SIGNMINUS*TRUE;
	     else if(argv[lar][0]=='+') *ptrflag[in+i]=SIGNPLUS*TRUE; */
	  *(atom[in+i].flag) = TRUE; 	/* otherwise FALSE */
	  if(atom[in].type == 3) *(atom[in+i].var.i)=1;
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
      if(atom[in].type%3 == 0) {l += fprintf(stderr," [%s]",atom[in].argname);}
      else l += fprintf(stderr," [%s <%s>]",atom[in].argname,atom[in].valname);
      for( i=0,cur=in; i<atom[in].numb[cur]-1; i++ ) 
	in++; // skip some indices 
    }
    Warning("");
  }
  
  if(flagv == 2) {
    Warning("Input parameters :\n==================\n");
    for(in=0;in<Narg;in++) WarningV("%s",atom[in].explain);
  }

  if(flagv != 1) ExitOK;

  return OK;
} // end of init_parser


/***************************************************************************/
int parse_param( int argc, char *argv[], int (*parsers[])(), int nparse ) {
/***************************************************************************/

  ParseATOM* atom;
  int narg;
  //  int nparse=NELEMS(parsers);
  
  atom = parser_alloc( NARGMAX );
  
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

