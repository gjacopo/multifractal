#ifndef   	_UTL_PARSE_H_
#define   	_UTL_PARSE_H_

#ifndef _PARSE_BASE_PARAMETERS_
#define _PARSE_BASE_PARAMETERS_
#endif

#ifndef _UTL_UTYPE_UNION_
#define _UTL_UTYPE_UNION_
typedef union utl_utype{
  float f;
  int i;
  double d;
  char c;
} UType;
#endif

enum {IS_STRING=-1, IS_FLAG, IS_INT, IS_FLOAT, IS_DOUBLE, IS_FLAGGED};
/* type of variable: 
 * 0: flag, 1: int; 2: float, +3 if flagged, -1: string
 */

/* Definition of the parse atom (for parse_argument routines) */
#ifndef _UTL_ParseATOM_STRUCT_
#define _UTL_ParseATOM_STRUCT_
typedef struct utl_parseatom{
  char  argname[MAXCHARLENGTH];    // name of the argument
  char  valname[MAXTEXTLENGTHMAX_NAME];    // name of its value
  char  explain[MAXTEXTLENGTH]; // explanation (as used by help)
  int   type;               // type of variable
  int   *flag;              // pointer to variable, if flag
  int   num;              // number of expected parameters
  UType *var;               // pointer to variable
  UType  minval;            // minimum value
  UType  maxval;            // maximum value
  char  *cvar;              // pointer to variable if string
} ParseATOM;

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  int parser_alloc( int Narg, char ***olarg, char***olval,char***olexp, 
		    double ***ptrvar, double***ptrval,
		    int ***iptrvar, int***iptrval, int ***ptrflag,
		    char ***cptrvar, int **type, int**olnumb );
  int parser_init( char **olarg, char **olval, char **olexp,
		   double **ptrvar, double **ptrval, 
		   int **iptrvar, int **iptrval, char **cptrvar,
		   int **ptrflag, int *type, int*olnumb );
  int parser_free( int Narg, char **olarg, char**olval,char**olexp, 
		   double **ptrvar, double**ptrval,
		   int **iptrvar, int**iptrval, int **ptrflag,
		   char **cptrvar, int *type, int*olnumb );
  int parser_load( int Narg, int argc, char *argv[],
		   char **olarg, char **olval, char **olexp,
		   double **ptrvar, double **ptrval, 
		   int **iptrvar, int **iptrval,
		   char **cptrvar,
		   int **ptrflag, int *type, int*olnumb );
  int parse_param( int argc, char *argv[], int (*parsers[])(), int nparse );
  int parser_select( const int choice[], int nopt );

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !_UTL_PARSE_H_ */
