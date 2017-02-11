#ifndef   	UTILS_H_
#define   	UTILS_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  /* Size and format related constants */

#ifndef LITTLE_INDIAN
#define LITTLE_INDIAN 1
#endif
  
#ifndef I_TBYTE 
#define I_TBYTE 0
typedef char TBYTE;
#endif
#ifndef I_TCHAR 
#define I_TCHAR 1
typedef char TCHAR;
#endif
#ifndef I_TUCHAR
#define I_TUCHAR 2
typedef unsigned char TUCHAR;
#endif
#ifndef I_TSHORT
#define I_TSHORT 3
typedef int TSHORT;
#endif
#ifndef I_TUSHORT
#define I_TUSHORT 4
typedef unsigned short TUSHORT;
#endif
#ifndef I_TINT 
#define I_TINT 5
typedef int TINT;
#endif
#ifndef I_TUINT
#define I_TUINT 6
typedef unsigned int TUINT;
#endif
#ifndef I_TLONG
#define I_TLONG 7
typedef long TLONG;
#endif
#ifndef I_TULONG
#define I_TULONG 8
typedef unsigned long TULONG;
#endif
#ifndef I_TFLOAT
#define I_TFLOAT 9
typedef float TFLOAT;
#endif
#ifndef I_TDOUBLE
#define I_TDOUBLE 10
typedef double TDOUBLE;
#endif

#ifndef SZ_TCHAR
#define SZ_TCHAR   sizeof(TCHAR)
#endif
#ifndef SZ_TUCHAR
#define SZ_TUCHAR  SZ_CHAR
#endif
#ifndef SZ_TBYTE
#define SZ_TBYTE   SZ_CHAR
#endif
#ifndef SZ_TSHORT
#define SZ_TSHORT  sizeof(TSHORT)
#endif
#ifndef SZ_TUSHORT
#define SZ_TUSHORT  SZ_TSHORT
#endif
#ifndef SZ_TINT
#define SZ_TINT    sizeof(TINT)
#endif
#ifndef SZ_TUINT
#define SZ_TUINT    SZ_TINT
#endif
#ifndef SZ_TLONG
#define SZ_TLONG   sizeof(TLONG)
#endif
#ifndef SZ_TULONG
#define SZ_TULONG   SZ_TLONG
#endif
#ifndef SZ_TFLOAT
#define SZ_TFLOAT  sizeof(TFLOAT)
#endif
#ifndef SZ_TDOUBLE
#define SZ_TDOUBLE sizeof(TDOUBLE)
#endif

#ifndef _ALLTYPES_UNION_
#define _ALLTYPES_UNION_
typedef union typeun{
  TCHAR c;
  TUCHAR uc;
  TSHORT s;
  TUSHORT us;  
  TINT i;
  TUINT ui;
  TLONG l;
  TULONG ul;
  TFLOAT f;
  TDOUBLE d;
} TypeUn;
#endif

#ifndef UNSUPPORTEDTYPE
#define UNSUPPORTEDTYPE Exit("Unsupported data type");
#endif
#ifndef CASETYPE
#define CASETYPE(itype,do_char,do_uchar,do_short,do_ushort, \
          do_int,do_uint,do_long,do_ulong,do_float,do_double) \
switch (itype) { \
 case I_TCHAR:  (do_char); break; \
 case I_TUCHAR: (do_uchar); break; \
 case I_TSHORT: (do_short); break; \
 case I_TUSHORT: (do_ushort); break; \
 case I_TINT: (do_int); break; \
 case I_TUINT: (do_uint); break; \
 case I_TLONG: (do_long); break; \
 case I_TULONG: (do_ulong); break; \
 case I_TFLOAT: (do_float); break; \
 case I_TDOUBLE: (do_double); break; \
 default:  UNSUPPORTEDTYPE; \
}
#endif

#ifndef SETUNIONVAL
#define SETUNIONVAL(v,itype,p) \
CASETYPE( itype, \
   (p).c = v, \
   (p).uc = v, \
   (p).s = v, \
   (p).us = v, \
   (p).i = v, \
   (p).ui = v, \
   (p).l = v, \
   (p).ul = v, \
   (p).f = v, \
   (p).d = v \
)
#endif

#ifndef GETUNIONVAL
#define GETUNIONVAL(p,itype,v) \
CASETYPE( itype, \
   v = (p).c, \
   v = (p).uc, \
   v = (p).s, \
   v = (p).us, \
   v = (p).i, \
   v = (p).ui, \
   v = (p).l, \
   v = (p).ul, \
   v = (p).f, \
   v = (p).d \
)
#endif



#ifndef MIN_NUM_RESOLUTION
#define MIN_NUM_RESOLUTION 1e-17 // Maximum allowed numerical resolution
#endif

  /* Useful constants for parametrisation */
#ifndef NARGMAX
#define NARGMAX 100 /* Maximum number of arguments that can be passed */
#endif
  
#ifndef NARG
#define NARG NARGMAX /* see option files */
#endif

#ifndef IO_NOPT
#define IO_NOPT 0
#endif
#ifndef OP_NOPT
#define OP_NOPT 0
#endif
#ifndef FEA_NOPT
#define FEA_NOPT 0
#endif

  /* */

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#ifndef ERROR
#define ERROR -1
#endif
#ifndef OK
#define OK 0
#endif

#ifndef YES
#define YES 1
#endif
#ifndef NO
#define NO 0
#endif

#ifndef MARKED
#define MARKED 1
#endif
#ifndef NMARKED
#define NMARKED (1-MARKED)
#endif

#ifndef MAXCHARLENGTH
#define MAXCHARLENGTH  90
#endif
#ifndef MAXTEXTLENGTH
#define MAXTEXTLENGTH  500
#endif
#ifndef MAXNAMELENGTH
#define MAXNAMELENGTH MAXCHARLENGTH
#endif
#ifndef MAXLINELENGTH
#define MAXLINELENGTH  70
#endif
#ifndef MAXCHARBUFFLENGTH
#define MAXCHARBUFFLENGTH  255
#endif

#ifndef M_PI
  //#define M_PI   3.1415926535997932
#define M_PI 3.1415926535897932384626433832795028841971 
#endif

#ifndef MAX_DOUBLE
#define MAX_DOUBLE +1.e17
#endif

#ifndef INT_MAX
#define INT_MAX 2147483647
#endif
#ifndef MAX_INT
#define MAX_INT INT_MAX
#endif

#ifndef WHITE
#define WHITE 1
#endif
#ifndef BLACK
#define BLACK 0
#endif

#ifndef SIGNMINUS
#define SIGNMINUS -1 
#endif
#ifndef SIGNZERO
#define SIGNZERO 0
#endif
#ifndef SIGNPLUS
#define SIGNPLUS 1
#endif

#ifndef TOL
#define TOL 1.e17
#endif
#ifndef PCTOL
#define PCTOL 20.
#endif
#ifndef PCTMAX
#define PCTMAX 60.
#endif

#ifndef DIM2D
#define DIM2D 2
#endif
#ifndef DIM1D
#define DIM1D 1
#endif


/** Some useful variables/routines **/

#ifndef SC2SIZE
#define SC2SIZE(s)      (ROUND(3.*s)) 
#endif

#ifndef SWAP
#define SWAP(i,j)       { int _T=i; i=j; j=_T;}
#endif
#ifndef ABS
#define ABS(x)	        (((x)<0) ? -(x) : (x))
#endif

#ifndef BASELOG
#define BASELOG ((double)10.)
#endif

#ifndef LOG
#define LOG(x) (log((x))/log(BASELOG))
#endif

#ifndef LOG2
#define LOG2(x)          log(x) / log(2.)
#endif

#ifndef LOG10
#define LOG10(x)          log(x) / log(10.)
#endif

#ifndef MIN
#define MIN(a,b)        ((a)<(b)) ? (a) : (b)   
#endif
#ifndef MAX
#define MAX(a,b)        ((a)>(b)) ? (a) : (b)   
#endif

#ifndef SIGN
#define SIGN(x) (((x)>0.)?(1):(-1))
#endif

#ifndef ROUND
#define ROUND(f)        ((f>0) ? (int)(f+0.5) : (int)(f-0.5))
#endif

#ifndef SQR
#define SQR(x) ( (x) * (x) )
#endif
#ifndef SQRT
#define SQRT(x) ( sqrt(x) )
#endif


#ifndef NELEMS
#define NELEMS(t) (sizeof(t) / sizeof *(t))
#endif

#ifndef Fill1D
#define Fill1D(I,x,V)   { int X; \
                            for (X=0; X<(x); X++) \
                                    I[X] = (V);}
#endif

#ifndef Fill2D
#define Fill2D(I,x,y,V) { int X, Y; \
                            for (X=0; X<(x); X++) \
                              for (Y=0; Y<(y); Y++) \
                                    I[X][Y] = (V);}
#endif

#ifndef Free
#define Free(x)   {if((x)!=NULL) {free(x); x=NULL;}}
#endif

  /* Error handling */
  
#ifndef MAXMSGERROR
#define MAXMSGERROR 200
#endif

#ifndef Warning
#define Warning(f) {fprintf(stderr,"\n %s",(f));}
#endif
  
#ifndef WarningV
#define WarningV(f,g) {		       	\
    char warnmsg[MAXMSGERROR];	       	\
    sprintf(warnmsg,(f),(g));	      	\
    Warning(warnmsg); }
#endif
  
#ifndef WarningVV
#define WarningVV(f,g,h) {	      	\
    char warnmsg[MAXMSGERROR];	      	\
    sprintf(warnmsg,(f),(g),(h));	\
    Warning(warnmsg); }
#endif
  
#ifndef WarningVVV
#define WarningVVV(f,g,h,i) {	       	\
    char warnmsg[MAXMSGERROR];	       	\
    sprintf(warnmsg,(f),(g),(h),(i));  	\
    Warning(warnmsg); }
#endif


#ifndef WarnMsg
#define WarnMsg(f) {			 \
    char warnmsg[MAXMSGERROR];	       	 \
    sprintf(warnmsg,"(in %s): %s",__FUNCTION__,(f));	\
    Warning(warnmsg); }
#endif

#ifndef WarnMsgV
#define WarnMsgV(f,g) {	    \
    char warnmsgv[MAXMSGERROR];		    \
    sprintf(warnmsgv,"(in %s): ",__FUNCTION__); \
    strcat(warnmsgv,(f));			\
    WarningV(warnmsgv,(g)); }
#endif

#ifndef WarnMsgVV
#define WarnMsgVV(f,g,h) {		    \
    char warnmsgv[MAXMSGERROR];		    \
    sprintf(warnmsgv,"(in %s): ",__FUNCTION__); \
    strcat(warnmsgv,(f));			\
    WarningVV(warnmsgv,(g),(h)); }
#endif

#ifndef Error
#define Error(f) { \
WarnMsg(f); \
return ERROR; }
#endif

#ifndef ErrorV
#define ErrorV(f,g) { \
WarnMsgV(f,g); \
return ERROR; }
#endif

#ifndef ExitOK
#define ExitOK { \
Warning("Exiting the program...\n"); \
exit (OK); }
#endif

#ifndef Exit
#define Exit(f) { \
WarnMsg(f); \
Warning("Exiting the program...\n"); \
exit (ERROR); \
}
#endif

#ifndef TrackError
#define    TrackError(f,g) {if((f) == ERROR) Error(g)}
#endif

#ifndef TrackErrorAlloc
#define TrackErrorAlloc(f) {TrackError(f,"Error pointer allocation");}
#endif

#ifndef TrackNull
#define     TrackNull(f,g) {if((f) == NULL) Error(g)}
#endif

#ifndef TrackNullAlloc
#define  TrackNullAlloc(f) {TrackNull(f,"Error pointer allocation");}
#endif

#ifndef TrackEOF
#define      TrackEOF(f,g) {if((f) == EOF) Error(g)}
#endif

#ifndef TrackZero
#define     TrackZero(f,g) {if((f) == 0) Error(g)}

#ifndef IsError
#define       IsError(f,g) {if((f)== ERROR) Exit(g)}
#endif

#ifndef IsNull
#define        IsNull(f,g) {if((f)== NULL) Exit(g)}
#endif

#ifndef IsNullAlloc
#define     IsNullAlloc(f) {{IsNull(f,"Error pointer allocation");}}
#endif

#ifndef IF
#define IF(f) if(f == TRUE)
#endif
#ifndef ELSE
#define ELSE else
#endif

#ifndef FLAG_VERBOSE 
#define FLAG_VERBOSE FALSE
#endif
#ifndef IFVERBOSE 
#define IFVERBOSE IF(flag_verbose)
#endif

#endif

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !UTILS_H_ */
