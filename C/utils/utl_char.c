
#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_char.h>


/***************************************************************************/
int get_base( char *name, char splitter, char *base ) {
/***************************************************************************/
  int longo;
  int ic,ic0,ic1;
  
  longo=strlen(name);
  for (ic=longo-1; (ic>=0)&&(name[ic]!=splitter); ic--) ;
  if (ic<=0) {
    ic0=longo-1; // if no splitter found, everything is taken
    ic1=0;
  }  else  {
    ic0=ic;
    /*  We search now for an eventual slashbar splitting the chain  */
    
    for (ic1=ic0; (ic1>0)&&(name[ic1]!='/'); ic1--) ;
    if (name[ic1]=='/') ic1++;
  }
  
  for (ic=ic1; ic<ic0; ic++) base[ic-ic1]=name[ic];
  base[ic-ic1]='\0';
  
  return(ic);
} // end of get_base


/***************************************************************************/
int repeat_char( char*strbar, char c, int l, int close ) {
/***************************************************************************/
 
  int i;
  for( i=0; i<l; i++ ) strbar[i]=c;
  if(close==TRUE) strbar[l] = '\0';

  return OK;
} // end of repeat_char


/***************************************************************************/
int default_name( char *in, char *def, int close ) {
  /***************************************************************************/
  
  int l, k=0;

  if((l=strlen(in)) == 0)     Error("Empty string");

  while(l>=0 && in[l]!='.') l--;
  if(l < 0) Error("No extension in input name (should be of the form '*.*')");
  while(k < l) {
    def[k] = in[k]; 
    k++;
  }
  if(close==TRUE) def[k] = '\0';
 
  return OK;
} // end of default_name


/***************************************************************************/
int define_extension( char *ext, char *c, int dim, int i, int flag, 
		      char *exteff ) {
/***************************************************************************/

  if(flag == TRUE)
    sprintf( exteff, "%s", ext );
  
  else   
    if(dim == 1)   sprintf( exteff, "%s", ext );
    else           sprintf( exteff, "%s-%s%02d", ext, c, i );
  
  return OK;
} // end of define_extension


/***************************************************************************/
int extract_extension( char *nombre, char separador, char *ext) {
  /***************************************************************************/

  int longo;
  int ic,ic0;

  longo=strlen(nombre);
  for(ic=longo-1;(ic>=0)&&(nombre[ic]!=separador);ic--) ;

  ic0=ic+1;
  for(ic=ic0;ic<longo;ic++) ext[ic-ic0]=nombre[ic];
  ext[ic-ic0]='\0';
	
  return(longo-ic0);
} // end of extract_extension
