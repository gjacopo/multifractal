#ifndef _CHARUTL_H
#define _CHARUTL_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

  int get_base( char *name, char splitter, char *base );
  int repeat_char(char*strbar, char c, int l, int close );
  int default_name( char *in, char *def, int close );
  int define_extension( char *ext, char *c, int dim, int i, 
			int flag, char *exteff );
  int extract_extension( char *nombre, char separador, char *ext);
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
