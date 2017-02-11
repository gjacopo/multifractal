#ifndef   	_VAR_IO_H_
#define   	_VAR_IO_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

#ifdef _IO_PARSE_H_  

  ParIO p0_io ={
    DSPACE,           // 1: dim_space
    "","","",         // 2: in, out, ext
    FALSE,            // 5: flag_window
    X0,Y0,X,Y,        // 6: x0,y0,x,y
    FALSE,            // 10: flag_res
    DEFZOOM, DEFZOOM, // 11: bin, bout
    FALSE, FALSE,     // 13: flag_foto flag_inr
    FALSE, FALSE,     // 15: flag_color, flag_video
    FALSE,            // 17: flag_visu
    FALSE,            // 18: flag_verbose
    FALSE,            // 19: flag_debug
  }; 

  ParIO *p_io = &p0_io; 

  
#endif

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !VAR_IO_H_ */
