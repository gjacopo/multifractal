/* ===================================
** sig_inout.h
** started on Tue Jan 30 08:52:43 2007 
** ===================================
*/

#ifndef   	_SIG_INOUT_H_
#define   	_SIG_INOUT_H_

/* Data reading function pointer prototype */
typedef int(*ReadSignal)(char* filename, DataWindow select, Signal* data);

/* Struct to define the geometry of the processed file */
typedef struct mydatastructure{
  ReadSignal Reader; // Pointer to the reading routine to be used
  int dimx;
  int dimy;
  int dimv;
  int dimz; // Dimensions of the file
} DataStructure;
#define _DataStructure {NULL,0,0,0,0} // DataStructure initializer


#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */


#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !_SIG_INOUT_H_ */

