/* ===================================
** inout.h
** ===================================
*/
#ifndef _INOUT_H_
#define _INOUT_H_

/** Various pseudo-constants **/

#ifndef DSPACE
#define DSPACE DIM2D
#endif

#ifndef DEFZOOM
#define DEFZOOM 1.
#endif

#ifndef BLOCKIN
#define BLOCKIN 1.
#endif
#ifndef BLOCKOUT
#define BLOCKOUT 1.
#endif

/** Various useful type definitions for input/output **/

/* Window to be extracted from a given data file */
#ifndef _DATAWINDOW_TYPE_
#define _DATAWINDOW_TYPE_
typedef struct mydatawindow {
  int ix0;   // X first index
  int dimx;  // X range
  int iy0;   // Y first index
  int dimy;  // Y range
  int iv0;   // V first index
  int dimv;  // V range
  int iz0;   // Z first index
  int dimz;  // Z range
} DataWindow;
#endif

/** Structures for which the Signal structure is required **/
#ifdef _SIGNAL_H_

/* Data reading function pointer prototype */
#ifndef _READSIGNAL_TYPE_
#define _READSIGNAL_TYPE_
typedef int(*ReadSignal)(char*, DataWindow, Signal*);
// char *filename, DataWindow select, Matrix* data
#endif

/* Struct to define the geometry of the processed file */
#ifndef _DATASTRUCTURE_TYPE_
#define _DATASTRUCTURE_TYPE_
typedef struct mydatastructure{
  ReadSignal Reader; // Pointer to the reading routine to be used
  int dimx;
  int dimy;
  int dimv;
  int dimz; // Dimensions of the file
} DataStructure;
#endif
#ifndef _DATASTRUCTURE0
#define _DATASTRUCTURE0 {NULL,0,0,0,0} // DataStructure initializer
#endif

#endif /*!_SIGNAL_H_*/

/* Other useful variables (mainly obsolete) */
#ifndef XMAXVH
#define XMAXVH 2048
#endif
#ifndef YMAXVH
#define YMAXVH 2048
#endif
  //#define XMAXVH 1536
  //#define YMAXVH 1024

#ifndef XMAXINRIA
#define XMAXINRIA 2000
#endif
#ifndef YMAXINRIA
#define YMAXINRIA 2000
#endif

//  Function pointer type for routines in io2D
#ifndef _READ2D_TYPE_
#define _READ2D_TYPE_
typedef int(*Read2D)(int,int,int,int,int,int,char*,double ***);
#endif
#ifndef _READ2D_MASK_TYPE_
#define _READ2D_MASK_TYPE_
typedef int(*Read2D_mask)(int,int,int,int,int,int,char*,double ***, 
			  char **);
#endif

// The prototype of color variable
#ifndef _COLOR_TYPE_
#define _COLOR_TYPE_
typedef struct {
    char Red;
    char Green;
    char Blue;
} color;
#endif

#endif /*!_INOUT_H_*/
