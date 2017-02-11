/* ===================================
** signal.h
** ===================================
*/
#ifndef   	_SIGNAL_H_
#define   	_SIGNAL_H_ 

#ifndef _IMAGE_STRUCT_
#define _IMAGE_STRUCT_
typedef struct myimage {
  char *name;
  char *flag_read;  /* for each frame, tell us if the image has 
		     * been already loaded */
  int bd, type;
  int foto;
  int color;
  int xdim, ydim, vdim;
  int zdim; /* Total number of frames in the sequence the loaded 
	     * image belongs to */
  int iz; /* index of the current loaded image in the full sequence 
	   * if zdim>1 */
  struct myimage *prec; 
  struct myimage *suiv; /* for sequence: previous and successive frames */
  union mysigtype ***image;
} Image;
#endif

#ifndef _IMAGE0_ 
#define _IMAGE0_ {NULL,NULL,0,0,0,0,0,0,0,0,0,NULL,NULL,NULL} 
#endif 	   

#ifndef _SEQUENCE_STRUCT_
#define _SEQUENCE_STRUCT_
typedef struct mysequence {
  Image *head;
  Image *tail;
  int nSeq;
} Sequence;
#endif

#ifndef _SIGNAL_STRUCT_
#define _SIGNAL_STRUCT_
typedef struct mysignal {
  char name[MAXNAMELENGTH];
  int xdim, ydim;
  union mysigtype *mms; /* extrema of the signal (mms[0], mms[1]) 
			 * and range of values (mms[2])*/
  union mysigtype **signal;
  union mysigtype *vsignal;      /* Inner pointer to vectorized data */
  int *flag;
  struct mygradient *grad;
  struct myhisto *hist;
  struct mymask *mask;
} Signal;
#endif

#ifndef _MASK_STRUCT_
#define _MASK_STRUCT_
typedef struct mymask {
  int xdim, ydim;
  char **mask;         // Mask array, in case data are masked
  char *vmask;            // Inner pointer to vectorized mask
} Mask;
#endif


#ifndef _SIGNAL0_ 
#define _SIGNAL0_ {"",0,0,NULL,NULL,NULL,NULL,NULL,NULL} // Signal initalizer
#endif 	   

#endif 	    /* !SIGNAL_H_ */
