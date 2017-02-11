#ifndef _PIXREG_H
#define _PIXREG_H

/* structure for a pixel */
typedef struct mypixel {
  int m, n; /* coordinates of the pixel */
  int *iix, *iiy; /* coordinates of the neighbors */
  int nneighb; /* number of neighors */
  int wnbor; /* Size of the window to consider neighbors */
  struct mypixel *next; 
} Pixel;


/* structure for an entire region */
typedef struct myregion {
  int xdim, ydim;
  int num_pixels;
  int num_edge;
  double mean;
  double std;
  double ssquares; /* Sum of squares */
  double ssvalues; /* Sum of values */
  Pixel *pixel_list;
  Pixel *edge_list;
  char **visit;
  char **member;
  struct myregion * next;
} Region;

#ifndef  WNBOR
#define WNBOR 1
#endif

#endif
