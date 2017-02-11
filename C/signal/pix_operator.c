#include <stdio.h>
#include <math.h>

/* Librairies  ROUTINES */

#include <utils.h>
#include <utl_alloc.h>      
#include <utl_stats.h>      

#include <pixreg.h>
#include <pix_alloc.h>
#include <pix_operator.h>


/****************************************************************/
Pixel* file2pixel( FILE *fseed, int *count ) {
  /****************************************************************/
  Pixel *S=NULL, *Shead=NULL, *tmpS;
  int m, n;
  int i, N;
  int len, read;
  char *line = NULL;

  fscanf(fseed,"%d",&N);
  *count = N;
  
  getline(&line, &len, fseed); /* read the end of the line giving the number of 
				* seeds, whatever is after */
  
  while(N>0 && fscanf(fseed,"%d %d",&n, &m) != EOF) {
    getline(&line, &len, fseed); /* read the end of the line giving the
				  * coordinates of the seeds, whatever is after */
    tmpS = create_pixel( m, n );
    N--;
    if (Shead == NULL)	{
      Shead = tmpS;
      S = tmpS;
    } else if (S->next == NULL)	{
      S->next = tmpS;
    } else 	{
      S = S->next;
      S->next = tmpS;
    }
  }
 
  return Shead; 
}


/****************************************************************/
Pixel * mask2pixel( double **m, int xdim, int ydim, double value,
	 int *count ) {
/****************************************************************/

  Pixel *S=NULL, *Shead=NULL, *tmpS;
  int i, j;

  for (i=0; i<xdim; i++)
    for (j=0; j<ydim; j++)
      if(m[j][i] == value) {
	(*count)++;
	tmpS = create_pixel( j, i );
	if (Shead == NULL)	{
	  Shead = tmpS;
	  S = tmpS;
	} else if (S->next == NULL)	{
	  S->next = tmpS;
	} else 	{
	  S = S->next;
	  S->next = tmpS;
	}
      }

  return Shead;  
}

/****************************************************************/
Pixel * cmask2pixel( char **m, int xdim, int ydim, char value,
		     int *count ) {
/****************************************************************/

  Pixel *S=NULL, *Shead=NULL, *tmpS;
  int i, j;
  
  for (i=0; i<xdim; i++)
    for (j=0; j<ydim; j++)
      if(m[j][i] == value) {
	(*count)++;
	tmpS = create_pixel( j, i );
	if (Shead == NULL)	{
	  Shead = tmpS;
	  S = tmpS;
	} else if (S->next == NULL)	{
	  S->next = tmpS;
	} else 	{
	  S = S->next;
	  S->next = tmpS;
	}
      }

  return Shead;  
}


/****************************************************************/
int pixel2cmask( Pixel * pixlist, char value, 
		 int dimx, int dimy, char **mask ) {
  /****************************************************************/
  
  Pixel* tmp = pixlist;

  while (tmp != NULL)    {
    mask[tmp->m][tmp->n] = value;
    tmp = tmp->next;
  }  
  
  return OK;
}

/****************************************************************/
int pixel2mask( Pixel * pixlist, double value, 
		 int dimx, int dimy, double **mask ) {
  /****************************************************************/
  
  Pixel* tmp = pixlist;

  while (tmp != NULL)    {
    mask[tmp->m][tmp->n] = value;
    tmp = tmp->next;
  }  
  
  return OK;
}


/****************************************************************/
int region2mask( Region* R, double value, 
		 int dimx, int dimy, double **mask) {
  /****************************************************************/
 
  while (R != NULL)    {
    pixel2mask( R->pixel_list, value, dimx, dimy, mask );
    R = R->next;
  }
  
  return OK;
}

/****************************************************************/
int region2cmask( Region* R, char value, 
		 int dimx, int dimy, char **mask) {
  /****************************************************************/
 
  while (R != NULL)    {
    pixel2cmask( R->pixel_list, value, dimx, dimy, mask );
    R = R->next;
  }
  
  return OK;
}
