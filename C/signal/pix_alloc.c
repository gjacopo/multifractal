#include <stdio.h>
#include <math.h>

/* Librairies  ROUTINES */

#include <utils.h>
#include <utl_alloc.h>      
#include <utl_stats.h>      

#include <pixreg.h>
#include <pix_alloc.h>

#ifndef  WNBOR
#define WNBOR 1
#endif

/****************************************************************/
Pixel *alloc_pixel( ) {
/****************************************************************/
  Pixel *temp;

  if( (temp=(Pixel*)malloc(sizeof(Pixel))) == NULL)
    return NULL;
  else return temp;
}


/****************************************************************/
Pixel * create_pixel(int m, int n) {
/****************************************************************/
  Pixel *temp;

  if((temp=alloc_pixel()) == NULL) return NULL;

  temp->nneighb = 0;
  temp->wnbor = WNBOR;
  temp->iix = temp->iiy = NULL;

  temp->m = m;
  temp->n = n;
  temp->next = NULL;

  return temp;
}


/****************************************************************/
int create_neighbors (Pixel *temp, int xdim, int ydim) {
/****************************************************************/

  int wnbor=temp->wnbor;

  if( (temp->iix=(int*)malloc((2*wnbor+1)*sizeof(int))) == NULL ||
      (temp->iiy=(int*)malloc((2*wnbor+1)*sizeof(int))) == NULL )
    return ERROR;
  
  vecindex(temp->iix,temp->n,wnbor, xdim, 0 /* no period */);
  vecindex(temp->iiy,temp->m,wnbor, ydim, 0 /* no period */);

 return OK;
}


/****************************************************************/
int free_pixel(Pixel *head){
/****************************************************************/

  Pixel *temp, *temp2;

  temp = head;
  while (temp != NULL) {
    temp2 = temp->next;
    Free(temp->iix);
    Free(temp->iiy);
    Free(temp);
    temp=temp2;
  }

  return OK;
}


/****************************************************************/
int display_pixel(Pixel *head){
/****************************************************************/

  Pixel *temp;

    fprintf(stderr,"\n == display pixel ==");
  temp = head;
  while (temp != NULL) {
    fprintf(stderr,"\ny=%d x=%d",temp->m,temp->n);
    temp = temp->next;
  }

  return OK;
}

/****************************************************************/
Pixel * remove_pixel(Pixel *list, int mv, int nv) {
  /****************************************************************
   * remove_pixel takes a pixel list and a pair of integers (which are 
   * indices).
   * It finds the pixel in the list with these m,n values and removes 
   * it, deallocating its space.  
   ****************************************************************/

  Pixel *temp, *temp2=NULL, *temp3=NULL;
  int count=0;

  temp = list;
  while ((temp != NULL)&&(!((temp->m == mv) && (temp->n == nv)))) {
  /* on recherche l'element avec la valeur val dans la liste */
  /* on conserve aussi dans temp2 l'element qui le */
  /* precende dans la liste.                                 */
    temp2 = temp;
    temp = temp->next;
    count++;
  }

  /* remove the pixel from the list */
  if (temp == NULL) 
    return list; /* l'element n'existe pas dans la liste   */
  
  else {
    if( temp2 == NULL) { /* pas d'element avant, on etait le permier 
			  * on conserve l'element suivant      */
      temp2 = temp->next;
      Free( temp );
      return temp2;
    } else { /* l'element d'avant indique maintenant le suivant          
	      * de l'element trouve (on cours circuite l'element trouve) */
      temp2->next = temp->next;
      Free(temp);
      return list;   
    }
  }
  
}


/****************************************************************/
Region *alloc_region() {
/****************************************************************/
  Region *temp;

  if( (temp=(Region*)malloc(sizeof(Region))) == NULL )
    return NULL;
  else return temp;
}


/****************************************************************/
Region *create_region(int xdim, int ydim)  {
  /****************************************************************/
    
  int i,j;
  Region *temp;

  if( (temp=alloc_region()) == NULL )  return NULL;
  
  temp->mean = temp->std = 0.;
  temp->ssquares = temp->ssvalues = 0.;
  
  temp->num_pixels = temp->num_edge = 0;
  
  temp->pixel_list = NULL;
  temp->edge_list = NULL;

  if( (temp->visit=cmatrix2D(ydim,xdim)) == NULL ||
      (temp->member=cmatrix2D(ydim,xdim)) == NULL )
    return NULL;

  cfill(xdim,ydim,(char)NMARKED,temp->visit,NULL);
  cfill(xdim,ydim,(char)NMARKED,temp->member,NULL);

  temp->ydim = ydim; 
  temp->xdim = xdim;
 
  temp->next = NULL;
  
  return temp;
}


/****************************************************************/
int free_region(Region *reg) {
/****************************************************************/

  Region *temp, *temp2;
  
#ifdef DEBUG
  fprintf(stderr,"Total number of pixels is %d\n", reg->num_pixels);
#endif
  
  temp = reg;
  while(temp != NULL) {
    temp2 = temp->next;
    free_pixel(temp->pixel_list);
    free_pixel(temp->edge_list);
    if(temp->visit != NULL) 
      free_cmatrix2D(temp->visit, temp->ydim);
    if(temp->member != NULL) 
      free_cmatrix2D(temp->member, temp->ydim);
    Free(temp);
    temp = temp2;
    temp = temp->next;
  }
  
  return OK;
}



