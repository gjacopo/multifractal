#ifndef PIX_ALLOC_H
#define PIX_ALLOC_H

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  Pixel * create_pixel(int m, int n);
  int free_pixel(Pixel *head);
  Pixel * remove_pixel(Pixel *list, int mv, int nv);
  int display_pixel(Pixel *head);
  int create_neighbors (Pixel *temp, int xdim, int ydim);
  
  Region *create_region(int xdim, int ydim); 
  int free_region(Region *reg);

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
