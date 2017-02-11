
#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_alloc.h>	
#include <utl_stats.h>

#include <flt_stats1d.h>



/***************************************************************************/
double cuantil1D( int dimx, double *data, double prob0) {
  /***************************************************************************/
  double *sort;
  double quant;
  int ix;

  sort=(double *) calloc(dimx,sizeof(double));

  for(ix=0;ix<dimx;ix++) sort[ix]=data[ix];
  quicksort1D(0,dimx-1,sort);

  ix=(int)((1.-prob0)*(double)(dimx-1));
  quant=sort[ix];

  free(sort);
  return(quant);
}


/***************************************************************************/
int quicksort1D( int low, int high, double *data ) {
  /***************************************************************************/
  int pivot;

  if(low<high) {   
    pivot=partition1D(low,high,data);
    quicksort1D(low,pivot-1,data);
    quicksort1D(pivot+1,high,data);
  }

  return OK;
}

/***************************************************************************/
int quicksort1D_ref( int low, int high, int *ref, double *data ) {
  /***************************************************************************/
  
  int pivot;
  
  if(low<high)    {   
    pivot=partition1D_ref(low,high,ref,data);
    quicksort1D_ref(low,pivot-1,ref,data);
    quicksort1D_ref(pivot+1,high,ref,data);
  }
 
  return OK;
}

/***************************************************************************/
int partition1D( int low, int high, double *data) {
  /***************************************************************************/
  double pivot_item;
  double buff;
  int left,right;

  pivot_item=data[low];
  left=low;
  right=high;
  while(left<right) {
    while((data[left]<=pivot_item)&&(left<=right)) left++;
    while((data[right]>=pivot_item)&&(left<=right)) right--;
    if(left<right)  {
      buff=data[left];
      data[left]=data[right];
      data[right]=buff;
    }
  }
  data[low]=data[right];
  data[right]=pivot_item;

  return(right);

}

/***************************************************************************/
int partition1D_ref( int low, int high, int *ref, double *data) {
  /***************************************************************************/
  double pivot_item;
  double buff;
  int pivot_item_i,buffi;
  int left,right;

  pivot_item=data[low];
  pivot_item_i=ref[low];
  left=low;
  right=high;

  while(left<right)    {
    while((data[left]<=pivot_item)&&(left<=right)) left++;
    while((data[right]>=pivot_item)&&(left<=right)) right--;
    if(left<right)      {
      buff=data[left];
      data[left]=data[right];
      data[right]=buff;
      buffi=ref[left];
      ref[left]=ref[right];
      ref[right]=buffi;
    }   
    }
 
  data[low]=data[right];
  ref[low]=ref[right];
  data[right]=pivot_item;
  ref[right]=pivot_item_i;

  return right;
}


/***************************************************************************/
double moda1D(int dimx, double *datos)/*moda_1D*/ {
/***************************************************************************/
  double *histo;
  double mm[2];
  double out;
  int Nhisto=dimx/30; // 30 events by bin in average
  int ix,ip,ip0;

  TrackNullAlloc( histo=(double*)calloc(Nhisto,sizeof(double)) );

  extrema1D(dimx,datos,&mm[0],NULL);
  for(ix=0;ix<dimx;ix++)    {
    ip=(int)(((double)Nhisto)*(datos[ix]-mm[0])/(mm[1]-mm[0]));
    if(ip>Nhisto-1) ip=Nhisto-1;
    histo[ip]+=1.;
  }

  out = moda1D_histo(Nhisto,mm,histo);
  /* for(ip0=0,ip=0;ip<Nhisto;ip++) if(histo[ip]>histo[ip0]) ip0=ip;
     out=mm[0]+(mm[1]-mm[0])*(0.5+(double)ip0)/((double)Nhisto);
  */

  free(histo);

  return out;
}


/****************************************************************************/
double moda1D_histo( int Nhisto, double *mm, double *histo ) /*moda_por_histo*/{
  /****************************************************************************/
  double out;
  int ip,ip0;
    
  for(ip0=0,ip=0;ip<Nhisto;ip++) 
    if(histo[ip]>histo[ip0]) ip0=ip;
  out = mm[0] + (mm[1]-mm[0]) *(0.5+(double)ip0) / ((double)Nhisto);
  
  return out;
}



/****************************************************************************/
int extrema1D( int dimx, double *cont, double *mm, char *m )/*extrema_lista*/ {
  /****************************************************************************/
  /* Computes the extrema of the image */
  /****************************************************************************/
  
  int ix;
  int first;
  mm[0] = 1.e+17, mm[1] = -1.e17;

  if(m == NULL) {
    mm[0]= mm[1]= cont[0];
    for( ix=0; ix<dimx; ix++ )      {
      mm[0] = fMin(mm[0],cont[ix]);
      mm[1] = fMax(mm[1],cont[ix]);
    }
    
  } else {
    for( ix=0; ix<dimx; ix++ )      
      if(m[ix] == (char)TRUE) {
	/*if(first == TRUE) {
	  mm[0] = mm[1] = cont[ix]; 
	  first = FALSE;
	  }	*/  
	mm[0] = fMin(mm[0],cont[ix]);
	mm[1] = fMax(mm[1],cont[ix]);
      }
  }
  
  return OK;
}


/****************************************************************************/
double media1D( int dimx, double *cont, char *m )/*media_lista*/ {
  /****************************************************************************/
  /* Computes the mean of the vector  */
  /****************************************************************************/  
  int ix, count=0;
  double suma=0.;
  
  if(m == NULL) {
    for(ix=0;ix<dimx;ix++)	suma += cont[ix];
    count = dimx;

  }  else
    for(ix=0;ix<dimx;ix++)
	if(m[ix] == (char)TRUE) {
	  suma += cont[ix];
	  count ++;
	}
  
  suma /= ((double) count);
  
  return suma;
}


/****************************************************************************/
double dispersion1D( int dim, double *cont, char *m ) /*dispersion_lista*/{
  /****************************************************************************/
  /* Computes the dispersion (standard deviation) of a signal  */
  /****************************************************************************/

  int i, count=0;
  double sumaxx=0.,sumax=0.;

  if(m == NULL) {
    for( i=0; i<dim; i++ ) {
      sumax += cont[i];
      sumaxx += cont[i] * cont[i];
    }
    count = dim; 
    
  } else
    for(i=0;i<dim;i++)  if(m[i] == (char)TRUE) {
      sumax += cont[i];
      sumaxx += cont[i] * cont[i];
      count ++;
    }
  
  sumax /= ((double) count);
  sumaxx = sumaxx/((double) count)-sumax*sumax;
  sumaxx = sqrt(fabs(sumaxx));

  return sumaxx;
}


/****************************************************************************/
double covariance1D( int dimx, double *x, double *y, char *m)/*covarianza_lista*/ {
/****************************************************************************/
  double sumaxy=0.,sumax=0.,sumay=0.;
  int ix, count=0;

  if(m == NULL) {
    for(ix=0;ix<dimx;ix++) {
      sumax+=x[ix];
      sumay+=y[ix];
      sumaxy+=x[ix]*y[ix];
    }
    count = dimx;

  } else
    for(ix=0;ix<dimx;ix++)  if(m[ix] == (char)TRUE) {
      sumax+=x[ix];
      sumay+=y[ix];
      sumaxy+=x[ix]*y[ix];
      count ++;
    }
  
  sumax = sumax/((double) count);
  sumay = sumay/((double) count);
  sumaxy = sumaxy/((double) count)-sumax*sumay;

  return(sumaxy);
}


/***************************************************************************/
double anorma1D( int dimx, int dimy, double *cont, char*m )/*anorma_lista*/ {
  /***************************************************************************/
  /* Substracts the mean of the vector to itself  */
  /***************************************************************************/

  double media=0.;
  int ix, count=0;
  
  if(m == NULL) {
      for(ix=0;ix<dimx;ix++)	media+=cont[ix];
      count = dimx*dimy;

  }   else
      for(ix=0;ix<dimx;ix++)
	if(m[ix] == (char)TRUE) {
	  media+=cont[ix];
	  count ++;
	}
  
  media /= ((double)count);
  
  if(m == NULL) 
      for(ix=0;ix<dimx;ix++)
	cont[ix] -= media;
  
  else
      for(ix=0;ix<dimx;ix++)
	if(m[ix] == (char)TRUE) cont[ix] -= media;
  
  return media;
}


/***************************************************************************/
double anorma11D( int dimx, double *cont, char*m )/*anorma1_lista*/ {
  /***************************************************************************/
  /* Divides the image by its absolute mean  */
  /***************************************************************************/

  double media=0.;
  int ix, count=0;
  
  if(m == NULL) {
    for(ix=0;ix<dimx;ix++)	media+=fabs(cont[ix]);
    count = dimx;
  }   else
    for(ix=0;ix<dimx;ix++)
      if(m[ix] == (char)TRUE) {
	media+=fabs(cont[ix]);
	count ++;
      }
  
  media /= ((double)count);
  
  if(m == NULL) 
    for(ix=0;ix<dimx;ix++)   cont[ix] /= media;
  
  else
    for(ix=0;ix<dimx;ix++)
      if(m[ix] == (char)TRUE) cont[ix] /= media;
  
  return media;
}

/***************************************************************************/
double anorma21D( int dimx, double *cont, char*m ) {
  /***************************************************************************/
  /* Divides the image by its L2 norma */
  /***************************************************************************/

  double media=0.;
  int ix, count=0;
  
  if(m == NULL) {
    for(ix=0;ix<dimx;ix++)
      media+=cont[ix]*cont[ix];
    count = dimx;
  }   else
    for(ix=0;ix<dimx;ix++)
      if(m[ix] == (char)TRUE) {
	media+=cont[ix]*cont[ix];
	count ++;
      }
  
  media = sqrt(media / ((double)count));
  
  if(m == NULL) 
    for(ix=0;ix<dimx;ix++) cont[ix] /= media;
  else
    for(ix=0;ix<dimx;ix++) if(m[ix] == (char)TRUE) cont[ix] /= media;
  
  return media;
}

/***************************************************************************/
int denorma1D( int dimx, double norma, double *cont, char*m ) /*denorma_lista*/ {
  /***************************************************************************/
  /* Adds the value norma (generally, the mean of the image) to the image */
  /***************************************************************************/
  int ix;
  
  if(m == NULL)
    for(ix=0;ix<dimx;ix++)
      cont[ix] += norma;
  
  else
    for(ix=0;ix<dimx;ix++)
      if(m[ix] == (char)TRUE) cont[ix] += norma;
  
  return OK;
}

