
#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_alloc.h>	
#include <utl_stats.h>

#include <flt_stats1d.h>
#include <flt_stats2d.h>


/***************************************************************************/
double cuantil2D( int dimx, int dimy, double **data, double prob0) {
  /***************************************************************************/
  
  double quant;

  if(dimx*dimy<10000) 
    quant = cuantil2D_sort(dimx,dimy,data,prob0);
  else 
    quant = cuantil2D_histo(dimx,dimy,data,prob0);

  return quant;
} // end of cuantil2D


/***************************************************************************/
double cuantil2D_sort( int dimx, int dimy, double **data, double prob0) {
  /***************************************************************************/
  
  double *sort;
  double quant;
  int ix,iy;
  
  sort=(double*)calloc(dimx*dimy,sizeof(double));
  
  for(iy=0;iy<dimy;iy++)    
    for(ix=0;ix<dimx;ix++)     
      sort[ix+dimx*iy]=data[iy][ix];
  
  quicksort1D(0,dimx*dimy-1,sort);
  
  ix=(int)((1.-prob0)*(double)(dimx*dimy-1));
  quant=sort[ix];
  
  free(sort);
  return(quant);
} // end of cuantil2D_sort


/***************************************************************************/
double cuantil2D_histo( int dimx, int dimy, double **data, double prob0) {
  /***************************************************************************/

  double *histo;
  double maxd,mind,norma;
  double quant;
  const int Nhisto=10000;
  int ix,iy,ip;
  
  histo = (double*)calloc(Nhisto+1,sizeof(double));
  norma = 1./((double)(dimx*dimy));
  
  mind=maxd=data[0][0];
  for(iy=0;iy<dimy;iy++)  
    for(ix=0;ix<dimx;ix++)	{
      mind=fMin(mind,data[iy][ix]);
      maxd=fMax(maxd,data[iy][ix]);
    }
  maxd-=mind;
  
  if(maxd>1e-30)	{
    for(iy=0;iy<dimy;iy++)	 
      for(ix=0;ix<dimx;ix++)	  {
	ip=(int)(Nhisto*(data[iy][ix]-mind)/maxd);
	histo[ip]+=norma;
      }
    norma=0.;
    for(ip=0;(ip<=Nhisto)&&(norma<(1.-prob0));ip++) norma+=histo[ip];
    quant=mind+maxd*((double)ip)/((double)Nhisto);
  }
  else quant=mind;
  
  free(histo);
  return(quant);
} // end of cuantil2D_histo


/***************************************************************************/
double moda(int dimx, int dimy, double **datos)/*moda_2D*/ {
  /***************************************************************************/
  double *histo;
  double mm[2];
  double out;
  int Nhisto;
  int ix,iy,ip,ip0;

  Nhisto=(dimx*dimy)/30; // 30 events by bin in average
  TrackNullAlloc( histo=(double*)calloc(Nhisto,sizeof(double)) );

  extrema(dimx,dimy,datos,&mm[0],NULL);

  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)	{
      ip=(int)(((double)Nhisto)*(datos[iy][ix]-mm[0])/(mm[1]-mm[0]));
      if(ip>Nhisto-1) ip=Nhisto-1;
      histo[ip]+=1.;
    }

  out = moda1D_histo(Nhisto,mm,histo);
  /* for(ip0=0,ip=0;ip<Nhisto;ip++) if(histo[ip]>histo[ip0]) ip0=ip;
     out=mm[0]+(mm[1]-mm[0])*(0.5+(double)ip0)/((double)Nhisto);
  */
  
  free(histo);

  return out;
} // end of moda


/***************************************************************************/
double moda_redef(int leff, double **signal ) {
  /***************************************************************************/

#ifdef _PARSE_FRACTAL_PARAMETERS_
  if(p_frac->dim_space == DIM1D) 
#else
    if(DSPACE == DIM1D) 
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */
      return moda1D(leff,signal[0]);
    else 
      return moda(leff,leff,signal);
} // end of moda_redef


/****************************************************************************/
int cextrema( int dimx, int dimy, char **cont, int *mm, char **m ) {
  /****************************************************************************/
  /* Computes the extrema of the image */
  /****************************************************************************/
  
  int ix,iy;
  int first=TRUE;
  
  
  if(m == NULL) {
    mm[0]= mm[1]= (int)cont[0][0];
    for( iy=0; iy<dimy; iy++ )
      for( ix=0; ix<dimx; ix++ )       {
	mm[0] = Min(mm[0],(int)cont[iy][ix]);
	mm[1] = Max(mm[1],(int)cont[iy][ix]);
      } 
    
  } else
    for( iy=0; iy<dimy; iy++ )
      for( ix=0; ix<dimx; ix++ )       
	if(m[iy][ix] == (char)TRUE) {
	  if(first == TRUE) {
	    mm[0] = mm[1] = (int)cont[iy][ix]; 
	    first = FALSE;
	  }
	  mm[0] = Min(mm[0],(int)cont[iy][ix]);
	  mm[1] = Max(mm[1],(int)cont[iy][ix]);
	}
  
  return OK;
} // end of cextrema


/****************************************************************************/
int iextrema( int dimx, int dimy, int **cont, int *mm, char **m ) {
  /****************************************************************************/
  /* Computes the extrema of the image */
  /****************************************************************************/

  int ix,iy;
  int first = TRUE;
  mm[0]=65536, mm[1]=-65536;

  if(m == NULL) {
  mm[0]= mm[1]= cont[0][0];
    for( iy=0; iy<dimy; iy++ )
      for( ix=0; ix<dimx; ix++ )       {
	mm[0] = Min(mm[0],cont[iy][ix]);
	mm[1] = Max(mm[1],cont[iy][ix]);
      }
  
  } else
    for( iy=0; iy<dimy; iy++ )
      for( ix=0; ix<dimx; ix++ )      
	if(m[iy][ix] == (char)TRUE) {
	  /* if(first == TRUE) {
	     mm[0] = mm[1] = cont[iy][ix]; 
	     first = FALSE;
	     }  */
	  mm[0] = Min(mm[0],cont[iy][ix]);
	  mm[1] = Max(mm[1],cont[iy][ix]);
	}
  
  return OK;
} // end of iextrema


/****************************************************************************/
int extrema( int dimx, int dimy, double **cont, double *mm, char **m ) {
  /****************************************************************************/
  /* Computes the extrema of the image */
  /****************************************************************************/
  
  int ix,iy;
  int first;
  mm[0] = 1.e+17, mm[1] = -1.e17;

  if(m == NULL) {
    mm[0]= mm[1]= cont[0][0];
    for( iy=0; iy<dimy; iy++ )
      for( ix=0; ix<dimx; ix++ )      {
	mm[0] = fMin(mm[0],cont[iy][ix]);
	mm[1] = fMax(mm[1],cont[iy][ix]);
      }
  
  } else {
    for( iy=0; iy<dimy; iy++ )
      for( ix=0; ix<dimx; ix++ )      
	if(m[iy][ix] == (char)TRUE) {
	  /*if(first == TRUE) {
	     mm[0] = mm[1] = cont[iy][ix]; 
	     first = FALSE;
	     }	*/  
	  mm[0] = fMin(mm[0],cont[iy][ix]);
	  mm[1] = fMax(mm[1],cont[iy][ix]);
	}
  }
  
  return OK;
} // end of extrema

/***************************************************************************/
int extrema_update( int dimx, int dimy, double **in, double *min, double *max, 
		    char **m ) /* dMinMax */{ 
  /***************************************************************************/
  
  int i, j;
  double val;
  
  if(m == NULL)
    for( i=0; i<dimx; i++ )
      for( j=0; j<dimy; j++ ) {
	val = in[j][i];
	*min = (*min<val) ? *min : val;      
	*max = (*max>val) ? *max : val;
      }
  
  else
    for( i=0; i<dimx; i++ )
      for( j=0; j<dimy; j++ ) 
	if(m[j][i] == (char)TRUE) {
	  val = in[j][i];
	  *min = (*min<val) ? *min : val;      
	  *max = (*max>val) ? *max : val;
	}

  return OK;
} // end of extrema_update


/****************************************************************************/
double media( int dimx, int dimy, double **cont, char **m ) {
  /****************************************************************************/
  /* Computes the mean of the image  */
  /****************************************************************************/  
  int ix, iy, count=0;
  double suma=0.;
  
  if(m == NULL) {
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++)
	suma += cont[iy][ix];
    count = dimx*dimy;

  }  else
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++)
	if(m[iy][ix] == (char)TRUE) {
	  suma += cont[iy][ix];
	  count ++;
	}

  suma /= ((double) count);
  
  return suma;
} // end of media


/****************************************************************************/
double media_expon( int dimx, int dimy, double factor, double **cont, char **m ) {
  /****************************************************************************/
  /* Computes the mean of the image  */
  /****************************************************************************/  
  int ix, iy;
  double mean=0.,norm=0.,weight;
  
  if(m == NULL) {
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++) {
	norm += (weight = pow(factor,cont[iy][ix]));
	mean += cont[iy][ix] * weight;
      }
    
  }  else
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++)
	if(m[iy][ix] == (char)TRUE) {
	  norm += (weight = pow(factor,cont[iy][ix]));
	  mean += cont[iy][ix] * weight;
	}
  
  mean /= ((double) count);
  
  return mean;
} // end of media_expon


/****************************************************************************/
double dispersion( int dimx, int dimy, double **cont, char**m ) {
  /****************************************************************************/
  /* Computes the dispersion (standard deviation) of the image  */
  /****************************************************************************/

  int ix, iy, count=0;
  double sumaxx=0.,sumax=0.;

  if(m == NULL) {
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++) {
	sumax += cont[iy][ix];
	sumaxx += cont[iy][ix] * cont[iy][ix];
      }
    count = dimx*dimy;
  } else
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++) 
	if(m[iy][ix] == (char)TRUE) {
	  sumax += cont[iy][ix];
	  sumaxx += cont[iy][ix] * cont[iy][ix];
	  count ++;
	}

  sumax /= ((double) count);
  sumaxx = sumaxx/((double) count)-sumax*sumax;
  sumaxx = sqrt(fabs(sumaxx));

  return sumaxx;
} // end of dispersion


/****************************************************************************/
int mediadisp( int dimx, int dimy, double **cont, double *media, double *disp, 
	       char**m ) {
  /****************************************************************************/
  /* Computes the mean and the dispersion of the image  */
  /****************************************************************************/
  
  int ix, iy, count=0;
  double sumaxx=0.,sumax=0.;
  
  if(m == NULL) {
  for(ix=0;ix<dimx;ix++)
    for(iy=0;iy<dimy;iy++) {
      sumax += cont[iy][ix];
      sumaxx += cont[iy][ix] * cont[iy][ix];
    }
    count = dimx*dimy;
  } else
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++) 
	if(m[iy][ix] == (char)TRUE) {
	  sumax += cont[iy][ix];
	  sumaxx += cont[iy][ix] * cont[iy][ix];
	  count ++;
	}

  sumax /= ((double) count);
  sumaxx = sumaxx/((double) count)-sumax*sumax;
  sumaxx = sqrt(fabs(sumaxx));

  *media = sumax;
  *disp = sumaxx;

  return OK;
} // end of mediadisp

/****************************************************************************/
double dispersion_vec( int dimx, int dimy, double **vx, double **vy, char**m ) {
  /****************************************************************************/

  double dx,dy,out;
  
  dx = dispersion(dimx,dimy,vx,m);
  dy = dispersion(dimx,dimy,vy,m);

  out = sqrt(dx*dx+dy*dy);

  return out;
} // end of dispersion_vec


/***************************************************************************/
double anorma( int dimx, int dimy, double **cont, char**m ) {
  /***************************************************************************/
  /* Substracts the mean of the image to the image  */
  /***************************************************************************/

  double media=0.;
  int ix, iy, count=0;
  
  if(m == NULL) {
    for(iy=0;iy<dimy;iy++)    
      for(ix=0;ix<dimx;ix++)
	media+=cont[iy][ix];
    count = dimx*dimy;

  }   else
    for(iy=0;iy<dimy;iy++)    
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) {
	  media+=cont[iy][ix];
	  count ++;
	}
  
  media /= ((double)count);
  
  if(m == NULL) 
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	cont[iy][ix] -= media;
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) cont[iy][ix] -= media;
  
  return media;
} // end of anorma


/***************************************************************************/
double anorma1( int dimx, int dimy, double **cont, char**m ) {
  /***************************************************************************/
  /* Divides the image by its absolute mean  */
  /***************************************************************************/

  double media=0.;
  int ix, iy, count=0;
  
  if(m == NULL) {
    for(iy=0;iy<dimy;iy++)    
      for(ix=0;ix<dimx;ix++)
	media+=fabs(cont[iy][ix]);
    count = dimx*dimy;

  }   else
    for(iy=0;iy<dimy;iy++)    
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) {
	  media+=fabs(cont[iy][ix]);
	  count ++;
	}
  
  media /= ((double)count);
  
  if(m == NULL) 
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	cont[iy][ix] /= media;
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) cont[iy][ix] /= media;
  
  return media;
} // end of anorma1


/***************************************************************************/
double anorma1_line( int dimx, int dimy, double **data, double theta, char**m )
  /*anorma1_lineas*/{
  /***************************************************************************/
  /* Version incompleta: solo vale para lineas horizontales  y verticales */
  double media=0.;
  int ix, iy;
  
  if(fabs(theta)<0.01) {
    for(iy=0;iy<dimy;iy++)  media += anorma11D( dimx, data[iy], m[iy] );
    media /= (double)dimy;
    
  } else if(fabs(theta-M_PI/2)<0.01) {
    for(ix=0;ix<dimx;ix++) media +=  anorma11D( dimy, data[ix], m[iy] );
    media /= (double)dimx;

  } else media=anorma1(dimx,dimy,data, m);

  return media;
} // end of anorma1_line


/***************************************************************************/
double anorma2( int dimx, int dimy, double **cont, char**m ) {
  /***************************************************************************/
  /* Divides the image by its L2 norma */
  /***************************************************************************/

  double media=0.;
  int ix, iy, count=0;
  
  if(m == NULL) {
    for(iy=0;iy<dimy;iy++)    
      for(ix=0;ix<dimx;ix++)
	media+=cont[iy][ix]*cont[iy][ix];
    count = dimx*dimy;
  }   else
    for(iy=0;iy<dimy;iy++)    
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) {
	  media+=cont[iy][ix]*cont[iy][ix];
	  count ++;
	}
  
  media = sqrt(media / ((double)count));
  
  if(m == NULL) 
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	cont[iy][ix] /= media;
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) cont[iy][ix] /= media;
  
  return media;
} // end of anorma2


/***************************************************************************/
int norma( int dimx, int dimy, double **cont, char** m) {
  /***************************************************************************/
  /* Normalize an image */
  /***************************************************************************/

  int ix,iy;
  double mm[2];

  /* equivalent to the following sequence of operations :
   *       mm=(double*)alloc(2*sizeof(double));
   *       extrema( dimx, dimy, cont, mm, m );
   *       op_shift( dimx, dimy, -mm[0], cont, m ); 
   *       op_scale( dimx, dimy, 1./(mm[1]-mm[0]), cont, m );
   */

  mm[0]= mm[1]= cont[0][0];

  if(m == NULL) {
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)      {
	mm[0]=fMin(mm[0],cont[iy][ix]);
	mm[1]=fMax(mm[1],cont[iy][ix]);
      }
    
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	cont[iy][ix]= (cont[iy][ix] - mm[0]) / (mm[1]-mm[0]);
    
  } else {
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)      
	if(m[iy][ix] == (char)TRUE) {
	  mm[0]=fMin(mm[0],cont[iy][ix]);
	  mm[1]=fMax(mm[1],cont[iy][ix]);
      }
    
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) 
	  cont[iy][ix]= (cont[iy][ix] - mm[0]) / (mm[1]-mm[0]);
  }

  return OK;
} // end of norma


/***************************************************************************/
int denorma( int dimx, int dimy, double norma, double **cont, char**m ) {
  /***************************************************************************/
  /* Adds the value norma (generally, the mean of the image) to the image */
  /***************************************************************************/
  int ix,iy;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	cont[iy][ix] += norma;
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) cont[iy][ix] += norma;
  
  return OK;
} // end of denorma


/************************************************************************/
double extrema_local( double **cont, int *iix, int *iiy, 
			      double mms[2], int winsize ) {
  /************************************************************************/
  /* Calcul des extrema locaux du signal dans la fenêtre dont les coordonnées
   * en x et en y sont respectivement données par iix et iiy.
   * Les coordonnees initialisées à -1 ne sont pas prises en compte dans le
   * calcul des extrema.
   * Remarque: ce code n'est absolument pas optimisé, il envisage simplement 
   * toutes les possibilités de calcul des attributs. En particulier certaines
   * options de calcul paraissent inappropriées. */
  /************************************************************************/
  
  int dx, dy, i=0, j=0;

  /* Recherche du premier indice de pixels que l'on considère (i.e. pour lequel
   * les coordonnées en x et en y sont diff'erents de -1).
   * Remarque: il y en a forcèment un (notamment si l'on considère dans le 
   * cas glcm, la fenêtre des pixels translatés par (tx,ty): ceux-ci peuvent 
   * tous être en dehors des frontières) */

  for( dx=0; dx<2*winsize+1; dx++ ) if( iix[dx]>=0) break;
  for( dy=0; dy<2*winsize+1; dy++ ) if( iiy[dy]>=0 ) break;

  if(dy<2*winsize+1 && dx<2*winsize+1) 
    mms[0] = mms[1] = cont[iiy[dy]][iix[dx]];
  else return 1;
  
  /* Calcul des extrema locaux */
  for( dy=0; dy<2*winsize+1; dy++ )	
    for( dx=0; dx<2*winsize+1; dx++ )     {
      if( iiy[dy]>=0 && iix[dx]>=0) {
	  /* au cas des pixels ne rentreraient pas en compte dans le calcul de 
	   * l'attribut: typ. le cas quand on demande à ce que l'image ne soit pas 
	   * périodisée. 
	   * Ralentit considérablement le calcul...   */
	mms[0] = fMin( mms[0], cont[iiy[dy]][iix[dx]] );
	mms[1] = fMax( mms[1], cont[iiy[dy]][iix[dx]] );
      }      
    }
  
  return (mms[1] - mms[0]);
} // end of extrema_local

