#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Librairies  ROUTINES */

#include <utils.h>
#include <utl_alloc.h>      
#include <utl_stats.h>	

#include <utl_operator.h>

/***************************************************************************/
int fill0( int dimx, int dimy, double **matriz, char**m )/*limpia*/ {
  /***************************************************************************/
  int ix,iy;

  if(m == NULL)  /* exactly limpia */
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++)
	matriz[iy][ix] = 0.;
  
  else  
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++)
	if(m[iy][ix] == (char)TRUE) matriz[iy][ix] = 0.;
  
  return OK;
}

/***************************************************************************/
int fill1D0( int dimx, double *vector, char*m )/*limpia_lista*/ {
  /***************************************************************************/
  int ix;

  if(m == NULL)  /* exactly limpia_lista */
    for(ix=0;ix<dimx;ix++) vector[ix]=0.;
  
  else  
    for(ix=0;ix<dimx;ix++) if(m[ix] == (char)TRUE) vector[ix] = 0.;
  
  return OK;
}

/***************************************************************************/
int fill( int dimx, int dimy, double d, double **matriz, char**m ) {
  /***************************************************************************/
  int ix,iy;

  if(m == NULL)  
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++)
	matriz[iy][ix] = d;
  
  else  
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++)
	if(m[iy][ix] == (char)TRUE) matriz[iy][ix] = d;
  
  return OK;
}

/***************************************************************************/
int fill1D( int dimx, double d, double *vector, char*m ) {
  /***************************************************************************/
  int ix;

  if(m == NULL)  
    for(ix=0;ix<dimx;ix++) vector[ix] = d;
  
  else  
    for(ix=0;ix<dimx;ix++) if(m[ix] == (char)TRUE) vector[ix] = d;
  
  return OK;
}


/***************************************************************************/
int ifill( int dimx, int dimy, int i, int **matriz, char**m )/*limpia_int*/ {
  /***************************************************************************/
  int ix,iy;
  
  if(m == NULL)
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++)
	matriz[iy][ix] = i;
  
  else
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++)
	if(m[iy][ix] == (char)TRUE) matriz[iy][ix] = i;
  
  return OK;
}


/***************************************************************************/
int ifill1D( int dimx, int i, int *vector, char*m ) {
  /***************************************************************************/
  int ix;
  
  if(m == NULL)
    for(ix=0;ix<dimx;ix++) vector[ix] = i;
  
  else
    for(ix=0;ix<dimx;ix++) if(m[ix] == (char)TRUE) vector[ix] = i;
  
  return OK;
}

/***************************************************************************/
int cfill( int dimx, int dimy, char val, char **matriz, char**m  )/*limpia_char*/ {
  /***************************************************************************/
  int ix,iy;
  
  if(m == NULL)
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++)
	matriz[iy][ix] = val;
  
  else
    for(ix=0;ix<dimx;ix++)
      for(iy=0;iy<dimy;iy++)
	if(m[iy][ix] == (char)TRUE) matriz[iy][ix] = val;
  
  return OK;
}


/***************************************************************************/
int cfill1D( int dimx, char val, char *vector, char*m  )
  /*limpia_char_lista*/ {
  /***************************************************************************/
  int ix;
  
  if(m == NULL)
    for(ix=0;ix<dimx;ix++) vector[ix] = val;
  
  else
    for(ix=0;ix<dimx;ix++)  if(m[ix] == (char)TRUE) vector[ix] = val;
  
  return OK;
}


/***************************************************************************/
int copy( int dimx, int dimy, double **fuente, double **destino, char**m )/*asigna*/ {
  /***************************************************************************/
  int ix,iy;
  
  if(m == NULL) 
    for( iy=0; iy<dimy;iy++ )
      for( ix=0; ix<dimx;ix++ ) 
	destino[iy][ix] = fuente[iy][ix];
     
 else 
    for( iy=0; iy<dimy; iy++ )
      for( ix=0; ix<dimx; ix++ )
	if(m[iy][ix] == (char)TRUE) destino[iy][ix] = fuente[iy][ix];

  return OK;
}

/***************************************************************************/
int copy1D( int dimx, double *fuente, double *destino, char*m )/*asigna_lista*/ {
  /***************************************************************************/
  int ix;
  
  if(m == NULL) 
    for( ix=0; ix<dimx;ix++ ) 
      destino[ix] = fuente[ix];
  
  else 
    for( ix=0; ix<dimx; ix++ )
      if(m[ix] == (char)TRUE) destino[ix] = fuente[ix];
  
  return OK;
}


/***************************************************************************/
int ccopy( int dimx, int dimy, char **fuente, char **destino, char**m ) 
     /*asigna_char*/{
  /***************************************************************************/
  int ix,iy;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	destino[iy][ix] = fuente[iy][ix];
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) destino[iy][ix] = fuente[iy][ix];
  
  return OK;
}

/***************************************************************************/
int coarseres( int dimx, int dimy, double block, double **source,
	       double **target ) /*coarse_resolution*/ {
  /***************************************************************************/
  double bsize;
  int beff,xsize,ysize;
  int ix,iy;
  int bx,by;
  
  if(block<1.) {
    beff=(int) (1./block);
    for(ix=0;ix<dimx;ix++) 
      for(iy=0;iy<dimy;iy++) 
	for(by=0;by<beff;by++) 
	  for(bx=0;bx<beff;bx++) 
	    target[iy*beff+by][ix*beff+bx]=
	      source[iy][ix];
    
  } else {
    beff=(int) block;
    bsize=(double)(beff*beff);
    xsize=dimx/beff;
    ysize=dimy/beff;
    for(iy=0;iy<ysize;iy++) 
      for(ix=0;ix<xsize;ix++) {
	target[iy][ix]=0.;
	for(by=0;by<beff;by++) 
	  for(bx=0;bx<beff;bx++) 
	    target[iy][ix]+=
	      source[beff*iy+by][beff*ix+bx];
		  
	target[iy][ix]=target[iy][ix]/bsize;
      }
		
  }
	
  return OK;
}
/***************************************************************************/
int coarseres1D( int dimx, double block, double *source,
		 double *target ) /*coarse_resolution_lista*/ {
  /***************************************************************************/
  double bsize;
  int beff,xsize;
  int ix, bx;
  
  if(block<1.) {
    beff=(int) (1./block);
    for(ix=0;ix<dimx;ix++) 
      for(bx=0;bx<beff;bx++)  target[ix*beff+bx]=source[ix];
    
  } else {
    beff=(int) block;
    bsize=(double)beff;
    xsize=dimx/beff;
    for(ix=0;ix<xsize;ix++) {
      target[ix]=0.;
      for(bx=0;bx<beff;bx++) target[ix]+= source[beff*ix+bx];
      target[ix]=target[ix]/bsize;
    } 
  }
	
  return OK;
}


/***************************************************************************/
int ccoarseres( int dimx, int dimy, double block, char **source,
		char **target ) /*coarse_resolution_char*/ {
  /***************************************************************************/
  double cumul;
  double bsize;
  int beff,xsize,ysize;
  int ix,iy;
  int bx,by;
  
  if(block<1.) {
    beff=(int) (1./block);
    for(ix=0;ix<dimx;ix++) 
      for(iy=0;iy<dimy;iy++) 
	for(by=0;by<beff;by++) 
	  for(bx=0;bx<beff;bx++) 
	    target[iy*beff+by][ix*beff+bx]=
	      source[iy][ix];
    
  } else {
    beff=(int) block;
    bsize=(double)(beff*beff);
    xsize=dimx/beff;
    ysize=dimy/beff;
    
    for(iy=0;iy<ysize;iy++)	
      for(ix=0;ix<xsize;ix++)	    {
	cumul=0.;
	for(by=0;by<beff;by++)		{
	  for(bx=0;bx<beff;bx++)		    
	    cumul += (double)source[beff*iy+by][beff*ix+bx]
	      +((source[beff*iy+by][beff*ix+bx]<0)?256.:0.);
	}
	target[iy][ix]=(char)(cumul/bsize);
      }
  }
  
  return OK;
}

/***************************************************************************/
int coarseres_mask( int dimx, int dimy, double block, int ivabs,
		     char **source, char **target) /*coarse_resolution_mask*/ {
  /***************************************************************************/

  double bsize;
  double cumul;

  int *voting;
  int beff,xsize,ysize;
  int ix,iy;
  int bx,by;
  int voting0;
  int iv,iv0;

  if(block<1.)    {
    beff=(int) (1./block);
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	for(by=0;by<beff;by++)
	  for(bx=0;bx<beff;bx++)
	    target[iy*beff+by][ix*beff+bx]=source[iy][ix];

  } else   {
    TrackNullAlloc( voting=(int*)calloc(256,sizeof(int)) );
    beff=(int) block;
    bsize=(double)(beff*beff);
    xsize=dimx/beff;
    ysize=dimy/beff;

    for(iy=0;iy<ysize;iy++)      
      for(ix=0;ix<xsize;ix++)	    {
	cumul=0.;
	for(iv=0;iv<256;iv++) voting[iv]=0;
	for(by=0;by<beff;by++)	 
	  for(bx=0;bx<beff;bx++)  {
	    iv=Mod((int)source[beff*iy+by][beff*ix+bx],256);
	    voting[iv]++;
	  }
	iv0=-1;
	voting0=0;
	if(ivabs>0 && voting[ivabs]>0) iv0=ivabs;
	
	if(iv0==-1)
	  for(iv=0;iv<256;iv++)
	    if(voting[iv]>voting0)		{
	      voting0=voting[iv];
	      iv0=iv;
	    }
	target[iy][ix]=(char)iv0;
      }
    
    free(voting);
  }
  
  return OK;
}

/***************************************************************************/
int blurres( int dimx, int dimy, int block, double **source, 
	     double **target ) /*blur_resolution*/ {
  /***************************************************************************/
  double bsize;
  int *iix,*iiy;
  int ix,iy,dx,dy;

  TrackNullAlloc( iix=(int*)calloc(2*block+1,sizeof(int)) );
  TrackNullAlloc( iiy=(int*)calloc(2*block+1,sizeof(int)) );

  bsize=(double)(2*block+1);
  bsize=bsize*bsize;
  for(iy=0;iy<dimy;iy++)  {
    for(dy=0;dy<2*block+1;dy++) iiy[dy]=Mod(iy+dy-block,dimy);
    for(ix=0;ix<dimx;ix++) {
      for(dx=0;dx<2*block+1;dx++)
	iix[dx]=Mod(ix+dx-block,dimx);
      target[iy][ix]=0.;
      for(dy=0;dy<2*block+1;dy++)
	for(dx=0;dx<2*block+1;dx++)
	  target[iy][ix]+=source[iiy[dy]][iix[dx]];
			
      target[iy][ix]=target[iy][ix]/bsize;
    }
  }

  free(iix);
  free(iiy);

  return OK;
}


/***************************************************************************/
int de4(int dimx, int dimy, double **signal ) /*de4*/ {
  /***************************************************************************/
  int ix,iy;
  double **aux;

  TrackNullAlloc( aux=matrix2D(dimy/2,dimx/2) );

  coarseres(dimx,dimy,2.,signal,aux);
  for(iy=0;iy<dimy/2;iy++) 
    for(ix=0;ix<dimx/2;ix++) {
      signal[iy][ix]=aux[iy][ix];
      signal[iy][ix+dimx/2]=aux[iy][ix];
      signal[iy+dimy/2][ix]=aux[iy][ix];
      signal[iy+dimy/2][ix+dimx/2]=aux[iy][ix];
    }

  free_matrix2D(aux,dimy/2);

  return OK;
}

/***************************************************************************/
int de4_mirror(int dimx, int dimy, double **signal ) /**/ {
  /***************************************************************************/
  int ix,iy;
  double **aux;

  TrackNullAlloc( aux=matrix2D(dimy/2,dimx/2) );

  coarseres(dimx,dimy,2.,signal,aux);
  for(iy=0;iy<dimy/2;iy++)
    for(ix=0;ix<dimx/2;ix++) {
      signal[iy][ix]=aux[iy][ix];
      signal[iy][ix+dimx/2]=aux[iy][dimx/2-1-ix];
      signal[iy+dimy/2][ix]=aux[dimy/2-1-iy][ix];
      signal[iy+dimy/2][ix+dimx/2]=aux[dimy/2-1-iy][dimx/2-1-ix];
    }

  free_matrix2D(aux,dimy/2);

  return OK;
}


/***************************************************************************/
int paddwindow( int dimx, int dimy, int ix0, int iy0, int dx, 
		int dy, double **data, double **window ) /*recorta_ventana*/ {
  /***************************************************************************/
  int *iix,*iiy;
  int ix,iy;

  TrackNullAlloc( iix=(int*)calloc(dx,sizeof(int)) );
  TrackNullAlloc( iiy=(int*)calloc(dy,sizeof(int)) );

  for(ix=0;ix<dx;ix++) iix[ix] = Mod(ix0-dx/2+ix,dimx);
  for(iy=0;iy<dy;iy++) iiy[iy] = Mod(iy0-dy/2+iy,dimy);

  for(iy=0;iy<dy;iy++)
    for(ix=0;ix<dx;ix++)
      window[iy][ix]=data[iiy[iy]][iix[ix]];

  free(iix);
  free(iiy);

  return OK;
}


/***************************************************************************/
int cpaddwindow( int dimx, int dimy, int ix0, int iy0, int dx, 
		  int dy, char **data, char **window)/*recorta_ventana_char*/ {
/***************************************************************************/
  int *iix,*iiy;
  int ix,iy;
  
  TrackNullAlloc( iix=(int*)calloc(dx,sizeof(int)) );
  TrackNullAlloc( iiy=(int*)calloc(dy,sizeof(int)) );

  for(ix=0;ix<dx;ix++) iix[ix]=Mod(ix0+ix,dimx);
  for(iy=0;iy<dy;iy++) iiy[iy]=Mod(iy0+iy,dimy);

  for(iy=0;iy<dy;iy++)
      for(ix=0;ix<dx;ix++)
	  window[iy][ix]=data[iiy[iy]][iix[ix]];

  free(iix);
  free(iiy);

  return OK;	  
}

/***************************************************************************/
int unpaddwindow( int dimx, int dimy, int ix0, int iy0, int ix1, 
		  int iy1, int tx, int ty, double **window, double **total ) 
     /*remete_ventana*/ {
  /***************************************************************************/
  int *iix,*iiy;
  int ix,iy;

  TrackNullAlloc( iix=(int*)calloc(tx,sizeof(int)) );
  TrackNullAlloc( iiy=(int*)calloc(ty,sizeof(int)) );

  for(ix=0;ix<tx;ix++) iix[ix]=Mod(ix0-tx/2+ix,dimx);
  for(iy=0;iy<ty;iy++) iiy[iy]=Mod(iy0-ty/2+iy,dimy);

  for(iy=0;iy<ty;iy++)
    for(ix=0;ix<tx;ix++)
      total[iiy[iy]][iix[ix]]=window[iy1+iy][ix1+ix];
	
  return OK;
}


/***************************************************************************/
int op_diff( int dimx, int dimy, double **fuente, double **destino, char**m ) 
     /*asigna_resta*/{
  /***************************************************************************/
  int ix,iy;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	destino[iy][ix] -= fuente[iy][ix];
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) destino[iy][ix] -= fuente[iy][ix];
  
  return OK;
}


/***************************************************************************/
int op_diff1D( int dimx, double *fuente, double *destino, char*m ) 
     /*asigna_resta_lista*/{
  /***************************************************************************/
  int ix;
  
  if(m == NULL)
      for(ix=0;ix<dimx;ix++)
	destino[ix] -= fuente[ix];
  
  else
      for(ix=0;ix<dimx;ix++)
	if(m[ix] == (char)TRUE) destino[ix] -= fuente[ix];
  
  return OK;
}


/***************************************************************************/
int op_add( int dimx, int dimy, double **fuente, double **destino, char**m ) 
     /*asigna_suma*/{
  /***************************************************************************/
  int ix,iy;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	destino[iy][ix] += fuente[iy][ix];
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) destino[iy][ix] += fuente[iy][ix];
  
  return OK;
}


/***************************************************************************/
int op_combine( int dimx, int dimy, double escalar, double **fuente,
		double **destino, char**m ) /*asigna_combina*/{
  /***************************************************************************/
  int ix,iy;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	destino[iy][ix] += escalar*fuente[iy][ix];
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) 
	  destino[iy][ix] += escalar*fuente[iy][ix];
  
  return OK;
}


/***************************************************************************/
int op_shift( int dimx, int dimy, double shift, double **funcion, char**m ) 
     /*asigna_desplaza*/{
  /***************************************************************************/
  
  int ix,iy;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	funcion[iy][ix] += shift;
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) funcion[iy][ix] += shift;
  
  return OK;
}

/***************************************************************************/
int op_scale( int dimx, int dimy, double scale, double **func, char**m ) 
     /*asigna_escala*/{
  /***************************************************************************/

  int ix,iy;

  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	func[iy][ix] *= scale;
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) func[iy][ix] *= scale;
  
  return OK;
}

/***************************************************************************/
int cop_scale( int dimx, int dimy, int scale, char **funcion, char**m ) {
  /***************************************************************************/

  int ix,iy;

  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	funcion[iy][ix] = (char)(scale*(int)funcion[iy][ix]);
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) 
	  funcion[iy][ix] = (char)(scale*(int)funcion[iy][ix]);
  
  return OK;
}



/***************************************************************************/
int op_fabs( int dimx, int dimy, double **fuente, double **destino, char**m ) {
  /***************************************************************************/
  int ix,iy;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	destino[iy][ix] = fabs(fuente[iy][ix]);
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) destino[iy][ix] = fabs(fuente[iy][ix]);
  
  return OK;
}


/***************************************************************************/
int op_divide( int dimx, int dimy, double **fuente, double **destino, 
	       char**m ) /*asigna_divide*/{
  /***************************************************************************/
  int ix,iy;

  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(fabs(fuente[iy][ix])>1e-20)
	  destino[iy][ix] /= fuente[iy][ix];
	else destino[iy][ix]=0.;
  
  else
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)
      if(m[iy][ix] == (char)TRUE) {
	if(fabs(fuente[iy][ix])>1e-20)
	  destino[iy][ix] /= fuente[iy][ix];
	else destino[iy][ix]=0.;
      }
  
  return OK;
}

/***************************************************************************/
int op_divide1D( int dimx, double *fuente, double *destino, 
	       char*m ) /*asigna_divide_lista*/{
  /***************************************************************************/
  int ix;

  if(m == NULL)
    for(ix=0;ix<dimx;ix++)
      if(fabs(fuente[ix])>1e-20) destino[ix] /= fuente[ix];
      else destino[ix]=0.;
  
  else
    for(ix=0;ix<dimx;ix++)
      if(m[ix] == (char)TRUE) {
	if(fabs(fuente[ix])>1e-20)  destino[ix] /= fuente[ix];
	else destino[ix]=0.;
      }
  
  return OK;
}

/***************************************************************************/
int op_multiply( int dimx, int dimy, double **fuente, double **destino, 
		 char**m ) /*asigna_multiplica*/{
  /***************************************************************************/
  int ix,iy;

  if(m == NULL)
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)
      destino[iy][ix]=destino[iy][ix]*fuente[iy][ix];
  
  else
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++)
      if(m[iy][ix] == (char)TRUE) 
	destino[iy][ix]=destino[iy][ix]*fuente[iy][ix];

  return OK;
}

/***************************************************************************/
int op_multiply1D( int dimx, double *fuente, double *destino, char*m ) /*asigna_multiplica_lista*/{
  /***************************************************************************/
  int ix;

  if(m == NULL)
    for(ix=0;ix<dimx;ix++)
      destino[ix]=destino[ix]*fuente[ix];
  
  else
    for(ix=0;ix<dimx;ix++)
      if(m[ix] == (char)TRUE) 
	destino[ix]=destino[ix]*fuente[ix];

  return OK;
}

/***************************************************************************/
int opvec_divide(int dimx, int dimy, double **fx, double **fy, 
		       double **dx, double **dy, char**m ) /*asigna_divide_vec*/ {
  /***************************************************************************/
  /* Computes the complex division of two vectorial field seen as components 
   * in the complex plane. */
  /***************************************************************************/
  
  double buffx,buffy,mod;
  int ix,iy;

  if(m == NULL)
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++) {
      mod=fx[iy][ix]*fx[iy][ix]+fy[iy][ix]*fy[iy][ix];
      if(mod>1e-30) {
	buffx=(fx[iy][ix]*dx[iy][ix]+
	       fy[iy][ix]*dy[iy][ix])/mod;
	buffy=(fx[iy][ix]*dy[iy][ix]-
	       fy[iy][ix]*dx[iy][ix])/mod;
	dx[iy][ix]=buffx;
	dy[iy][ix]=buffy;
      } else 
	dx[iy][ix]= dy[iy][ix]= 0.;
    }
  
  else
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++) 
      if(m[iy][ix] == (char)TRUE) {
	mod=fx[iy][ix]*fx[iy][ix]+fy[iy][ix]*fy[iy][ix];
	if(mod>1e-30) {
	  buffx=(fx[iy][ix]*dx[iy][ix]+
		 fy[iy][ix]*dy[iy][ix])/mod;
	  buffy=(fx[iy][ix]*dy[iy][ix]-
		 fy[iy][ix]*dx[iy][ix])/mod;
	  dx[iy][ix]=buffx;
	  dy[iy][ix]=buffy;
	} else 
	  dx[iy][ix]= dy[iy][ix]= 0.;
      }
  
  return OK;
}

/***************************************************************************/
int opvec_multiply( int dimx, int dimy, double **fx, double **fy, 
		    double **dx, double **dy, char**m ) /*asigna_multiplica_vec*/ {
  /***************************************************************************/
  double buffx,buffy;
  int ix,iy;

  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)	{
	buffx = fx[iy][ix]*dx[iy][ix] - fy[iy][ix]*dy[iy][ix];
	buffy = fx[iy][ix]*dy[iy][ix] + fy[iy][ix]*dx[iy][ix];
	dx[iy][ix] = buffx;
	dy[iy][ix] = buffy;
    }
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)	
	if(m[iy][ix] == (char)TRUE) {
	  buffx = fx[iy][ix]*dx[iy][ix] - fy[iy][ix]*dy[iy][ix];
	  buffy = fx[iy][ix]*dy[iy][ix] + fy[iy][ix]*dx[iy][ix];
	  dx[iy][ix] = buffx;
	  dy[iy][ix] = buffy;
	}
  
  return OK;
}

/***************************************************************************/
int op_modulus( int dimx, int dimy, double **vx, double **vy,
		double **mod, char**m ) /*asigna_modulo*/ {
  /***************************************************************************/
  int ix,iy;

  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	mod[iy][ix] = sqrt(vx[iy][ix]*vx[iy][ix]
			   +vy[iy][ix]*vy[iy][ix]);
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) mod[iy][ix] = sqrt(vx[iy][ix]*vx[iy][ix]
						       +vy[iy][ix]*vy[iy][ix]);
  
  return OK;
}

/***************************************************************************/
int opvec_norma( int dimx, int dimy, double **vx, double **vy, char**m ) 
     /*normaliza_vector*/ {
  /***************************************************************************/
  double mod;
  int ix,iy;

  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) {
	mod=sqrt(vx[iy][ix]*vx[iy][ix]+vy[iy][ix]*vy[iy][ix]);
	if(mod>1e-30) {
	  vx[iy][ix]=vx[iy][ix]/mod;
	  vy[iy][ix]=vy[iy][ix]/mod;
	} else
	  vx[iy][ix]= vy[iy][ix]= 0.;
      }
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) 
	if(m[iy][ix] == (char)TRUE) {
	  mod=sqrt(vx[iy][ix]*vx[iy][ix]+vy[iy][ix]*vy[iy][ix]);
	  if(mod>1e-30) {
	    vx[iy][ix]=vx[iy][ix]/mod;
	    vy[iy][ix]=vy[iy][ix]/mod;
	  } else
	    vx[iy][ix]= vy[iy][ix]= 0.;
	}
  
  return OK;
}

/***************************************************************************/
int opvec_rotate( int dimx, int dimy, double angle, double **vx, double **vy, 
		  char**m ) /*gira*/ {
  /***************************************************************************/
  int ix,iy;
  double normax,normay;
  double buffx,buffy;

  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)  {
	buffx=vx[iy][ix]*cos(angle)-vy[iy][ix]*sin(angle);
	buffy=vx[iy][ix]*sin(angle)+vy[iy][ix]*cos(angle);
	vx[iy][ix]=buffx;
	vy[iy][ix]=buffy;
      }
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)  
	if(m[iy][ix] == (char)TRUE) {
	  buffx=vx[iy][ix]*cos(angle)-vy[iy][ix]*sin(angle);
	  buffy=vx[iy][ix]*sin(angle)+vy[iy][ix]*cos(angle);
	  vx[iy][ix]=buffx;
	  vy[iy][ix]=buffy;
	}
  
  return OK;
}


/***************************************************************************/
int op_affine( int dimx, int dimy, double scale, double escalar, 
	       double **funcion, char**m ) 
     /*asigna_afin*/ {
  /***************************************************************************/
  
  int ix,iy;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	funcion[iy][ix] = scale * funcion[iy][ix] + escalar;
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)
	if(m[iy][ix] == (char)TRUE) 
	  funcion[iy][ix] = scale * funcion[iy][ix] + escalar;
  
  return OK;
}


/***************************************************************************/
int op_addsq( int dimx, int dimy, double **fuentes, double **destino, char**m ) 
     /*asigna_multiplica_cuadra*/ {
  /***************************************************************************/
  int ix,iy;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) 
	destino[iy][ix] = fuentes[iy][ix]*fuentes[iy][ix]+ 
	  destino[iy][ix]*destino[iy][ix];
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) 
	if(m[iy][ix] == (char)TRUE) 
	  destino[iy][ix] = fuentes[iy][ix]*fuentes[iy][ix]+ 
	    destino[iy][ix]*destino[iy][ix];
  
  return OK;
}

/***************************************************************************/
int op_sqrtaddsq( int dimx, int dimy, double **fuentes, double **destino,
		  char**m )      /*asigna_multiplica_racina*/ {
  /***************************************************************************/
  int ix,iy;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) 
	destino[iy][ix] = sqrt( fuentes[iy][ix]*fuentes[iy][ix]+ 
				destino[iy][ix]*destino[iy][ix] );
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) 
	if(m[iy][ix] == (char)TRUE) 
	  destino[iy][ix] = sqrt( fuentes[iy][ix]*fuentes[iy][ix]+ 
				  destino[iy][ix]*destino[iy][ix] );
  
  return OK;
}

/***************************************************************************/
int op_square( int dimx, int dimy, double **fuentes, double **destino, 
	       char**m )      /*asigna_cuadra*/ {
  /***************************************************************************/
  int ix,iy;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) 
	destino[iy][ix] = fuentes[iy][ix]*fuentes[iy][ix];
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) 
	if(m[iy][ix] == (char)TRUE) 
	  destino[iy][ix] = fuentes[iy][ix]*fuentes[iy][ix];
  
  return OK;
}



/***************************************************************************/
int op_logaddsq( int dimx, int dimy, double **fuentes, double **destino, 
		 char**m )      /*asigna_log_multiplica_cuadra*/ {
  /***************************************************************************/
  int ix,iy;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) {
	destino[iy][ix] = fuentes[iy][ix]*fuentes[iy][ix]+ 
	  destino[iy][ix]*destino[iy][ix] ;
	if(destino[iy][ix] > 1.e-17 )	destino[iy][ix] = log(destino[iy][ix]);
	else destino[iy][ix] = 0.;
      }
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) 
	if(m[iy][ix] == (char)TRUE) {
	  destino[iy][ix] = fuentes[iy][ix]*fuentes[iy][ix]+ 
	    destino[iy][ix]*destino[iy][ix] ;
	  if(destino[iy][ix] > 1.e-17 )	destino[iy][ix] = log(destino[iy][ix]);
	  else destino[iy][ix] = 0.;
	}
  
  return OK;
}



/***************************************************************************/
int op_log( int dimx, int dimy, double **fuentes, double **destino, char**m ) 
     /*asigna_log*/ {
  /***************************************************************************/
  int ix,iy;
  int nul=0;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) 
	if(destino[iy][ix] > 1.e-17 )	destino[iy][ix] = log(fuentes[iy][ix]);
	else {
	  destino[iy][ix] = 0.;
	  nul++;
	}
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) 
	if(m[iy][ix] == (char)TRUE) {
	  if(destino[iy][ix] > 1.e-17 )	
	    destino[iy][ix] = log(fuentes[iy][ix]);
	else {
	  destino[iy][ix] = 0.;
	  nul++;
	}
	}
  
  return OK;
}



/***************************************************************************/
int opvec_sqdiff( int dimx, int dimy, 
		  double **ux, double **uy, double **vx, double **vy,
		  double **destino, char**m ) /*asigna_cuadradiff_vec*/ {
  /***************************************************************************/
  int ix,iy;

  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) {
	destino[iy][ix]= (vx[iy][ix]-ux[iy][ix])*(vx[iy][ix]-ux[iy][ix])
	  + (vy[iy][ix]-uy[iy][ix])*(vy[iy][ix]-uy[iy][ix]);
	destino[iy][ix]=sqrt(destino[iy][ix]);
      }
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++) 
	if(m[iy][ix] == (char)TRUE) {
	  destino[iy][ix]= (vx[iy][ix]-ux[iy][ix])*(vx[iy][ix]-ux[iy][ix])
	    + (vy[iy][ix]-uy[iy][ix])*(vy[iy][ix]-uy[iy][ix]);
	  destino[iy][ix]=sqrt(destino[iy][ix]);
	}
  
  return OK;
}


/***************************************************************************/
int opvec_orient( int dimx, int dimy, double **ax, double **ay,
		  double **gx, double **gy, double **orient, char **m ) 
     /*asigna_orientation*/{
  /***************************************************************************/
  
  int ix,iy;
  double prod, gmod, mod;
  
  if(m == NULL)
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)	{
	/* Calcul des normes de (gx,gy) et (ax,ay) */
	gmod = sqrt(gx[iy][ix]*gx[iy][ix] + gy[iy][ix]*gy[iy][ix]);
	mod = sqrt(ax[iy][ix]*ax[iy][ix] + ay[iy][ix]*ay[iy][ix]) * gmod;
	/* Calcul du produit scalaire de ces deux vecteurs */
	prod = gx[iy][ix]*ax[iy][ix] + gy[iy][ix]*ay[iy][ix];
	if(fabs(mod)>1.e-17 ) prod /= mod;
	/* Calcul de l'angle entre ces deux vecteurs */
	orient[iy][ix] = acos( prod );
	/* The acos() function returns the arc cosine in radians and the value is  
	 * mathematically defined to be between 0 and PI (inclusive).  */ 
      }
  
  else
    for(iy=0;iy<dimy;iy++)
      for(ix=0;ix<dimx;ix++)	
	if(m[iy][ix] == (char)TRUE) {
	  gmod = sqrt(gx[iy][ix]*gx[iy][ix] + gy[iy][ix]*gy[iy][ix]);
	  mod = sqrt(ax[iy][ix]*ax[iy][ix] + ay[iy][ix]*ay[iy][ix]) * gmod;
	  prod = gx[iy][ix]*ax[iy][ix] + gy[iy][ix]*ay[iy][ix];
	  if(fabs(mod)>1.e-17 ) prod /= mod;
	  orient[iy][ix] = acos( prod );
	}
  
  return OK;
}


/***************************************************************************/
int op_binary( int dimx, int dimy, double **cont, char** mask,
	       int sign, double thres ) {
  /***************************************************************************/
  /* if sign>0, check:      cont>thres
   * i.e: if cont>thres, then mask=TRUE
   *      else                mask=FALSE
   * if sign<0, check:      -cont>-thres <=> cont<thres
   * i.e. the opposite    
   */

  int ix,iy;
  double prod, gmod, mod;
  int N=0;
  
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++) 
      if((double)sign*cont[iy][ix] > (double)sign*thres) {
	mask[iy][ix] = (char)TRUE; 
	N++;
      } else 	mask[iy][ix] = (char)FALSE;
  
  return N;
}


/***************************************************************************/
int op2ptr_binary( int dimx, int dimy, double **cont, char* c,
		   int sign, double thres ) {
  /***************************************************************************/
  
  int ix,iy;
  double prod, gmod, mod;
  int N=0;
  
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++) 
      if((double)sign*cont[iy][ix] > (double)sign*thres)  {
	N++;
	c[iy*dimx+ix] = (char)TRUE;
      } else c[iy*dimx+ix] = (char)FALSE;
  
  return N;
}


/***************************************************************************/
int op2ptr_threshold( int dimx, int dimy, double **cont, double* c,
		      int sign, double thres ) {
/***************************************************************************/
  
  int ix,iy;
  double prod, gmod, mod;
  int N=0;
  
  for(iy=0;iy<dimy;iy++)
    for(ix=0;ix<dimx;ix++) 
      if((double)sign*(c[N]=cont[iy][ix]) > (double)sign*thres)   N++;

  return N;
}

/***************************************************************************/
int opbin_and( int dimx, int dimy, int sign, 
	       char **fuen, char **dest ) {
  /***************************************************************************/
  /* if sign>0, check when both conditions are TRUE
   * i.e: if c1=c2=TRUE, then c2=TRUE
   *      else                c2=FALSE
   * if sign<0, check when one of the 2 conditions is FALSE
   * i.e: if c1=FALSE or c2=FALSE, then c2=TRUE
   *      else                          c2=FALSE
   */
  int N=0;
  char res;
  int ix, iy;

  if(sign>0) res=(char)TRUE; /* 1 */
  else       res=(char)FALSE; /* 0 */
  
  for( iy=0; iy<dimy; iy++ )
    for( ix=0; ix<dimx; ix++ ) 
      if((char)(fuen[iy][ix]*dest[iy][ix]) == res) {
	dest[iy][ix] = (char)TRUE;
	N++; 
      } else dest[iy][ix] = (char)FALSE;

  return N;
}


/***************************************************************************/
int op_rdiffratio( double **i1, double**i2, double **i3,
		    int dimx, int dimy,
		    double **o1, double **o2, double **o3 ) {
  /***************************************************************************/
  
  int ix, iy;
  double a, b;
  int flag_warning=TRUE;

  for( ix=0; ix<dimx; ix++ )
    for( iy=0; iy<dimy; iy++ ) {
      o1[iy][ix] = (i1[iy][ix] - i2[iy][ix]) / (i1[iy][ix] + i2[iy][ix]);
      o2[iy][ix]  = (i1[iy][ix] - i3[iy][ix]) / (i1[iy][ix] + i3[iy][ix]);
      if(i3[iy][ix] < 1.e-17) {
	flag_warning=FALSE;
	o3[iy][ix] = -1.1;
      } else {
	a = i1[iy][ix] / i3[iy][ix];
	b = i2[iy][ix] / i3[iy][ix];
	o3[iy][ix] = (a - b) / (a + b);
      }
    }

  return flag_warning;
}



/***************************************************************************/
int cop_find( int dimx, int dimy, char **m, int nvalues, int *mm) {
  /***************************************************************************/
  
  int ix, iy;
  int c=0, i, flag;
  
  mm[0] = m[0][0];

  for( iy=0; iy<dimy; iy++ )
    for( ix=0; ix<dimx; ix++ ) {
      flag=FALSE;
      for( i=0; i<=c; i++ ) 
	if(m[iy][ix] == (char)mm[i]) { flag=TRUE; break; }
      if(flag == FALSE)  { 
	c++; /* one more value encountered */
	if(c >= nvalues) return FALSE; /* we found more than NVALUES values
				       * in the buffer */
	else  mm[c] = (int)m[iy][ix];
      }
    }
  
  return c+1;  /* we found NVALUES values or less in the buffer */
}


/***************************************************************************/
int op_find( int dimx, int dimy, double **cont, int nvalues, double *mm) {
  /***************************************************************************/
  
  int ix, iy;
  int c=0, i, flag;
  
  mm[0] = cont[0][0];
  
  for( iy=0; iy<dimy; iy++ )
    for( ix=0; ix<dimx; ix++ ) {
      flag=FALSE;
      for( i=0; i<=c; i++ ) 
	if(cont[iy][ix] == mm[i]) { flag=TRUE; break; }
      if(flag == FALSE)  { 
	c++;
	if(c >= nvalues) return FALSE; /* we found more than NVALUES values
				       * in the buffer */
	else  mm[c] = cont[iy][ix];
      }
    }
  
  return (c+1);  /* we found (c+1)<=NVALUES values in the buffer */
}


/***************************************************************************/
int cop_change( int dimx, int dimy, char **cont, int i0, int i1, char **m ) {
  /***************************************************************************/
  int ix, iy;
  int c=0;

  if(m != NULL) {
    for( iy=0; iy<dimy; iy++ )
      for( ix=0; ix<dimx; ix++ ) 
	if(m[iy][ix]==TRUE && cont[iy][ix]==(char)i0) cont[iy][ix] = (char)i1;
  
  } else
    for( iy=0; iy<dimy; iy++ )
      for( ix=0; ix<dimx; ix++ )
	if(cont[iy][ix] == (char)i0) cont[iy][ix] = (char)i1;
  
  return OK;
}


/***************************************************************************/
int op_change( int dimx, int dimy, double **cont, double i0, double i1, char **m ) {
  /***************************************************************************/
  int ix, iy;
  int c=0;
  
  if(m != NULL) {
    for( iy=0; iy<dimy; iy++ )
      for( ix=0; ix<dimx; ix++ ) 
	if(m[iy][ix]==TRUE && cont[iy][ix]==i0) cont[iy][ix] = i1;
    
  } else
    for( iy=0; iy<dimy; iy++ )
      for( ix=0; ix<dimx; ix++ )
	if(cont[iy][ix] == i0) cont[iy][ix] = i1;
  
  return OK;
}
