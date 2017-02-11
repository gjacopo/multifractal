#include <stdio.h> 
#include <stdlib.h>
#include <math.h>

/* Librairies  ROUTINES */

#include <utils.h>
#include <utl_alloc.h>      
#include <utl_operator.h>

#include <flt_deriva2d.h>
#include <flt_stats1d.h>
#include <flt_stats2d.h>

#include <signal.h>
#include <sig_alloc.h>
#include <sig_operator.h>
#include <sig_gradient.h>


/***************************************************************************/
Gradient *alloc_gradient () {
  /***************************************************************************/
  Gradient *grad=NULL;
  
  if( (grad=(Gradient*)malloc(sizeof(Gradient))) == NULL)
    return NULL;
  else return grad;
} // end of alloc_gradient


/***************************************************************************/
Gradient *create_gradient( Signal* signal ) {
  /***************************************************************************/
  Gradient *grad=NULL;

  /* Note: here the flag parsed to create_signal needs to be set to NULL
   * to avoid infinite recursive call to this function */
  if( (grad=alloc_gradient()) == NULL||
      (grad->gx=create_signal(signal->xdim,signal->ydim,NULL)) == NULL ||
      (grad->gy=create_signal(signal->xdim,signal->ydim,NULL)) == NULL ||
      (grad->normg=create_signal(signal->xdim,signal->ydim,NULL)) == NULL)
    return NULL;
  else 
    return grad;
} // end of create_gradient


/***************************************************************************/
int free_gradient( Gradient* grad ) {
  /***************************************************************************/

  free_signal(grad->gx); grad->gx=NULL;
  free_signal(grad->gy); grad->gy=NULL;
  free_signal(grad->normg); grad->normg=NULL;
  /*
    if(image->gradient->signal != NULL)
    liberar_matriz(signal->gradient->signal,signal->gradient->ydim);    
    Free(signal->gradient);
  */

  return OK;
} // end of free_gradient


/***************************************************************************/
int copy_gradient( Gradient *gin, Gradient* gdest, char **m ) {
  /***************************************************************************/

  if(gdest == NULL)
    TrackNullAlloc( gdest=create_gradient(gin->gx) );

  if(gin->grad != NULL) {
    copy_signal(gin->gx,gdest->gx,m);
    copy_signal(gin->gy,gdest->gy,m);
    copy_signal(gin->normg,gdest->normg,m);
  }
  
  return OK;
} // end of copy_gradient


/***************************************************************************/
int compnaive_gradient( Signal* signal, Gradient *grad ) {
  /***************************************************************************/
  int xdim=signal->xdim, ydim=signal->ydim;

  TrackError( copy(xdim,ydim,signal->signal/*source*/,grad->normg->signal,NULL),
	      "Error copying signal" );
  TrackError( modgradderiva( xdim, ydim, 
		 grad->normg->signal, grad->gx->signal, grad->gy->signal ),
	      "Error deriving signal" );
  
  return OK;
} // end of compnaive_gradient


/***************************************************************************/
int compnaive_normgrad( Signal* signal, Gradient *grad ) {
  /***************************************************************************/
  int xdim=signal->xdim, ydim=signal->ydim;

  TrackError( copy( xdim, ydim, signal->signal/*source*/, grad->normg->signal, NULL ),
	      "Error copying signal" );  
  TrackError( modderiva( xdim, ydim, grad->normg->signal ),
	      "Error deriving signal" );
  
  return OK;
} // end of compnaive_normgrad


/***************************************************************************/
int compute_gradient( Signal* signal, Gradient *grad ) {
  /***************************************************************************/
  double **gx, **gy;
  int xdim=signal->xdim, ydim=signal->ydim;
  int xeff=dimensiona(signal->xdim), yeff=dimensiona(signal->ydim);

  if( (gx=matrix2D(yeff,xeff)) == NULL ||
      (gy=matrix2D(yeff,xeff)) == NULL)
    return ERROR;

  /* fill0(xeff,yeff,gx,NULL);
     fill0(xeff,yeff,gy,NULL);
  */
  /* Compute directional derivatives */
  copy( xdim, ydim, signal->signal/*source*/,gx/*destination*/, NULL);
  deriva( xeff, yeff, gx, gy );
  
  /* Compute modulus */
  copy( xdim, ydim, gx, grad->gx->signal, NULL);
  copy( xdim, ydim, gy, grad->gy->signal, NULL);
  op_sqrtaddsq( xdim, ydim, gx, gy, NULL );
  copy( xdim, ydim, gy, grad->normg->signal, NULL );

  free_matrix2D(gx,yeff);
  free_matrix2D(gy,yeff);

  return OK;
} // end of compute_gradient
