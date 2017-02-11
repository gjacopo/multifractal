#ifndef   	SIG_GRADIENT_H_
#define   	SIG_GRADIENT_H_

typedef struct mygradient {
  struct mysignal *gx; 
  struct mysignal *gy;
  struct mysignal *normg;
} Gradient;

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */
  
  Gradient *alloc_gradient ();
  Gradient *create_gradient( Signal* signal );
  int free_gradient( Gradient *grad );
  int copy_gradient( Gradient *gin, Gradient* gdest, char **m );
  
  int compute_gradient( Signal* signal, Gradient *grad );
  int compnaive_gradient( Signal* signal, Gradient *grad );
  int compnaive_normgrad( Signal* signal, Gradient *grad );
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !SIG_GRADIENT_H_ */
