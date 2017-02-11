#ifndef   	SIG_OPERATOR_H_
#define   	SIG_OPERATOR_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

  int compare_dimsignal( Signal *s1, Signal *s2);
  int compare_dimimage( Image *im1, Image *im2);
  
  int copy_signal( Signal *sig, Signal* signal, char **m );
  
  int affine_signal( Signal* s, double scale, double shift, char**m );
  int substract_signal( Signal* sin, Signal* sout, char**m );
  int add_signal( Signal* sin, Signal* sout, char**m );
  int divide_signal( Signal* sin, Signal* sout, char**m );
  int cuadra_signal( Signal* s, char**m );
  
  int norma_signal( Signal* s, char**m );
  int fill_signal( Signal* s, double val, char**m );
  
  /* For all the following tests, the test depends on the sign of the
   * variable sign:
   * if sign>0, check:      sig>thres
   * i.e: if cont>thres, mask=TRUE
   *      else           mask=FALSE
   * if sign<0, check:      -sig>-thres <=> sig<thres
   * i.e. the opposite     */
  int threshold2ptr_signal( Signal* sig, int sign, double thres, double* s );
  int mask2ptr_signal( Signal* sig, int sign, double thres, char* s );
  int mask_signal( Signal* sig, int sign, double thres, char**m );
  
  int vratio_signal( Image* im, Signal *s, int ic1, int ic2, 
		     double c, char**m );
  int vdiffratio_signal( Image* im, Signal *s, int ic1, int ic2, 
			 double c, char**m );
  int rvdiffratio_signal( Image* im, Signal *s, 
			  int ic1, int ic2, int ic3, 
			  double c, char**m );
  
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !SIG_OPERATOR_H_ */
