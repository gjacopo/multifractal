#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_alloc.h>
#include <utl_stats.h>


/****************************************************************************/
int eigen2x2( double *m, double *val, double *vec ) {
/****************************************************************************
 * Returns the eigenvalues and eigenvectors of a 2x2 matrix
 * Inputs:
 *    - matrix 2x2 M = (m[0] m[1]
 *                      m[2] m[3])
 * Outputs:
 *    - vector 2x1 of eigenvalues VAL = (val[0] val[1]) where
 *         val[0] = 0.5 * ( sum - det )
 *         val[1] = 0.5 * ( sum + det )
 *      and val[0] <= val[1]
 *      with:
 *            sum = m[0]+m[1]
 *            det = sqrt( 4*m[1]*m[3] + (m[0]-m[1])^2 )
 *    - matrix 2x2 of eigenvectors VEC = (vec[0] vec[1]  where
 *                                        vec[2] vec[3])
 *         vec[0] = cos(theta1) 
 *         vec[1] = sin(theta1)
 *         vec[2] = cos(theta2)
 *         vec[3] = sin(theta2)
 *      with:
 *            theta1 = atan2( val[0]-m[0], m[1])
 *            theta2 = atan2( val[1]-m[0], m[1])
 * 
/****************************************************************************/

  double a, b;
  double theta1, theta2;
    
  a = m[0] + m[3];
  b = a*a - 4(m[0]*m[3] - m[1]*m[2]);
  val[0] = 0.5 * (a+b);
  val[1] = 0.5 * (a-b);

  theta1 = atan2( val[0] - m[0], m[1]);
  theta2 = atan2( val[1] - m[0], m[1]);

  vec[0] = cos(theta1);
  vec[1] = sin(theta1);
  vec[2] = cos(theta2);
  vec[3] = sin(theta2);

  return OK;
}


/****************************************************************************/
int localtensor() {
/****************************************************************************/
 

  return OK;
}

/****************************************************************************/
int ( double alpha, double sigma ) {
  /****************************************************************************/
  
  double A;
  double s_u, s_v;
  double val[2];
  double vec[4];

  A = (val[1] - val[0]) / (val[0] + val[1]);

  /* directional scales s_u and s_v of the anisotropic gaussian kernel */
  s_u = alpha / (alpha+A) * sigma;
  s_v = (alpha+A) / alpha * sigma; 

  /* anisotropic gaussian kernel */
  
  return OK;
}


