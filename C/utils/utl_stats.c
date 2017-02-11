
#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_stats.h>



/***************************************************************************/
double fit( double *x, double *y, int n, 
	    double *a, double *b, double *corr ) {
  /***************************************************************************
   * Least square linear regresssion of y with respect to x 
   ***************************************************************************/ 
  double sumx,sumy;
  double sx,sy,sxy;
  int i;

  sumx = sumy= 0.;
  sumxx = sumxy= 0.;
  sumyy = 0.;
  for( i=0; i<n; i++ ) {
    sumx += x[i];
    sumy += y[i];
    sx += x[i]*x[i];
    sxy += x[i]*y[i];
    sy += y[i]*y[i];
  }
   sumx/=(double)n;
   sumy/=(double)n;
   sxy/=(double)n;
   sx/=(double)n;
   sy/=(double)n;

  sxy -= sumx*sumy;
  sx = sqrt(fabs(sx-sumx*sumx));
  sy = sqrt(fabs(sy-sumy*sumy));
  
  if(sx > MIN_NUM_RESOLUTION) {
    *a = sxy / (sx*sx);
    *b = sumy - (*a)*sumx;
    if(sy<MIN_NUM_RESOLUTION) *corr=1.;
    else                      *corr = sxy / (sx*sy);
    
  } else {
    *a = 0.;
    *b = sumy;
    *corr = 1.;
  }
  
  return sy*sy-a[0]*a[0]*sx*sx;
} // end of fit


/***************************************************************************/
double fit_best(double *x, double *y, int n, int hmin, 
		double **repx, double **repy, int *N) {
  /***************************************************************************/
  
  double *histoc,*rep1,*repe2;
  double mm1[2];
  double norm;
  double m2,me2,var2,vare2,regr;
  double histac,rep1ac,repe2ac;
  int ix,ip,ip0;
  
  *N = n;
  if( (histoc = (double*)calloc(n,sizeof(double))) == NULL ||
      (rep1 = (double*)calloc(n,sizeof(double))) == NULL ||
      (repe2 = (double*)calloc(n,sizeof(double))) == NULL )
    return ERROR; // note that regr cannot be <0
  
  mm1[0] = mm1[1] = x[0];
  for( ix=0; ix<n; ix++ )    {
    mm1[0] = fMin(mm1[0],x[ix]);
    mm1[1] = fMax(mm1[1],x[ix]);
  }
  mm1[1] -= mm1[0];
  
  m2 = var2 = 0.;
  for( ix=0; ix<n; ix++ )  {
    ip = (int) (((double)*N)*(x[ix]-mm1[0])/mm1[1]);
    if(ip >= *N) ip = *N - 1;
    histoc[ip] += 1.;
    rep1[ip] += x[ix];
    repe2[ip] += y[ix];
    m2 += y[ix];
    var2 += y[ix]*y[ix];
  }
  m2 /= (double)n;
  var2 = var2/((double)n) - m2*m2;
  
  histac = 0.;
  rep1ac = repe2ac = 0.;
  ip0 = 0;
  for( ip=0; ip<*N; ip++ )    {
    if(histac < hmin)      {
      histac += histoc[ip];
	  rep1ac += rep1[ip];
	  repe2ac += repe2[ip];
    }  else {
      histoc[ip0] = histac;
      rep1[ip0] = rep1ac;
      repe2[ip0] = repe2ac;
      histac = histoc[ip];
      rep1ac = rep1[ip];
      repe2ac = repe2[ip];
      ip0 ++;
    }
  }

  histoc[ip0] = histac;
  rep1[ip0] = rep1ac;
  repe2[ip0] = repe2ac;
  *N = ip0 + 1;
	
  for( ip=0; ip<*N; ip++ )
    {
      rep1[ip] /= histoc[ip];
      repe2[ip] /= histoc[ip];
    }

  me2 = vare2 = norm = 0.;
  for( ip=0; ip<*N; ip++ )    {
    norm += histoc[ip];
    me2 += repe2[ip] * histoc[ip];
    vare2 += repe2[ip] * repe2[ip] * histoc[ip];
  }
  me2 /= norm;
  vare2 = vare2/norm - me2*me2;
  
  if(var2 > 1e-30) regr = vare2 / var2;
  else regr = 1.;
  
  if( (*repx=(double*)calloc(*N,sizeof(double))) == NULL ||
      (*repy=(double*) calloc(*N,sizeof(double))) == NULL)
    return ERROR;
  
  for( ip=0; ip<*N; ip++ )    {
    repx[0][ip] = rep1[ip];
    repy[0][ip] = repe2[ip];
  }
 
  free(histoc);
  free(rep1);
  free(repe2);

  return(regr);
} // end of fit_best


/***************************************************************************/
int fit_line( double *yValues, double *xValues, int dimx, double *slope )  {
  /***************************************************************************
   * Fit a signal with a straight line by regression, i.e. it computes the
   * slope a of the line y=ax+b that best fits the data [xValues,yValues].
   *
   * Parameters: 
   *      - xValues : the list of indexes of the signal,
   *      - yValues : the values of the signal to fit, 
   *      - dimx : lenght of the signal,
   *      - slope : slope a of the approximation line, what we want
   * Returns the number of scales used in the approximation
   *
   * Note: in LastWave toolbox, this operation is mainly realized by the 
   * function:
   *      LineFitSig
   * in file signal_function.c (package package_signal), except that we
   * consider only signals with regular time intervals.
   ***************************************************************************/

  int i;
  double t,sxoss,sx=0.,sy=0.,st2=0.;
  double a = 0.;
  /* a : the equation line is y = a*x+b */
  double x, y;
  int ss=0;

  for( i=0; i<dimx; i++ )   {
    sx += xValues[i];
    sy += yValues[i];
    ss++;
  }
  sxoss = sx/(double)ss;
  
  for( i=0; i<dimx; i++ )   {
    t = xValues[i] - sxoss;
    st2 += t*t;
    a += t*yValues[i];
  }
  a /= st2;
  *slope = a;
  
  return ss;
} // end of fit_line


/***************************************************************************/
int dimensiona( int dim) {
  /***************************************************************************/

  int out;

  out=1;
  while(out<dim) out=2*out;

  return(out);
} // end of dimensiona


/***************************************************************************/
int adimensiona( int size) {
  /***************************************************************************/

  int out,s0;

  for(out=0,s0=size;s0>1;s0=s0/2) out++;

  return(out);
} // end of adimensiona


/***************************************************************************/
int adimensiona_pos( int size ) { /* strictement positif */
/***************************************************************************/

  int out,s0;
  for( out=0,s0=size; s0>0; s0=s0/2 )	out++;

  return out;
} // end of adimensiona_pos


/***************************************************************************/
double fMax(double a,double b) {
  /***************************************************************************/
  if(a > b) return a;     else return b;
} // end of fMax


/***************************************************************************/
double fMin(double a,double b) {
  /***************************************************************************/
  if(a<b) return a;       else return b;
} // end of fMin


/***************************************************************************/
int Max(int a,int b) {
  /***************************************************************************/
  if(a>b) return a;       else return b;
} // end of Max


/***************************************************************************/
int Min(int a,int b) {
  /***************************************************************************/
  if(a<b) return a;  else return b;
} // end of Min


/***************************************************************************/
int Mod(int a, int b) {
  /***************************************************************************/
  int output = a/b;
  output = a - output*b;
  if(output < 0) output += b;
 
  return output;
} // end of Mod


/***************************************************************************/
int Round( double a) {
  /***************************************************************************/
  int out  = (int)(a+0.5);
  if(a < -0.5) out--;

  return out;
} // end of Round


/***************************************************************************/
void C_mult( double a1, double b1, double a2, double b2, 
	     double *a, double *b ) {
  /***************************************************************************/
  *a = a1*a2 - b1*b2;
  *b = a1*b2 + a2*b1;
} // end of C_mult


/***************************************************************************/
void C_sqrt( double a0, double b0, double *a, double *b) {
  /***************************************************************************/
  double mod = sqrt(a0*a0+b0*b0);

  *a = sqrt(0.5*(mod+a0));
  if(b0<0) *b = -sqrt(0.5*(mod-a0));
  else *b = sqrt(0.5*(mod-a0));
} // end of C_sqrt


/***************************************************************************/
double angulo( double dx, double dy) {
  /***************************************************************************/
  /*  Computes the orientation of any vector */
   /***************************************************************************/

  double salida;
  
  if(fabs(dx)>1e-30)    {
    salida=atan(dy/dx);
    if(dx<0) salida+=M_PI;
    if(salida<0) salida+=2*M_PI;
    if(salida>2*M_PI) salida-=2*M_PI;
  }
  else if(dy>0) salida=M_PI/2;
  else salida=3*M_PI/2;
  
  return(salida);
} // end of angulo


/************************************************************************/
int resolution( int winsize, int limit )  {
  /************************************************************************/
  /* Defining the resolution parameter taking into account the size of the 
   * window. 
   * In our application limit=NBLIMIT=6 */
  /************************************************************************/
  
  int res;
  
  res = (int)(2.*log(2.*winsize+1.)/log(2.)) - 1.;
  /* limiting the resolution anyway to limit bits */
  
  if(limit > 0) res = Min( res, limit ); 
  res = Max( res, 0 );
  res = (int) pow(2.,(double)res);
  
  return res;  
} // end of resolution


/************************************************************************/
int vecindex( int *vecii, int i, int wsize, int dim, int per ) {
  /************************************************************************/
  
  int d, count=0;

  IF(per) {
    /* Version image periodisee: hypothese de cyclite des bords */
    for( d=0; d<2*wsize+1; d++ ) vecii[d] =  /*(i+d-winsize) % dim;*/
				     Mod( i+d-wsize, dim ); 
    count = dim;
    
  } ELSE {
    /* Version image non periodisee */
    for( d=0; d<2*wsize+1; d++ ) 
      if( i-wsize+d<0 || i-wsize+d>=dim)  vecii[d] = -1; 
      else  {
	vecii[d] = i+d-wsize; 
	count ++;
      }
  }
  
  return count; /* Nombre de pixels qui vont entrer dans le calcul */
} // end of vecindex
