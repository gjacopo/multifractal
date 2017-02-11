#include <stdio.h>
#include <math.h>

#include <utils.h>
#include <utl_alloc.h>
#include <utl_stats.h>
#include <utl_op.h>
		
#include <flt_.h>		
#include <flt_fft2d.h>		
#include <flt_deriva2d.h>		
#include <flt_wavelet2d.h>

#include <inout.h>  	
#include <io_grafic.h>  	

#include <mf_manifold.h>

#ifdef _PARSE_FRACTAL_PARAMETERS_
#include <mf_parse.h>
extern ParFRAC *p_frac;
#endif/* !_PARSE_FRACTAL_PARAMETERS_ */


/* outside
   TrackNullAlloc( prob=(double*)calloc(nbox+1,sizeof(double)) );  
   for(ip=0;ip<=nbox;ip++) prob[ip]=0.;
   expon_density(...,prob,...);
   write_expon_density(...,prob,...);
   outside  free(prob);
*/
/***************************************************************************/
int expon_density( int dimx, int dimy, double *mm, double **expon, 
		   double *prob, double *h_inf, double *delta_h ) {
  /***************************************************************************/
  
  double h,h_l, h_h, quantil;
  int ix,iy,ip;  
  int  N=0;
  
  int nbox;
  double hmin, hmax;
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  nbox = p_frac->nbox;
  mm[0] = hmin = p_frac->hmin, mm[1] = hmax = p_frac->hmax;
#else
  nbox = NBOX;
  mm[0] = hmin = HMIN, mm[1] = hmax = HMAX;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
  
  
  for(iy=0;iy<dimy;iy++)    
    for(ix=0;ix<dimx;ix++)      {
      h = (expon[iy][ix]-hmin)*((double)nbox)/(hmax-hmin);
      if((h>=-0.5) && (h<=nbox))	{
	N++;
	ip = (int) h;
	prob[ip] += 1.;
      }
    }
  
  if(N >= 1)    {
    h_l = h_h = hmin - 1.;
    quantil = 0.;
    for( ip=0; (ip<=nbox)&&(h_h==hmin-1); ip++ )      {
      quantil += prob[ip]/(double)N;
      if((quantil>MSMMIN) && (h_l==hmin-1)) 
	h_l = hmin+((double)ip)/nbox*(hmax-hmin);
      if((quantil>MSMMAX) && (h_h==hmin-1)) 
	h_h = hmin+((double)ip)/nbox*(hmax-hmin);
    }
    if(h_l < h_h)      {
      *h_inf = (h_l+h_h)/2.;
      *delta_h = (h_h-h_l)/2.;
    }    else      {
#ifdef _PARSE_FRACTAL_PARAMETERS_
      *h_inf = p_frac->hinf;
      *delta_h = p_frac->dh;
#else
      *h_inf = HINF;
      *delta_h = DH;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
    }
    
    IFVERBOSE 
      WarningVV("MSM at h=%f with a dispersion of %f\n\n",*h_inf,*delta_h);
    
  }  else
    Error("The chosen interval of exponents is empty!!");
  
  return N;
} // end of expon_density


  /* outside
     TrackNullAlloc( prob=(double*)calloc(nbox+1,sizeof(double)) );  
     for(ip=0;ip<=nbox;ip++) prob[ip]=0.;
     expon_density_dif(...,prob,...);
     write_expon_density(...,prob,...);
     outside  free(prob);
  */
/***************************************************************************/
int expon_density_dif( int dimx, int dimy, double *mm, double **expon,
		       int i, double *valor, double *prob, double *disp) {
  /***************************************************************************/
  double h, h_l, h_h, quantil;
  double valesp= 0.,dis= 0.;
  int ix,iy,ip;
  int N=dimx*dimy;
  int nbox;
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  nbox = p_frac->nbox;
  hmin = p_frac->hmin, hmax = p_frac->hmax;
#else
  nbox = NBOX;
  hmin = HMIN, hmax = HMAX;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
  
  mm[0] = mm[1] = expon[0][0];
  for(iy=0;iy<dimy;iy++)    {
    for(ix=0;ix<dimx;ix++)	{
      mm[1] = fMax(maximo,expon[iy][ix]);
      mm[0] = fMin(minimo,expon[iy][ix]);
    }
  }
  
  for(iy=0;iy<dimy;iy++)    {
    for(ix=0;ix<dimx;ix++)      {
      h = (expon[iy][ix] - mm[0])*nbox/(mm[1] - mm[0]);
      valesp += expon[iy][ix];
      dis = valesp + expon[iy][ix]*expon[iy][ix];
      ip = (int) h;
      prob[ip] += 1.;
    }
  }
  valesp /= (dimx*dimy);
  dis = sqrt(dis/(dimx*dimy) - valesp*valesp);
  
  valor[i] = valesp;
  disp[i] = dis;
  
  return N;
} // end of expon_density_dif


  /* outside  
     char name[MAXCHARLENGTH];
     sprintf(name,"Density_dif_h.%s.%s",ext,nombre_wv);
     sprintf(name,"Density_h.%s.%s",ext,name_in);
  */


  /* outside:	Wavelet parameters	
     orden[0]=1.;  a0[0]=.33;  quanto[0]=1.25; 
     orden[1]=1.5;  a0[1]=.75;  quanto[1]=1.25;
     orden[2]=2.; a0[2]=1.;  quanto[2]=1.25;
     orden[3]=3.;  a0[3]=2.5;  quanto[3]=1.05;  
     iwav=1; 
     sc0=a0[iwav], quanto=quanto[iwav], orden=orden[iwav], ord_der=0
     TrackNullAlloc( histoh=matrix2D(3,NBOX+1) );
     [...]
     multifractal_histo( dimx, dimy, sc0, quanto, ord_der, orden, 
     msm, signal, expon, histoh );
     [...]
     free_matrix2D(histoh,3);                                   x  */
/***************************************************************************/
int multifractal_histo( int dimx, int dimy,
			double sc0, double quanto, int ord_der, double orden, 
			char **msm, double **signal, double **expon, 
			double **histoh ) {
  /***************************************************************************/
  /* Computes the multifractal exponents with new wavelets */
  /***************************************************************************/
  double **muR,**muI;
  double orden[4],a0[4],quanto[4];
  double scale, a,b=0.,corr;
  double maxh=-1e30, minh=1e30;
  double hmsm= 0., dhmsm= 0.;
  int Nev=0;
  int ix,iy,ia,ip,iwav;
  int degree;

  int npoints, nbox;
  
#ifdef _PARSE_FRACTAL_PARAMETERS_
  nbox = p_frac->nbox, npoints = p_frac->npoints;
#else
  nbox = NBOX, npoints = NPOINTS;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
  
  TrackNullAlloc( muR=matrix2D(dimy,dimx) );
  TrackNullAlloc( muI=matrix2D(dimy,dimx) );
  
  copy(dimx,dimy,signal,muR,NULL);
  modderiva(dimx,dimy,muR);
  FFFT2D(dimx,dimy,muR,muI,-1);
  
  maxh=-1e30;  minh=1e30;
  N = wt2d_transform( dimx, dimy, sc0, quanto, orden, ord_der, npoints,
		      muR, muI, expon, &minh, &maxh );
  /*  minh=HMIN; maxh=HMAX;  */
  
  IFVERBOSE 
    WarningVV("Minimal exponent: %f ; maximal exponent: %f\n", minh,maxh);
  
  /*	Computation of the associated histograms of exponents	*/
  for(iy=0;iy<dimy;iy++)   
    for(ix=0;ix<dimx;ix++)	{
      ip = (int)(Nbox * (expon[iy][ix]-minh) / (maxh-minh));
      if(ip<0) ip=0;
      else if(ip>nbox)  ip = nbox;
      histoh[0][ip] += 1.;
      if(msm[iy][ix] != C0)	{
	Nev++;
	histoh[1][ip] += 1.;
      }
    }
  
  for(ip=0;ip<=nbox;ip++)  {
    histoh[0][ip] = histoh[0][ip] / ((double)(dimx*dimy));
    histoh[1][ip] = histoh[1][ip] / ((double)Nev);
    if(histoh[0][ip] > 1.e-30) 
      histoh[2][ip] = histoh[1][ip] / histoh[0][ip];
    else histoh[2][ip] = 0.;
    b += histoh[2][ip];
  }
  
  for(ip=0;ip<=nbox;ip++)    {
    histoh[2][ip] = histoh[2][ip] / b;
    a = minh + (maxh-minh) / ((double)nbox)*((double)ip);
    hmsm +=  a * histoh[2][ip];
    dhmsm += a * a * histoh[2][ip];
  }
  dhmsm = sqrt(dhmsm-hmsm*hmsm);
  IFVERBOSE WarningVV("Effective MSM at: %f +- %f\n",hmsm,dhmsm);
  
  return OK;
} // end of write_expon_density



