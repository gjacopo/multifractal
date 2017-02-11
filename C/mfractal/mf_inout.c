/* ===================================
** mf_inout.c
** started on Wed Jan 24 19:11:25 2007 
** ===================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <inout.h>  	
#include <io_grafic.h>  	

#include <mfractal.h>
#include <mf_inout.h>


/***************************************************************************/
int Dh_read( char *name_in, double *h, double *Dh, double *errDh ) {
  /***************************************************************************/
  double **series;
  int dims,dimt;
  int it, Nh;
  FILE *canal;
  
  /* first open and read the dimension of the file */
  TrackNull( canal=fopen(name_in,"rt"), "Error opening file in 'r' mode" );
  dims = column_serie_temp(canal);
  dimt = line_serie_temp(canal);

  if((h==NULL) || (Dh==NULL) || (errDh==NULL)) {
    /* We just want to return the dimension */
    if((dims>3) && (dimt>0))  Nh = dimt;
    else                      Nh = ERROR;
    
  } else {
    /* We do read the data */
    TrackNullAlloc( series=matrix2D(dims,dimt) );
    read_serie_temp( canal, dims, dimt, series );
    
    if((dims>3) && (dimt>0))  {
      Nh = dimt;
      for( it=0; it<Nh; it++ )	{
	h[it] = series[0][it];
	Dh[it] = series[1][it];
	errDh[it] = series[2][it];
      }
    }  else   Nh = ERROR;
  }
  
  /*      Memory release and end     */
  fclose(canal);
  if((h!=NULL) && (Dh!=NULL) && (errDh!=NULL))
    free_matrix2D(series,dims);
  
  return Nh;
} // end of load_Dh


/***************************************************************************/
int Dh_write( char *name, int Nr, double *h, double *Dh, double *errDh ) {
  /***************************************************************************/
  FILE *canal;
  int ip;
  
  TrackNull( canal=fopen(name,"wt"), "Error opening file in 'w' mode" );
  for( ip=0; ip<Nr; ip++ )      {
    Dh_th = theoretical_Dh(h[ip]);
    fprintf(canal,"%f  %f  %f  %f\n",h[ip],Dh[ip],errDh[ip],Dh_th);
  }
  
  fclose(canal);

  return OK;
}


/***************************************************************************/
int line_serie_temp( FILE *canal ) {
  /***************************************************************************/
  int pass=TRUE;
  int dimt;
  
  TrackEOF( fseek(canal,0,SEEK_SET), "Error with fseek positioning");
  for( dimt=0; pass; dimt++ ) {
    column_serie_temp(canal);
    if(fgetc(canal) == EOF) pass = FALSE;
  }
  TrackEOF( fseek(canal,0,SEEK_SET), "Error with fseek positioning");

  return dimt;
} // end of line_serie_temp


/***************************************************************************/
int column_serie_temp( FILE *canal ) {
  /***************************************************************************/
  char camp[MAXNAMELENGTH];
  int status, pass=TRUE;
  int dims;
  
  for( dims=0; pass; dims++ ) 
    if((fscanf(canal,"%s",camp)==EOF) ||
       ((status=fgetc(canal))==EOF) || 
       (status==0x0a)) 
      pass = FALSE;       
  
  return dims;
} // end of column_serie_temp


/***************************************************************************/
int read_serie_temp( FILE *canal, int dims, int dimt, double **series ) {
  /***************************************************************************/
  char camp[MAXNAMELENGTH];
  float data;
  int it,is,pass=TRUE;

  for( it=0; (it<dimt)&&pass; it++ )
    for( is=0; (is<dims)&&pass; is++ ) {
      if(fscanf(canal,"%s",camp) == EOF) pass=FALSE;
      IF(pass)	{
	sscanf(campo,"%f",&data);
	series[0][is][it] = (double)data;
      }
    }
  
  return OK;
} // end of read_serie_temp


/***************************************************************************/
int write_geomap ( char *name, double **Q ) {
  /***************************************************************************/
  FILE *canal;
  int ileffs,inums;
  
  TrackNull( canal=fopen(name,"wt"), "Error opening file in 'w' mode" );
  for(inums=0;inums<Nnums;inums++)    {
    for(ileffs=0;ileffs<Nleffs;ileffs++)
      fprintf(canal,"%f  ",Q[inums][ileffs]);
    fprintf(canal,"\n");
  }
  fclose(canal);

  return OK;
} // end of write_geomap


/***************************************************************************/
int write_typemap ( char *name, double **Q ) {
  /***************************************************************************/
  FILE *canal;
  int dim1,dim2;
  int i1= 0,i2;
  
  TrackNull( canal=fopen(name,"wt"), "Error opening file in 'w' mode" );
#ifdef _PARSE_FRACTAL_PARAMETERS_
  if(p_frac->type_mfsim == TYPLOGPOISSON) 
#else
    if(TYPE_MFSIM == TYPLOGPOISSON) 
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
      {
	dim1 = 1,	dim2 = NLPs;
	fprintf(canal,"%f  ",Q[i1][0]);
 	for(i2=1;i2<dim2;i2++)	  {
	  if(codinfs[i2-1 ]>codinfs[i2]) fprintf(canal,"\n");
	  fprintf(canal,"% f  ",Q[i1][i2]);
	} 
	
      }  else	{
	dim1 = Nmeans,  dim2 = Nsigmas;
	for(i1=0;i1<dim1;i1++)	  {
	  for(i2=0;i2<dim2;i2++)
	    fprintf(canal,"%f  ",Q[i1][i2]);
	  fprintf(canal,"\n");
	}
      }
  
  fclose(canal);
  
  return OK;
} // end of write_typemap


/***************************************************************************/
int visualise_gris_pieces( int dimx, int dimy, char *dext,
			   double m, double M, double delta_h, 
			   double **data) {
  /***************************************************************************/
  char name[MAXNAMELENGTH];
  char ***manifold,***excluded;
  double singmax;
  int nn;
  int ix,iy,ip;
  double blout;
  
#ifdef _PARSE_IO_PARAMETERS_
  blout=p_io->block_out;
#else
  blout=BLOCKOUT;
#endif /*!_PARSE_IO_PARAMETERS_*/

  nn = (M-m) / delta_h;

  TrackNullAlloc( manifold=cmatrix3D(nn,dimy,dimx) );
  TrackNullAlloc( excluded=cmatrix3D(2,dimy,dimx) );

  for(iy=0;iy<dimy;iy++)    
    for(ix=0;ix<dimx;ix++)      {
      for(ip=0;ip<nn;ip++) manifold[ip][iy][ix]= (char)255; // C0;
      excluded[0][iy][ix]= (char)255; // C0;
      excluded[1][iy][ix]= (char)255; // C0;
      
      ip=(int) floor((data[iy][ix]-m)/delta_h);
      //      ip=(int) floor(0.5*(data[iy][ix]-m)/delta_h+.5);
      if((ip>=0)&&(ip<nn)) manifold[ip][iy][ix] = (char)0; //  CP;
      if(data[iy][ix]<m-delta_h) excluded[0][iy][ix] = (char)0; //  CP;
      if(data[iy][ix]>M-delta_h) excluded[1][iy][ix] = (char)0; //  CP;
    }
    
  /*	The manifolds are represented	*/  
  for(ip=0;ip<nn;ip++)    {
    sprintf(name,"manifold.%02d.%s.gif",ip,dext); // manifold
    write_mask_block(dimx,dimy,blout,name,manifold[ip]);
  }
  
  sprintf(name,"excluded-.%s.gif",dext); // excluded-
  write_mask_block(dimx,dimy,blout,name,excluded[0]);
  sprintf(name,"excluded+.%s.gif",dext); // excluded+
  write_mask_block(dimx,dimy,blout,name,excluded[1]);
  
  /*     Freeing memory before closign the loop    */
  
  free_cmatrix3D(excluded,2,dimy);
  free_cmatrix3D(manifold,nn,dimy);

  return OK;
} // end of visualise_gris_pieces


/********** A VOIR *********/

/***************************************************************************/
int lecture_image( int dimx, int dimy, int dimxinrim, int dimyinrim, 
		    int dimv, int dimz, int iz, int bd, char *name_in, 
		    char *ext, double ***signal, double *med_cr) {
  /***************************************************************************/
  char name[90];
  double **aux,***auxt;
  int levels;
  int ix,iy,ic;
  
  if(FOTO)  {
    if(XMAX*YMAX>0)	{
      auxt=matrix3D(dimv,dimyinrim,dimxinrim);
      aux=matrix2D(YMAX,XMAX);
      if(dimv==3) levels=read_color_block(dimxinrim,dimyinrim,1,name_in,auxt,med_cr);
      else levels=read_foto_gris(dimxinrim,dimyinrim,name_in,aux);
      for(ic=0;ic<dimv;ic++) {
	paddwindow(dimxinrim,dimyinrim,X0,Y0,XMAX,YMAX,auxt[ic],aux);
	coarseres(dimx,dimy,BLOCK,aux,signal[ic]);
	med_cr[ic]=anorma(dimx,dimy,signal[ic],NULL);
      }
      free_matrix3D(auxt,dimv,dimyinrim);
      free_matrix2D(aux,YMAX);
      
    }    else  {
      if(dimv==3)	{
	for(ic=0;ic<dimv;ic++) fill0(dimx,dimy,signal[ic],NULL);
	levels=read_color_block(dimx,dimy,BLOCK,name_in,signal,med_cr);
      } else if(dimv==1){
	fill0(dimx,dimy,signal[0],NULL);
	aux=matrix2D(dimyinrim,dimxinrim);
	levels=read_foto_gris(dimxinrim,dimyinrim,name_in,aux);
	coarseres(dimxinrim,dimyinrim,BLOCK,aux,signal[0]);
	free_matrix2D(aux,dimyinrim);
	med_cr[0]=anorma(dimx,dimy,signal[0],NULL);
      }
    }
    
  }  else {
    if(XMAX*YMAX>0)      {
      auxt=matrix3D(dimv,YMAX,XMAX);
      levels=read_inrimage_window(dimxinrim,dimyinrim,dimv,dimz,iz,bd,
				  X0,Y0,XMAX,YMAX,name_in,auxt);
      for(ic=0;ic<dimv;ic++)
	coarseres(XMAX,YMAX,BLOCK,auxt[ic],signal[ic]);
      free_matrix3D(auxt,dimv,YMAX);
    }    else	{
      auxt=matrix3D(dimv,dimyinrim,dimxinrim);
      levels=read_inrimage(dimxinrim,dimyinrim,dimv,dimz,iz,bd,name_in,auxt);
      for(ic=0;ic<dimv;ic++)
	coarseres(dimxinrim,dimyinrim,BLOCK,auxt[ic],signal[ic]);
      free_matrix3D(auxt,dimv,dimyinrim);
    }
    for(ic=0;ic<dimv;ic++)
      med_cr[ic]=anorma(dimx,dimy,signal[ic],NULL);
  }
  
  if (ORIGVISU) /* Visualize normalized original signal */ 
    write_foto(VIDEO,BLOUT,dimx,dimy,dimv,dimz,iz,"norma",ext,signal);
  
  return(levels);
} // end of lecture_image


/***************************************************************************/
void write_corte( int dimx, int dimy, int dimv, int xeff, int yeff,  
		  char *ext, int *levels, double *med_cr, double **meang, 
		  char ***msm, double ***gx, double ***gy) {
  /***************************************************************************/
  /* Called by the main              */
  /***************************************************************************/
  
  FILE *canal;
  char name[90];
  int ic;
  
  for(ic=0;ic<dimv;ic++)    {
    if(dimv==1) sprintf(name,"mux.%s.dat",ext);
    else sprintf(name,"mux.%s-V%02d.dat",ext,ic);
    write(xeff,yeff,name,gx[ic]);
    
    if(dimv==1) sprintf(name,"muy.%s.dat",ext);
    else sprintf(name,"muy.%s-V%02d.dat",ext,ic);
    write(xeff,yeff,name,gy[ic]);
    
    if(dimv==1) sprintf(name,"mean.%s.dat",ext);
    else sprintf(name,"mean.%s-V%02d.dat",ext,ic);
    canal=fopen(name,"wb");
    fwrite(levels,sizeof(int),1,canal);
    fwrite(&(med_cr[ic]),sizeof(double),1,canal);
    fwrite(&(meang[ic][0]),sizeof(double),1,canal);
    fwrite(&(meang[ic][1]),sizeof(double),1,canal);
    fclose(canal);
    
    /* It is necessary to rewrite the file 'contour' to fit the correct 
     * size */
    if(dimv==1) sprintf(name,"contour.%s.gif",ext,ic);
    else sprintf(name,"contour.%s-V%02d.gif",ext,ic);
    write_mask_block(dimx,dimy,1.,name,msm[ic]);
  }
  
} // end of write_corte


/***************************************************************************/
void read_corte( int dimx, int dimy, int dimv, int xeff, int yeff, 
		char *ext, int *levels, double *med_cr, double **meang, 
		char ***msm, double ***gx, double ***gy) {
  /***************************************************************************/
  /* Called by the main              */
  /***************************************************************************/

  FILE *canal;
  char name[90];
  int ic;
  
  for(ic=0;ic<dimv;ic++) {
    
    if(dimv==1) sprintf(name,"mux.%s.dat",ext);
    else sprintf(name,"mux.%s-V%02d.dat",ext,ic);
    read(xeff,yeff,name,gx[ic]);
    
    if(dimv==1) sprintf(name,"muy.%s.dat",ext);
    else sprintf(name,"muy.%s-V%02d.dat",ext,ic);
    read(xeff,yeff,name,gy[ic]);
    
    if(dimv==1) sprintf(name,"mean.%s.dat",ext);
    else sprintf(name,"mean.%s-V%02d.dat",ext,ic);
    canal=fopen(name,"rb");
    fread(levels,sizeof(int),1,canal);
    fread(&(med_cr[ic]),sizeof(double),1,canal);
    fread(&(meang[ic][0]),sizeof(double),1,canal);
    fread(&(meang[ic][1]),sizeof(double),1,canal);
    fclose(canal);
    
    if(dimv==1) sprintf(name,"contour.%s.gif",ext);
    else sprintf(name,"contour.%s-V%02d.gif",ext,ic);
    read_gif_rgb(dimx,dimy,0,name,msm[ic],msm[ic],msm[ic]);
    
  }
} // end of read_corte

/***************************************************************************/
int write_expon_density( double *prob, double *mm, int N,
			 char *ext, char *name_in ) {
  /***************************************************************************/
  double h;
  int ip;
  FILE *canal; 
  
  int nbox =
#ifdef _PARSE_FRACTAL_PARAMETERS_
    p_frac->nbox;
#else
  NBOX;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/
  
  TrackNullAlloc( canal=fopen(name,"wt") );
  
  for( ip=0; ip<=nbox; ip++ )      {
    h = mm[0] + ((double)ip) / nbox * (mm[1]-mm[0]);
    prob[ip] = prob[ip] * nbox / ((double)N*(mm[1]-mm[0]));
    fprintf(canal,"%f	%f\n",h,prob[ip]);
  }
  
  fclose(canal);
  
  return OK;
} // end of write_expon_density


/***************************************************************************/
void phase_color( double mx, double my, char *Red, char *Green, char *Blue) {
  /***************************************************************************/
  /* Called by prepare_sources_fase_block 
   * Color representation of the space of phases:
   *                    
   *                   \pi/2
   *               \     |     /
   *                \    |    /
   *                 \  BLUE /
   *                  \  |  /
   *                   \ | /
   *    \pi--- GREEN  ---0---  RED --- 0
   *                   / | \
   *                  /  |  \ 
   *                 /       \
   *                / YELLOW  \
   *               /           \
   * 
   */
  /***************************************************************************/
  
  double theta;
  
  if((fabs(mx)>1e-30)||(fabs(my)>1e-30))    {
    theta=anguphase_colorlo(mx,my);
    if((theta<=M_PI/4.)||(theta>=7.*M_PI/4.))   {
      *Red=(char)255;
      *Green=(char)0;
      *Blue=(char)0;
    }    else if(theta<=3.*M_PI/4.)      {
      *Red=(char)0;
      *Green=(char)0;
      *Blue=(char)255;
    }    else if(theta<=5.*M_PI/4)      {
      *Red=(char)0;
      *Green=(char)255;
      *Blue=(char)0;
    }    else	{
      *Red=(char)255;
      *Green=(char)255;
      *Blue=(char)0;
    }
  }  else    {
    *Red=(char)0;
    *Green=(char)0;
    *Blue=(char)0;
  }
} // end of phase_color


/***************************************************************************/
void write_multifractal( int dimx, int dimy, int dimv, int dimz, 
		    int iz, /* int xeff, int yeff, */char*ext, 
		    char ***msm, double ***signal ) {
/***************************************************************************/
  
  double ***expon;
  FILE *canal;
  char name[90];
  int ic;
  
  if( (expon=matrix3D( dimv, dimy, dimx )) == NULL)  
    exit(-1);
  
  if(dimz>1) sprintf(name, "histo_mf_%s-Z%d.fi", ext, iz);
  else sprintf(name, "histo_mf_%s.fi", ext );
  canal=fopen(name,"wt");
  for(ic=0;ic<dimv;ic++)
    multimf_histo(dimx,dimy, canal,
		       msm[ic],signal[ic], expon[ic]); 
  fclose(canal);
  
  write_foto(VIDEO,BLOUT,dimx,dimy,dimv,dimz,iz,"sing_fractal",ext,expon);
  
  free_matrix3D(expon,dimv, dimy);
} // end of write_multifractal


/***************************************************************************/
void vuelcaresult( int dimx, int dimy,  double sc0, double qs,
		   double minimo, double maximo, int N, 
		   const char *name_wv, double **expon, char *ext) {
/***************************************************************************/
  double per;
  
  per=100.*((double)N)/((double)dimx*dimy);
  WarningV("Analysis obtained for %s wavelet\n",name_wv);
  WarningVV("minimum scale %f scale quantum %f\n\n",sc0,qs);
  Warning("Results:\n\n");
  WarningVV("Minimal exponent %f; Maximal exponent %f:\n\n",
	 minimo,maximo);
  WarningV("Percentage of good points:  %f %%\n",per);
} // end of vuelcaresult


/***************************************************************************/
void process_sources( int dimx, int dimy, int xeff, int yeff, 
		      double **gx, double **gy) {
  /***************************************************************************/

  double **mod;
  int ix,iy;
  
  mod=matrix2D(yeff,xeff);
  
  for(iy=0;iy<dimy;iy++)  {
    for(ix=0;ix<dimx;ix++)	{
      mod[iy][ix]=sqrt(gx[iy][ix]*gx[iy][ix]+gy[iy][ix]*gy[iy][ix]);
      if(mod[iy][ix]>1e-30) mod[iy][ix]=log(mod[iy][ix]);
      else mod[iy][ix]=-30.*log(10.);
    }
  }
  filtro(xeff,yeff,0.,2.,mod);
  write_foto_4(xeff,yeff,"prueba.gif",mod);
  
  free_matrix2D(mod,yeff);
} // end of process_sources


/***************************************************************************/
void write_histogram( char *name, double minimo, double maximo, 
		       double *prob ) {
  /***************************************************************************/

  FILE *canal;
  double x;
  int ip;
  
  canal=fopen(name,"wt");
  for(ip=0;ip<=NBOX;ip++)   {
    x=minimo+((double)ip)/NBOX*(maximo-minimo);
    fprintf(canal,"%f   %f\n",x,prob[ip]);
  }
  fclose(canal);
  
} // end of write_histogram


/***************************************************************************/
void represent_sources( int dimx, int dimy, int dimv, int dimz, int iz, char *ext, 
			 double ***gx, double ***gy) {
  /***************************************************************************/
  /* Represents the vectorial field of sources by its norm and its orientation.
   * =========================================================================
   * Called by compute_sources   */ 
  /***************************************************************************/

  char ***Red,***Green,***Blue;
  double ***mod;
  int ic;
  
  /* Calcul et representation du mode */
  mod=matrix3D(dimv,dimy,dimx);
  
  for(ic=0;ic<dimv;ic++)
    prepare_sources_mod(dimx,dimy,1,gx[ic],gy[ic],mod[ic]);
  write_foto(VIDEO,BLOUT,dimx,dimy,dimv,dimz,iz,"sources",ext,mod);
  free_matrix3D(mod,dimv,dimy);
  
  /* Calcul et representation de l'orientation */
  Red=cmatrix3D(dimv,dimy,dimx);
  Green=cmatrix3D(dimv,dimy,dimx);
  Blue=cmatrix3D(dimv,dimy,dimx);
  
  for(ic=0;ic<dimv;ic++)
    prepare_sources_fase_block(dimx,dimy,1,gx[ic],gy[ic],
			       Red[ic],Green[ic],Blue[ic]);
  write_foto_RGB(VIDEO,BLOUT,dimx,dimy,dimv,dimz,iz,"sources_vec",ext,Red,Green,Blue);
  
  /* Liberation des memoires */
  free_cmatrix3D(Red,dimv,dimy);
  free_cmatrix3D(Green,dimv,dimy);
  free_cmatrix3D(Blue,dimv,dimy);
} // end of represent_sources


/***************************************************************************/
void represent_sources_mod( int dimx, int dimy, char *ext, int silog,
			     double **mux, double **muy) {
  /***************************************************************************/
  /* Represents the norm of the vectorial field of sources
   * =========================================================================
   * Called by the main */
  /***************************************************************************/

  char name[90];
  double **aux;
  double buffx,buffy;
  double maximo,minimo;
  int ix,iy;

  aux=matrix2D(dimy,dimx);

  maximo=-1e30;
  minimo=1e30;
  for(ix=0;ix<dimx;ix++)   
    for(iy=0;iy<dimy;iy++)   {
      buffx=mux[iy][ix];
      buffy=muy[iy][ix];
      aux[iy][ix]=sqrt(buffx*buffx+buffy*buffy);
      if(silog) {
	if(aux[iy][ix]>1e-30) {
	  maximo=fMax(maximo,log(aux[iy][ix]));
	  minimo=fMin(minimo,log(aux[iy][ix]));
	}
      }
    }
  
  if(silog)   
    for(ix=0;ix<dimx;ix++)	
      for(iy=0;iy<dimy;iy++)	{
	if(aux[iy][ix]>1e-30)
	  aux[iy][ix]=log(aux[iy][ix])-minimo;
	else aux[iy][ix]=0.;
      }
  
  strcpy(name,"sources.");
  strcat(name,ext);
  strcat(name,".gif");
  write_foto_block(dimx,dimy,BLOUT,name,aux);

  free_matrix2D(aux,dimy);
} // end of represent_sources_mod


/***************************************************************************/
void prepare_sources_mod( int dimx, int dimy, int silog, 
			  double **mux, double **muy, double **mod ) {
  /***************************************************************************/
  /* Computes the norm of a vectorial field (eventually log-representtion). 
   * Used for the representtion of the sources.
   * =========================================================================
   * Called by represent_sources   */
  /***************************************************************************/
  
  double buffx,buffy;
  double maximo,minimo;
  int ix,iy;
  
  for(iy=0;iy<dimy;iy++)  
    for(ix=0;ix<dimx;ix++)	{
      mod[iy][ix]=sqrt(mux[iy][ix]*mux[iy][ix]+muy[iy][ix]*muy[iy][ix]);
      if(silog)	    {
	if(mod[iy][ix]>1e-30)
	  mod[iy][ix]=log(mod[iy][ix]);
	else mod[iy][ix]=-30*log(10.);
      }
    }
} // end of prepare_sources_mod


/***************************************************************************/
void represent_sources_fase( int dimx, int dimy, char *ext,
			      double **mux, double **muy) {
  /***************************************************************************/
  /* Represents the orientation of the vectorial field of sources
   * =========================================================================
   * Called by  the main */
  /***************************************************************************/

  char name[90];
  char **Red,**Green,**Blue;
  int dimxf,dimyf;

  dimxf=BLOUT*dimx;
  dimyf=BLOUT*dimy;
  
  Red=cmatrix2D(dimyf,dimxf);
  Green=cmatrix2D(dimyf,dimxf);
  Blue=cmatrix2D(dimyf,dimxf);
  
  prepare_sources_fase_block(dimx,dimy,BLOUT,mux,muy,Red,Green,Blue);
  
  strcpy(name,"sources_vec.");
  strcat(name,ext);
  strcat(name,".gif");
  
  write_RGB_block(dimxf,dimyf,1.,name,Red,Green,Blue);
  
  free_cmatrix2D(Red,dimyf);
  free_cmatrix2D(Green,dimyf);
  free_cmatrix2D(Blue,dimyf);
} // end of represent_sources_fase


/***************************************************************************/
void prepare_sources_fase_block( int dimx, int dimy, double block, double **mux, 
				 double **muy, char **Red, char **Green, char **Blue) {
  /***************************************************************************/
  /* Computes the orientation of a vectorial field. Used for the representtion 
   * of the sources.
   * =========================================================================
   * Called by  the represent_sources_fase */
  /***************************************************************************/
  
  char aR,aG,aB;
  double mx,my;
  int beff;
  int xmax,ymax;
  int ix,iy;
  int ibx,iby;

  if(block>=1.)    {
    beff=(int) block;
    for(ix=0;ix<dimx;ix++)      
      for(iy=0;iy<dimy;iy++)	{
	phase_color(mux[iy][ix],muy[iy][ix],&aR,&aG,&aB);
	for(iby=0;iby<beff;iby++)	  
	  for(ibx=0;ibx<beff;ibx++)	    {
	    Red[iy][ix]=aR;
	    Green[iy][ix]=aG;
	    Blue[iy][ix]=aB;
	  }
      }
    
  }  else    {
    beff=(int)(1./block);
    xmax=dimx/beff;
    ymax=dimy/beff;
    
    for(iy=0;iy<ymax;iy++)      
      for(ix=0;ix<xmax;ix++)	{
	mx= my= 0.;
	for(iby=0;iby<beff;iby++)    
	  for(ibx=0;ibx<beff;ibx++)   {
	    mx+=mux[iy][ix];
	    my+=muy[iy][ix];
	  }
	mx=mx/((double)(beff*beff));
	my=my/((double)(beff*beff));
	phase_color(mux[iy][ix],muy[iy][ix],&(Red[iy][ix]),&(Green[iy][ix]),
		    &(Blue[iy][ix])); 
	/* ou:
	phase_color(mx,my,&(Red[iy][ix]),&(Green[iy][ix]), &(Blue[iy][ix])); 
	* ??? */
      }
    
  }
} // end of prepare_sources_fase_block

/***************************************************************************/
int write_multifractal_histo( int dimx, int dimy, /* int xeff, int yeff, */
			      FILE *canal, char **msm, double **signal, 
			      double **expon ) {
  /***************************************************************************/
  
  int nbox;

#ifdef _PARSE_FRACTAL_PARAMETERS_
  nbox = p_frac->nbox;
#else
  nbox = NBOX;
#endif /*!_PARSE_FRACTAL_PARAMETERS_*/

  if(canal != NULL) {
    for(ip=0;ip<=nbox;ip++)    {
      a = minh + (maxh-minh) / ((double)nbox)*((double)ip);
      fprintf(canal,"%f  %f  %f  %f\n",a,histoh[0][ip],
	      histoh[1][ip],histoh[2][ip]);
    }
    fprintf(canal,"\n");
  }


} // end of write_multifractal_histo


/***************************************************************************/
void save_sources( int dimx, int dimy, int dimv, int dimz, int iz, char *ext, 
			 double ***gx, double ***gy) {
  /***************************************************************************/
  /* Saves the vectorial field of sources.
   * =========================================================================
   * Called by the main program.   */ 
  /***************************************************************************/

  char name[90];
  int ic;
  
  write_foto(VIDEO,BLOUT,dimx,dimy,dimv,dimz,iz,"sources_x",ext,gx);
  write_foto(VIDEO,BLOUT,dimx,dimy,dimv,dimz,iz,"sources_y",ext,gy);
} // end of save_sources


/***************************************************************************/
void represent_sources_angulo( int dimx, int dimy, int dimv, int dimz, int iz, char *ext, 
			      double ***gx, double ***gy) {
  /***************************************************************************/
  /* Represents the vectorial field of sources by its norm and its orientation.
   * =========================================================================
   * Called by compute_sources   */ 
  /***************************************************************************/
  
  char name[90];
  double ***mod, ***ang;
  double mm[2];
  int ic;
  
  mod=matrix3D(dimv,dimy,dimx);
  
  for(ic=0;ic<dimv;ic++)
    prepare_sources_mod(dimx,dimy,1,gx[ic],gy[ic],mod[ic]);
  write_foto(VIDEO,BLOUT,dimx,dimy,dimv,dimz,iz,"sources",ext,mod);
  free_matrix3D(mod,dimv,dimy);
  
  ang=matrix3D(dimv,dimy,dimx);
  for(ic=0;ic<dimv;ic++) {
    prepare_sources_angulo( dimx, dimy, 1, gx[ic], gy[ic], ang[ic] );
    /* extrema(dimx, dimy,ang[ic] ,mm, NULL); */
  }
  write_foto(VIDEO,BLOUT,dimx,dimy,dimv,dimz,iz,"sources_vec",ext,ang);
  free_matrix3D(ang,dimv,dimy);
  
} // end of represent_sources_angulo


/***************************************************************************/
void prepare_sources_angulo( int dimx, int dimy, double block, 
			     double **mux, double **muy, double **ang ) {
  /***************************************************************************/
  
  double mx,my;
  int beff;
  int xmax,ymax;
  int ix,iy;
  int ibx,iby;
  double a;
  
  if(block>=1.)    {
    beff=(int) block;
    for(ix=0;ix<dimx;ix++)      
      for(iy=0;iy<dimy;iy++)	{
	a=angulo(mux[iy][ix],muy[iy][ix]);
	for(iby=0;iby<beff;iby++)	  
	  for(ibx=0;ibx<beff;ibx++)	 
	    ang[iy][ix]=a;
      }
    
  }  else    {
    beff=(int)(1./block);
    xmax=dimx/beff;
    ymax=dimy/beff;
    
    for(iy=0;iy<ymax;iy++)      
      for(ix=0;ix<xmax;ix++)	{
	mx= my= 0.;
	for(iby=0;iby<beff;iby++)    
	  for(ibx=0;ibx<beff;ibx++)   {
	    mx+=mux[iy][ix];
	    my+=muy[iy][ix];
	  }
	mx=mx/((double)(beff*beff));
	my=my/((double)(beff*beff));
	ang[iy][ix]=angulo(mux[iy][ix],muy[iy][ix]);
	/* ou 
	   ang[iy][ix]=angulo(mx,my);	   
	   * ??? */
      }
    
  } 
} // end of prepare_sources_angulo


/***************************************************************************
void visualize_color( int dimx, int dimy, char *dext, double **expon, 
		      const char *name_wv) {
  char name[90];
  char **Red,**Green,**Blue;
  double h0,h1,rep;
  int ix,iy;
  
  h0=-0.5;
  h1=1.;
  
  Red=cmatrix2D(dimy,dimx);
  Green=cmatrix2D(dimy,dimx);
  Blue=cmatrix2D(dimy,dimx);
  
  for(iy=0;iy<dimy;iy++)  {
    for(ix=0;ix<dimx;ix++)	{
      rep=(expon[iy][ix]-h0)/(h1-h0);
      if(rep<0.) rep=0.;
      if(rep>1.) rep=1.;
      paleta(rep,&(Red[iy][ix]),&(Green[iy][ix]),&(Blue[iy][ix]));
    }
  }
  sprintf(name,"sing_total.%s.%s.ppm",dext,name_wv);
  write_ppm(dimx,dimy,1,name,Red,Green,Blue);
  
  free_cmatrix2D(Red,dimy);
  free_cmatrix2D(Green,dimy);
  free_cmatrix2D(Blue,dimy);
}
***************************************************************************/

/***************************************************************************
void visualize_gris( int dimx, int dimy, char *dext,double **expon, 
		     const char *name_wv) {
  char name[90];
  double **aux;
  double h0,h1;
  int ix,iy;

  h0=-0.5;
  h1=.5;
  
  aux=matrix2D(dimy,dimx);
  
  for(iy=0;iy<dimy;iy++)	  {
    for(ix=0;ix<dimx;ix++)      {
      aux[iy][ix]=-expon[iy][ix];
      if(aux[iy][ix]>-h0) aux[iy][ix]=-h0;
      if(aux[iy][ix]<-h1) aux[iy][ix]=-h1;
    }
  }
  sprintf(name,"sing_total.%s.%s.gif",dext,name_wv);
  write_foto_4(dimx,dimy,name,aux);
  
  free_matrix2D(aux,dimy);
}
***************************************************************************/



/***************************************************************************
double visualize_gris_pieces( int dimx, int dimy, double h_inf, double delta_h, 
			      char *dext,
			      double **expon, const char *name_wv) {
 
  char name[90];
  char ***manifold,***excluded;
  double singmax;
  int ix,iy,ip;
  
  singmax=h_inf+NSING*2*delta_h;
  fprintf(stderr,"First represented singularity: %f\n",h_inf);
  fprintf(stderr,"First excluded singularity: %f\n",singmax);
  
  manifold=cmatrix3D(NSING,dimy,dimx);
  excluded=cmatrix3D(2,dimy,dimx);
  
  for(iy=0;iy<dimy;iy++)    {
    for(ix=0;ix<dimx;ix++)        {
      for(ip=0;ip<NSING;ip++) manifold[ip][iy][ix]=C0;
      excluded[0][iy][ix]=C0;
      excluded[1][iy][ix]=C0;
      
      ip=(int) floor(0.5*(expon[iy][ix]-h_inf)/delta_h+.5);
      if((ip>=0)&&(ip<NSING)) manifold[ip][iy][ix]=CP;
      if(expon[iy][ix]<h_inf-delta_h) excluded[0][iy][ix]=CP;
      if(expon[iy][ix]>singmax-delta_h) excluded[1][iy][ix]=CP;
    }
  }
  
  //      The manifolds are represented           
  
  for(ip=0;ip<NSING;ip++)    {
    sprintf(name,"manifold_%02d.%s.%s.gif",ip,dext,name_wv);
    write_binary_block(dimx,dimy,BLOUT,name,manifold[ip]); 
  }
  
  sprintf(name,"excluded-.%s.%s.gif",dext,name_wv);
  write_binary_block(dimx,dimy,BLOUT,name,excluded[0]);
  sprintf(name,"excluded+.%s.%s.gif",dext,name_wv);
  write_binary_block(dimx,dimy,BLOUT,name,excluded[1]);
  
//    Freeing memory before closign the loop    
  
  free_cmatrix3D(excluded,2,dimy);
  free_cmatrix3D(manifold,NSING,dimy);

  return singmax;
}
  ***************************************************************************/
 
