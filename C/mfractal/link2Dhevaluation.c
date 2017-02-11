
int METHOD=15; // Composite variable, defining the method/s to be used 
              // in the analysis
              // 1*GH+2*GMWP+4*MOMENTS+8*WTMM
              // WTMM is called via LastWave, if LASTWAVE flag is set

int main(int argc, char *argv[]) {

  /*	DATA			*/
  char path[90],base[90];
  
  /*	QUANTITIES TO BE COMPUTED	*/
  double **Qg,**Qw,**Qm,**Ql;
  
  /*	AUXILIAR VARIABLES	*/
  Quality Qtot;
  char nombre[90];
  char GHname[90],GMWPname[90],Momname[90],WTMMname[90];
  int ileffs,inums;
  int iLPs,imeans,isigmas;
  
  /*		Program		*/
  parse_arguments(argc,argv);
  
  /*        Re-fitting parameters      */
  
  /*
    To avoid problems, the parameter DERIVA_MODO is shortcircuited  
    to a plain rightwards finite difference
  */
  DERIVA_MODO=1;
  
  /*     LEFF is converted to the minimum power of 2 larger than it   */
  LEFF=dimensiona(LEFF);
  
  /*        Inputting the path to the data      */
  if(VERBOSE) printf("Please, give me the path to the data files\n");
  scanf("%s",path);
  
  /*       Generating the base for the MF names     */
  if(VERBOSE)
    printf("Processing for %d output %dD signals\n",NSERIES,D_space);
  
  /*       Analyzing series with the methods      */  
  if(GEO_MAP+TYPE_MAP==0) analiza_series(path);
  
  if(GEO_MAP)	{
    
    /*        Initialization    */
    genera_base(0,0,base);
    Qg=reservar_matriz(Nnums,Nleffs);
    sprintf(GHname,"GeoMap-GH-%s.dat",base);
    /*        Obtaining the quality matrices       */
    for(inums=0;inums<Nnums;inums++)	    {
      NSERIES=nums[inums];
      for(ileffs=0;ileffs<Nleffs;ileffs++)	{
	LEFF=leffs[ileffs];
	Qtot=analiza_series(path);
	Qg[inums][ileffs]=Qtot.g;
	graba_geomapa(GHname,Qg);
      }      
    }
    
    /*      Memory release before end    */
    liberar_matriz(Qg,Nnums);
  }
  
  if(TYPE_MAP)    {
    NSERIES=1;
    LEFF=65536;
    switch(TYPE)      {
    case 0:
      /*        Memory reservation for the quality matrices    */      
      Qg=reservar_matriz(1,NLPs);
      /*        Obtaining the quality matrices       */      
      for(iLPs=0;iLPs<NLPs;iLPs++)	{
	HINF=hinfs[iLPs];
	CODINF=codinfs[iLPs];
	Qtot=analiza_series(path);
	Qg[0][iLPs]=Qtot.g;
      }
      /*       Recording the results         */      
      sprintf(base,"Log-Poisson");
      if(D_space==1) strcat(base,"_1D");
      else strcat(base,"_2D");
      sprintf(nombre,"TypeMap-GH-%s.dat",base);
      graba_typemap(nombre,Qg);
      /*      Memory release before end    */      
      liberar_matriz(Qg,1);
      break;
      
    case 1:
      /*        Memory reservation for the quality matrices    */
      Qg=reservar_matriz(Nmeans,Nsigmas);
      /*        Obtaining the quality matrices       */
      for(imeans=0;imeans<Nmeans;imeans++)	{
	MU=means[imeans];
	for(isigmas=0;isigmas<Nsigmas;isigmas++)			{
	  SIGMA=sigmas[isigmas];
	  Qtot=analiza_series(path);
	  Qg[imeans][isigmas]=Qtot.g;
	}
      }      
      /*       Recording the results         */      
      sprintf(base,"Log-Normal");
      if(D_space==1) strcat(base,"_1D");
      else strcat(base,"_2D");
      sprintf(nombre,"TypeMap-GH-%s.dat",base);
      graba_typemap(nombre,Qg);
      /*      Memory release before end    */
      liberar_matriz(Qg,Nmeans);
      break;
      
    default:
      printf("Map not implemented for this type, sorry!\n");
      break;
    }
    
  }

  exit(0);
}


/*      genera_base( verbose, leff, base );
 * replaced by
 *      genera_base( leff, base );
 */

/*      Quality analiza_series( char *path)
 * replaced by
 *      analiza_multifractal
 */

/*      media_exp( dimx, sc0, expon );
 * replaced by
 *      media_expon( dimx, (D_space==1)?1:dimx, sc0, expon, NULL );
 */

/*      estima_Dh_g( path, base, &Nh_g, &h_g, &Dh_g, &errDh_g );
 * replaced by
 *      int nbox=
 *      #ifdef _PARSE_FRACTAL_PARAMETERS_
 *          p_frac->nbox;
 *      #else
 *        NBOX;
 *      #endif // !_PARSE_FRACTAL_PARAMETERS_
 *      // we need to allocate memory with more than necessary 
 *      h=(double*)calloc(nbox,sizeof(double));
 *      Dh=(double*)calloc(nbox,sizeof(double));
 *      errDh=(double*)calloc(nbox,sizeof(double));
 *      sprintf(name,"%s%s",path,base);
 *      Dh_estima_g( name, leff, shift_g, h_g, Dh_g, errDh_g );
 */

/*      estima_Dh_w( path, base, &Nh_w, &h_w, &Dh_w, &errDh_w );
 * replaced by
 *      int nbox=
 *      #ifdef _PARSE_FRACTAL_PARAMETERS_
 *          p_frac->nbox;
 *      #else
 *        NBOX;
 *      #endif // !_PARSE_FRACTAL_PARAMETERS_
 *      // we need to allocate memory with more than necessary 
 *      h=(double*)calloc(nbox,sizeof(double));
 *      Dh=(double*)calloc(nbox,sizeof(double));
 *      errDh=(double*)calloc(nbox,sizeof(double));
 *      sprintf(name,"%s%s",path,base);
 *      Dh_estima_gmwp( name, leff, &shift_w, &shift_g,
 *                      h_w, Dh_w, errDh_w );
 */

/*      estima_Dh_mom( path, base, Nh_m, &h_m, &Dh_m, &errDh_m );
 * replaced by
 *      N_m=nmoms-1;
 *      h_m=(double*)calloc(N_m,sizeof(double));
 *      Dh_m=(double*)calloc(N_m,sizeof(double));
 *      errDh_m=(double*)calloc(N_m,sizeof(double));
 *      sprintf(name,"%s%s",path,base);
 *      Dh_estima_moments( name, dimx, Moms, nmoms, Dist, ndists, 
 *      		   h_m, Dh_m, errDh_m );
 */

/*      estima_Dh_wtmm( path, base, Nh_wtmm, &h_wtmm, &Dh_wtmm, &errDh_wtmm );
 * replaced by the sequence:
 *      h_wtmm=(double*)calloc(nmoms,sizeof(double));
 *      Dh_wtmm=(double*)calloc(nmoms,sizeof(double));
 *      errDh_wtmm=(double*)calloc(nmoms,sizeof(double));
 *      Nh_wtmm=nmoms;
 *      sprintf(name,"%s%s",path,base);
 *      Dh_estima_wtmm( name, dimx, Moms, Nmom, h_wtmm, Dh_wtmm, errDh_wtmm );
 */

/*      WTMM_PF_compute( signal, dimx, wav, order, nq, qArray, 
 *                            nsc, scArray, ExtWTlis, Z, n_ext );
 * replaced by
 *      wtmm1D_pf( signal, dimx, wav, order, nq, qArray, nsc, scArray,
 *                 ExtWTlis, Z, n_ext );
 */

/*      WTMM_find_extrema( wtrans, wtrans_ind, dimx, sc );
 * replaced by
 *      wtmm1D_find_extrema( wtrans, wtrans_ind, dimx, sc );
 */

/*      PF_extrema_track( wtrans, wtrans_ind, maxsig, maxsig_ind, jj, n_max,
 *                    	  nq, qArray, nsc, scArray, ExtWTlis, Z );
 * replaced by
 *      pf1D_extrema_track ( wtrans, wtrans_ind, maxsig, maxsig_ind, jj, n_max,
 *                    	  nq, qArray, nsc, scArray, ExtWTlis, Z );
 */

/*      WT_transform( signal, dimx, wav, order, sc, wtrans );
 * replaced by
 *      wt1D_transform( signal, dimx, wav, order, sc, wtrans );
 */

/*      W_projection( signal, filt, dimx, sc, wtrans );
 * replaced by
 *      wt1D_projection( signal, filt, dimx, sc, wtrans );
 */

/*      directSpecCompute( Z, dimx, nq, qArray, nws, wsArray, Tauq, H, Dh );
 * replaced by
 *      spectrum1D_direct( Z, dimx, nq, qArray, nws, wsArray, Tauq, H, Dh );
 */

/*      canonMeanCompute( ExtWTlis, dimx, n_ext, nq, qArray, nws, wsArray,
 *		          sTq, sTqLogT ) 
 * replaced by
 *      mean1D_canon( ExtWTlis, dimx, n_ext, nq, qArray, nws, wsArray,
 *		      sTq, sTqLogT )
 */

/*     canonSpecCompute( sTq, sTqLogT, dimx, nq, qArray, nws, wsArray,
 *		         Tauq, H, Dh );
 * replaced by
 *      spectrum1D_canon( sTq, sTqLogT, dimx, nq, qArray, nws, wsArray,
 *		          Tauq, H, Dh );
 */

/*      acumula_momentos( signal, moments );
 * replaced by
 *      moments_accumula( leff, ((dim_space==DIM1D)?1:leff), signal, 
*                         moms, nmoms, dist, ndists, moments ); 
 */

/*      calcula_taup( moments, taup );
 * replaced by
 *      moments_taup( moments, nmoms, dist, ndists, taup );
 */

/*      N_m=simple_Legendre_transform( taup, &h_m, &Dh_m, &errDh_m );
 * replaced by the sequence:
 *      N_m=Nmom-1;
 *      h_m=(double*)calloc(N_m,sizeof(double));
 *      Dh_m=(double*)calloc(N_m,sizeof(double));
 *      errDh_m=(double*)calloc(N_m,sizeof(double));
 *      legendre_transform( moms, N_m, taup, h_m, Dh_m, errDh_m );
 */

/*      genera_expon( leff, signal, expon );
 * replaced by
 *      expon_genera( leff, signal, expon ); 
 */

/*      genera_expon_1D( leff, serie, expon );
 * replaced by
 *      expon1D_genera( leff, serie, expon );
 */

/*      genera_expon_2D( leff, image, expon );
 * replaced by
 *      expon2D_genera( leff, image, expon );
 */

/*      acumula_histograma( leff, mh, expon, histo );
 * replaced by
 *      histogram_accumula( leff, ((dim_space==DIM1D)?1:leff), mh, expon, histo );
 */

/*      calcula_Dh_histo( sc0, mh, histo, h, Dh, errDh );
 * replaced by the sequence:
 *        int nbox=
 *      #ifdef _PARSE_FRACTAL_PARAMETERS_
 *          p_frac->nbox;
 *      #else
 *        NBOX;
 *      #endif // !_PARSE_FRACTAL_PARAMETERS_
 *      histo_r=(double*)calloc(nbox,sizeof(double));
 *      h_r=(double*)calloc(nbox,sizeof(double));
 *      Nh = Dh_filter( histo, mh, histo_r, h_r );
 *      // The number of exponents is now known; 
 *      // we allocate memory accordingly 
 *      (*h)=(double*)calloc(Nh,sizeof(double));
 *      (*Dh)=(double*)calloc(Nh,sizeof(double));
 *      (*errDh)=(double*)calloc(Nh,sizeof(double));
 *      // create a mute variable 
 *      width=(int*)calloc(nbox,sizeof(int)) ;
 *      ifill1D(nbox,1,width,NULL); 
 *      Dh_errorbar_weight(sc0, Nh, histo_r, h_r, width, h, Dh, errDh );
 *      free(histo_r);
 *      free(h_r);
 *      free(width); 
 * See also function Dh_histo
 */

/*      calcula_Dh_histo_no_ponderado( sc0, mh, histo, h, Dh, errDh );
 * replaced by the sequence:
 *        int nbox=
 *      #ifdef _PARSE_FRACTAL_PARAMETERS_
 *          p_frac->nbox;
 *      #else
 *        NBOX;
 *      #endif // !_PARSE_FRACTAL_PARAMETERS_
 *      histo_r=(double*)calloc(nbox,sizeof(double));
 *      h_r=(double*)calloc(nbox,sizeof(double));
 *      width=(int*)calloc(nbox,sizeof(int));
 *      Nh = Dh_filter_weight( histo, mh, histo_r, h_r, width );
 *      // The number of exponents is now known; 
 *      // we allocate memory accordingly 
 *      (*h)=(double*)calloc(Nh,sizeof(double));
 *      (*Dh)=(double*)calloc(Nh,sizeof(double));
 *      (*errDh)=(double*)calloc(Nh,sizeof(double));
 *      Dh_errorbar_weight(sc0, Nh, histo_r, h_r, width, h, Dh, errDh );
 *      free(histo_r);
 *      free(h_r);
 *      free(width); 
 * See also function Dh_histo_weight
 */

/*      int carga_Dh( char *name_in, int *Nh, h, 
 *                    double **double **Dh, double **errDh );
 * replaced by
 *      int load_Dh( char *name_in, double *h, double *Dh, double *errDh ) ;
 * with Nh returned
 */

/*      lee_serie_temp( name_in, &dims, &dimt, &series );
 * replaced by the sequence
 *      canal=fopen(name_in,"rt");
 *      dims = columnas_serie_temp(canal);
 *      dimt = lineas_serie_temp(canal);
 *      TrackNullAlloc( series=matrix2D(dims,dimt) );
 *      read_serie_temp( canal, dims, dimt, series );
 *      fclose(canal);
 */

/*      columnas_serie_temp( canal );
 * replaced by
 *      column_serie_temp( canal );
 */

/*      lineas_serie_temp( canal );
 * replaced by
 *      
 */

/*      registra_Dh( graba, Nr, nombre, sc0, h, Dh, errDh );
 * replaced by the sequence:
 *      if(graba)  {
 *         canal=fopen(nombre,"wt");
 *         for(ip=0;ip<Nr;ip++)      {
 *           Dh_th=theoretical_Dh(h[ip]);
 *           fprintf(canal,"%f  %f  %f  %f\n",h[ip],Dh[ip],errDh[ip],Dh_th);
 *         }
 *         fclose(canal);
 *      }
 *      registra_Dh( Nr, sc0, h, Dh, errDh );
 */

/*      theoretical_Dh( double h)
 * unchanged            */

/*      theoretical_Deltah( double *hmin, double *hmax)
 * unchanged
 */

/*      lee_datos( leff, name_in, signal );
 * replaced by
 *      read_data( leff, (D_space==1)?1:leff, name_in, signal );
 */

/*      graba_datos( leff, nombre_in, signal );
 * replaced by
 *      write_data( leff, , (D_space==1)?1:leffnombre_in, signal );
 */

/*      lee_datos_float( leff, nombre_in, signal );
 * replaced by
 *      read_data_infloat( leff, (D_space==1)?1:leff, nombre_in, signal );
 */ 

/*      graba_datos_float( leff, nombre_in, signal );
 * replaced by
 *      write_data_infloat( leff, (D_space==1)?1:leff, nombre_in, signal );
 */

/*      graba_serie( leff, nombre_in, datos );
 * replaced by
 *      write_serie( leff, nombre_in, datos );
 */

/*      graba_imagen( leff, nombre_in, datos );
 * replaced by
 *      write_image_foto( leff, leff, nombre_in, datos );{
 */

/*      graba_geomapa( nombre, Q );
 * replaced by: 
 *      write_geomap( name, Q );
 */

/*      graba_typemap( nombre, Q );
 * replaced by                                        
 *      write_typemap( name, Q );
 */
 
/*      double moda( leff, signal );
 * replaced by
 *      moda_redef( leff, signal );
 */

/*      moda_1D( dimx, datos );
 * replaced by
 *      moda1D( dimx, datos );
 */

/*      moda_2D( dimx, dimy, datos );
 * replaced by
 *      moda( dimx, dimy, datos );
 */

/*      moda_por_histo( Nhisto, mm, histo );
 * replaced by
 *      moda1D_histo( Nhisto, mm, histo );
 */
