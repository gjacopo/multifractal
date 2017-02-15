
/*      main
 * replaced by
 *      create_multifractal
 */

/*      genera_levi( void)
 * replaced by
 *      genera_levi( prob_levi )
 */

/*      genera_multifractal( int leff, double **signal)
 * replaced by
 *      genera_multifractal( leff, wav, ord_der, prob_exp, signal)
 */

/*      genera_multifractal_1D( leff, serie );
 * replaced by
 *      genera1D_multifractal( leff, wav, sc0, ord_der, NULL,serie  );
 */

/*      genera_multifractal_2D( leff, image );
 * replaced by
 *      genera2D_multifractal( leff, wav, ord_der, NULL, image );
 */

/*      genera_alpha( int dimeta, int Nwav, double *alpha);
 * replaced by
 *      genera_alpha(dimeta, Nwav, (prob_exp=NULL), alpha0 );
 */

/*      genera_alpha_1D( dimeta, Nwav, alpha0 );
 * replaced by
 *      genera1D_alpha(dimeta, Nwav, (prob_exp=NULL), alpha0 );
 */

/*      genera_alpha_2D( dimeta, Nwav, alpha0 );
 * replaced by
 *      genera2D_alpha(dimeta, Nwav, (prob_exp=NULL), alpha0 );
 */

/*      genera_h( );
 * unchanged
 */

/*      poisson( lambda );
 * unchanged
 */

 /*      normal_standard( );
 * unchanged
 */

 /*      levi_standard();
 * unchanged
 */

/*      binomial( h0, h1 );
 * unchanged
 */

/*      wavelet_1D( leff, sc, wave );
 * replaced by
 * wavelet1D_define_unit( leff, sc, WAVBASE, DERWAVBASE, 
 *                        D0[WAVBASE], wave );
 */

/*      wavelet_2D( leff, sc, wave );
 * replaced by
 * wavelet2D_define_unit( leff, leff, sc, WAVBASE, DERWAVBASE, 
 *                        D0[WAVBASE], wave );
 */

/*      lee_datos( leff, nombre_in, signal );
 * replaced by
 *      read_data( leff, (D_space==1)?1:leff, nombre_in, signal );
 */

/*      graba_datos( leff, nombre_in, signal );
 * replaced by
 *      write_data( leff, (D_space==1)?1:leff, nombre_in, signal );
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

/*      moda( leff, signal );
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

