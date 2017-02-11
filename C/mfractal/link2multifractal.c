


/*      calcula_multifractal_1D( dimx, verb, signal, expon, mm_h );
 * replaced by
 *      calcula1D_multifractal( dimx, signal, expon, mm_h );
 */

/*      escala_wavelet_1D( dimx )
 * with :
 *  SC0_1D[ORD_DER][WAV]
 *  ORD_DER 
 *  EXP_WL_1D[WAV] 
 * replaced by
 *      wavelet1D_escala( dimx, SC0_1D[ORD_DER][WAV] )
*/

/*      genera_wavelet_1D( dimx, sc, ord_der, expon, wave );
 * replaced by
 *      wavelet1D_genera( dimx, sc, wave );
*/

/*      define_wavelet_1D( dimx, sc, ord_der, expon, wave );
 * with :
 *  expon <=0 => gaussian
 *  expon >0 => lorenzian (with exponent=expon)
 * replaced by
 *      wavelet1D_define( dimx, sc, wav, (ord_der=0), wave );
 * with :
 *  exponent = wav/2
 */

/*      normaliza_wavelet_1D( dimx, sc, ord_der, wave );
 * replaced by
 *      wavelet1D_normaliza( dimx, sc, ord_der, wave );
 */

/*      escala_lineal_2D( dimx, dimy );
 * replaced by
 *      escala2D_lineal( dimx, dimy );
 */

/*      escala_wavelet_2D( dimx, dimy )
 * replaced by
 *      wavelet2D_escala( dimx, dimy, SC0_2D[ORD_DER][WAV] );
 */

/*      escala_wavelet_2D_lineas( dimx, dimy )
 * with :
 *  SC0_2D[ORD_DER][WAV]
 *  THETAU
 * replaced by
 *      wavelet2D_escala_line( dimx, dimy, SC0_2D[ORD_DER][WAV] );
 */

/*      modderiva( dimx, dimy, modg );
 * unchanged
 */

/*       modderiva_1D( dimx, dimy, modg );
 * replaced by
 *       modderiva_line( dimx, dimy, modg, THETAU );
 */

/*       calcula_multifractal_2D( dimx, dimy, verb, signal, expon, mm_h );
 * replaced by
 *       calcula2D_multifractal( dimx, dimy, signal, expon, mm_h );
 */

/*       genera_wavelet_2D( dimx, dimy, sc, wave)
 * replaced by
 *       wavelet2D_genera( dimx, dimy, sc, wave );
 */

/*       define_wavelet_2D( dimx, dimy, sc, wave );
 * with :
 *  EXP_WL_2D[4] = { -1., 1., 1.5 , 2. }; 
 *  WAV = 0 => gaussian
 *  WAV \in [1,3] => lorenzian (with exponent=WAV/2+0.5)
 * replaced by
 *       wavelet2D_define( dimx, dimy, sc, WAV, (ord_der=0), wave );
 */

/*       void genera_wavelet_1D_2D( dimx, dimy, sc, wave );
 * replaced by
 *       wavelet2D_genera_line( dimx, dimy, sc, wave );
*/

/*       define_wavelet_1D_2D( dimx, dimy, sc, wave);
 * with :
 *  EXP_WL_2D[4] = { -1., 1., 1.5 , 2. }; 
 *  WAV = 0 => gaussian
 *  WAV \in [1,3] => lorenzian (with exponent=WAV/2+0.5)
 * replaced by
 *       wavelet2D_define_line( dimx, dimy, sc, THETA0, WAV, (ord_der=0), wave );
 */

/*       normaliza_wavelet_2D( idimx, int dimy, sc, wave );
 * replaced by
 *       wavelet2D_normaliza( dimx, dimy, sc, ORD_DER, wave );
*/

/*       anorma1_lineas( dimx, dimy, data);
 * replaced by
 *        anorma1_line( dimx, dimy, data, THETA0, NULL );
*/

/*       visualiza_color( dimx, dimy, blout, dext, h0, h1, expon );
 * replaced by
 *       sprintf(name,"sing_total.%s.ppm",dext );
 *       visualise_color( dimx, dimy, name, h0, h1, expon );
*/


/*       visualiza_gris_trozos( dimx, dimy, blout, dext, h0, h1, delta_h,
 *                              expon );
 * replaced by
 *       visualise_gris_pieces( dimx, dimy, dext, h0, h1, delta_h, expon );
*/
