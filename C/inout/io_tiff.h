/* ===================================
** io_tiff.h
** started on Fri Feb 23 11:02:47 2007 
** ===================================
*/

#ifndef   	_IO_TIFF_H_
#define   	_IO_TIFF_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

  ERROR_TYPE read_image_data(FILE *fp, IMAGE *im, int pc);

IMAGE *read_image(char *fn)


#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !_IO_TIFF_H_ */

