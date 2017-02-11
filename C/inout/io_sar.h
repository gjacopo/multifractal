#ifndef SAR_H
#define SAR_H


#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

//#define XMAXSAR 5632
//#define X0SAR 192

#define X0SAR 384
#define XMAXSAR 5536
#define YMAXSAR 26624
#define Y0SAR 1


void read_sar( int dimx, int dimy, int x0, int y0, char *nombre, double ***data);
void write_sar( int dimx, int dimy, int x0, int y0, char *nombre, double ***data);

#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif
