/* ===================================
** filter.h
** started on Fri Jan 19 17:49:29 2007 
** ===================================
*/

#ifndef   	_FILTER_H_
#define   	_FILTER_H_

/* Default values for derivative related flags */
#ifndef FORWARD
#define FORWARD 1
#endif
#ifndef HALFPIXEL
#define HALFPIXEL 0
#endif
#ifndef MODE_DERIVA
#define MODE_DERIVA HALFPIXEL
#endif

#ifndef FREQSIN
#define FREQSIN 1
#endif
#ifndef FREQID
#define FREQID 0
#endif
#ifndef MODE_FREQ
#define MODE_FREQ FREQSIN
#endif

#ifndef THETA0
#define THETA0 0.
#endif

/* Default values for FFT related flags */
#ifndef FLAG_MEMORY
#define FLAG_MEMORY FALSE
#endif

/* Default values for WT related flags */
#ifndef WAVMORLET
#define WAVMORLET -2
#endif
#ifndef WAVGAUSS
#define WAVGAUSS -1
#endif
#ifndef WAVHAAR
#define WAVHAAR 0
#endif
#ifndef SWAVLORENTZ 
#define SWAVLORENTZ 1
#endif
#ifndef WAV 
#define WAV WAVGAUSS
#endif
#ifndef ORDDER
#define ORDDER 2
#endif

#ifndef S0
#define S0 0.33 
#endif

#ifndef MINSCALE
#define MINSCALE 1. 
#endif
#ifndef MAXSCALE
#define MAXSCALE (MINSCALE)
#endif
#ifndef SCTIME
#define SCTIME 5.       
#endif
#ifndef NVOICES
#define NVOICES 10
#endif

#ifndef SCRATIO
#define SCRATIO (1./16.)
#endif
#ifndef SCSTEP
#define SCSTEP 2. //Scale step
#endif

#ifndef WAVRANGE
#define WAVRANGE 64.
#endif

#endif 	    /* !_FILTER_H_ */

