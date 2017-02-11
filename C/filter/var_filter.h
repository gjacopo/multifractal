#ifndef   	VAR_OPERATOR_H_
#define   	VAR_OPERATOR_H_

#ifdef __cplusplus
extern "C"
{
#endif	/* __cplusplus */

  /* ...and a default initialisation of it */  
  ParFILT p0_fil = {
    /** Variables related to FFT analysis */
    FLAG_MEMORY,               // 1: flag_memory
    /** Variable related to derivative filters */
    MODE_DERIVA,          // 2: flag_deriva
    MODE_FREQ,          // 3: mode_freq
    /** Variables related to WT analysis */
    WAV,                  // 4: wav
    ORDDER,              // 5: ord_der
    S0,                   // 6: scale0
    THETA0,               // 7: theta
    {1.,  1.,  1., .5},   // 7: D0
    MINSCALE,             //  : minscale
    MAXSCALE,              // : maxscale
    NVOICES,              //  : nvoices
    SCTIME,               //  : sctime
    SCRATIO               //  : scratio
  };

  ParFILT *p_fil=p_fil0; /* default initialisation of structure */

  /* Variables used in WTMM estimation */
  double qDefArray[]={-4.,-3.6,-3.2,-3.,-2.8,-2.6,-2.4,-2.2,-2.,-1.8,-1.6,-1.4,-1.2,-1.1,-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.05,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.3,2.6,3.,3.5,4.,5.,6.,7.,8.};
  int nqDef=65;
  
  /* Minimum scale for each wavelet.in order: haar, morlet, gaussian and lorentzian*/
  const double D0[]= {1.,  1.,  1., .5}; 
 
#ifdef __cplusplus
}
#endif		/* __cplusplus */

#endif 	    /* !VAR_OPERATOR_H_ */
