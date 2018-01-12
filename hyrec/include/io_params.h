#ifndef IO_PARAMS_H
#define IO_PARAMS_H

/*** Effective rate tables and associated parameters ***/

#define ALPHA_FILE  "hyrec/data/Alpha_inf.dat"     /* Effective recombination coefficients to 2s and 2p */
#define RR_FILE     "hyrec/data/R_inf.dat"         /* Effective transfer rate R_{2p,2s} */
#define TR_MIN 0.004                         /* Minimum Tr in eV */
#define TR_MAX 0.4                           /* Maximum Tr in eV */
#define NTR    100                           /* Number of Tr values */
#define TM_TR_MIN 0.1                        /* Same thing for ratio Tm/Tr*/
#define TM_TR_MAX 200.0
#define NTM 100             

/*** Tables and parameters for radiative transfer calculation ***/

#define TWOG_FILE "hyrec/data/two_photon_tables.dat" 
#define NSUBLYA  140
#define NSUBLYB  271
#define NVIRT    311
#define NDIFF    80

/* Higher-resolution tables  */
/* #define TWOG_FILE "data/two_photon_tables_hires.dat" */
/* #define NSUBLYA  408 */
/* #define NSUBLYB  1323 */
/* #define NVIRT    1493 */
/* #define NDIFF    300 */
/* #define DLNA    8.47e-5 */

#define L2s1s     8.2206               /* 2s -> 1s two-photon decay rate in s^{-1} (Labzowsky et al 2005) */

/************ SWITCHES FOR RADIATIVE TRANSFER. ALL SWITCHES SET TO 1 ARE THE DEFAULT MODEL  ************/

#define EFFECT_A    1    /* 2s-->1s stimulated two-photon decays and non-thermal absorptions */
#define EFFECT_B    1    /* Sub-Lyman alpha two-photon transitions 3s/3d<--> 1s and 4s/4d<-->1s */
#define EFFECT_C    1    /* Super-Lyman alpha two-photon transitions 3s/3d<--> 1s and 4s/4d<-->1s */
#define EFFECT_D    1    /* Raman scattering from 2s and 3s/3d */
#define DIFFUSION   1    /* Lyman alpha frequency diffusion */

#endif
