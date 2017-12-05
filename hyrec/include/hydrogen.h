#ifndef HYDROGEN_H
#define HYDROGEN_H
/********************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                        */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                                     */
/*                                                                                                      */
/*         hydrogen.h: all functions related to Hydrogen recombination                                  */
/*                                                                                                      */
/*         Units used: cgs + eV (all temperatures in eV)                                                */
/*                                                                                                      */
/*         Version: December 2014                                                                       */
/*                                                                                                      */
/*         Revision history:                                                                            */
/*            - Written November 2010.                                                                  */
/*            - January 2011: updated value of 2s--1s decay rate,                                       */
/*                            changed temperature range for effective rates                             */
/*            - May 2012:   - Using the photon distortion instead of absolute value of radiation field  */
/*                          - Accounting for explicit dependence on alpha and m_e                       */
/*                          - Some definitions moved to header file history.h                           */  
/*            - December 2014: - Accounts for additional energy injection                               */                 
/********************************************************************************************************/ 

#include "io_params.h"
/* CLASS common.h for all good CLASS stuff */
#include "common.h"

/* Definition of different recombination models  */

#define PEEBLES   0    /* Peebles's effective three-level atom */
#define RECFAST   1    /* Effective three-level atom for hydrogen with fudge factor F = 1.14 */
#define EMLA2s2p  2    /* Correct EMLA model, with standard decay rates from 2s and 2p only (accounts for nmax = infinity, l-resolved) */
#define FULL      3    /* All radiative transfer effects included. Additional switches in header file hydrogen.h */

/* When the probability of being ionized from n=2 becomes lower than PION_MAX, 
   switch off radiative transfer calculation as it becomes irrelevant */
#define PION_MAX  1e-2      


/****** CONSTANTS IN CGS + EV UNIT SYSTEM *******/

#define EI   13.598286071938324              /* Hydrogen ionization energy in eV, reduced mass, no relativistic corrections */

/* Energy differences between excited levels of hydrogen -- used often */
#define E21  10.198714553953742
#define E31  12.087365397278509
#define E41  12.748393192442178
#define E32  1.8886508433247664
#define E42  2.5496786384884356

#define hPc       1.239841874331e-04   /* hc in eV cm */
#define mH        0.93878299831e9      /* Hydrogen atom mass in eV/c^2 */ 
#define kBoltz    8.617343e-5          /* Boltzmann constant in eV/K */


#ifdef __cplusplus
extern "C" {
#endif

/*********** EFFECTIVE 3-LEVEL A LA PEEBLES ***************/ 
double SAHA_FACT(double fsR, double meR);
double LYA_FACT(double fsR, double meR);
double L2s_rescaled(double fsR, double meR);
void rescale_T(double *T, double fsR, double meR);

double alphaB_PPB(double TM, double fsR, double meR);
double rec_TLA_dxHIIdlna(double xe, double xHII, double nH, double H, double TM, double TR,
                         double Fudge, double fsR, double meR, double dEdtdV_dm, double dEdtdV_pbh,
                         double f_ion, double f_exc);

#ifdef __cplusplus
}
#endif


/************* EFFECTIVE MULTI LEVEL ATOM *******************/

#define DLNA     8.49e-5    /* Timestep. Maximum compatible with these tables is 8.49e-5 */

/**** Structure containing all atomic data for hydrogen ****/

typedef struct {
  /* Tables of effective rates */
  double *logTR_tab;
  double *TM_TR_tab;
  double **logAlpha_tab[2];
  double *logR2p2s_tab;
  double DlogTR, DTM_TR;

  /* Tables of 2-photon rates */
  double *Eb_tab;       /* Energies of the virtual levels in eV */
  double *A1s_tab;      /* 3*A2p1s*phi(E)*DE */ 
  double *A2s_tab;      /* dLambda_2s/dE * DeltaE if E < Elya dK2s/dE * Delta E if E > Elya */
  double *A3s3d_tab;    /* (dLambda_3s/dE + 5*dLambda_3d/dE) * Delta E for E < ELyb, Raman scattering rate for E > ELyb */
  double *A4s4d_tab;    /* (dLambda_4s/dE + 5*dLambda_4d/dE) * Delta E */
  
} HYREC_ATOMIC;


/**** Structure containing all radiative transfer tables ****/

typedef struct {

  double z0;               // first redshift at which radiation fields are stored
  double **Dfminus_hist;
  double **Dfnu_hist;
  double **Dfminus_Ly_hist;

} RADIATION;

#ifdef __cplusplus
extern "C" {
#endif

void allocate_radiation(RADIATION *rad, long int Nz);
void free_radiation(RADIATION *rad);


void allocate_atomic(HYREC_ATOMIC *atomic);
void normalise_atomic(HYREC_ATOMIC *atomic);
int allocate_and_read_atomic(HYREC_ATOMIC *atomic, ErrorMsg error_message);

void free_atomic(HYREC_ATOMIC *atomic);
int interpolate_rates(double Alpha[2], double DAlpha[2], double Beta[2], double *R2p2s,
                      double TR, double TM_TR, HYREC_ATOMIC *atomic, double fsR, double meR,
                      ErrorMsg error_message);
int rec_HMLA_dxHIIdlna(double xe, double xHII, double nH, double H, double TM, double TR,
                       HYREC_ATOMIC *atomic, double fsR, double meR, double dEdtdV_dm, double dEdtdV_pbh,
                       double f_ion, double f_exc, double *dxHIIdlna, ErrorMsg error_message);
void populate_Diffusion(double *Aup, double *Adn, double *A2p_up, double *A2p_dn, 
                        double TM, double Eb_tab[NVIRT], double A1s_tab[NVIRT]);
int populateTS_2photon(double Trr[2][2], double *Trv[2], double *Tvr[2], double *Tvv[3],
                       double sr[2], double sv[NVIRT], double Dtau[NVIRT],
                       double xe, double xHII, double TM, double TR, double nH, double H, HYREC_ATOMIC *atomic,
                       double Dfplus[NVIRT], double Dfplus_Ly[],
                       double Alpha[2], double DAlpha[2], double Beta[2], double fsR, double meR, double dEdtdV_dm,
                       double dEdtdV_pbh, double f_exc, ErrorMsg error_message);
void solveTXeqB(double *diag, double *updiag, double *dndiag, double *X, double *B, unsigned N);
void solve_real_virt(double xr[2], double xv[NVIRT], double Trr[2][2], double *Trv[2], double *Tvr[2], 
                     double *Tvv[3], double sr[2], double sv[NVIRT]);
int interp_Dfnu(double lna_start, double dlna, double *ytab, unsigned int iz, double lna,
                double* Dfnu_interp, ErrorMsg error_message);
int fplus_from_fminus(double Dfplus[NVIRT], double Dfplus_Ly[], double **Dfminus_hist, double **Dfminus_Ly_hist,
                      double TR, double zstart, unsigned iz, double z, double Eb_tab[NVIRT], ErrorMsg error_message);
int rec_dxHIIdlna(int model, double xe, double xHII, double nH, double H, double TM, double TR,
                  HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz, double z,
                  double fsR, double meR, double dEdtdV_dm, double dEdtdV_pbh, double f_ion, double f_exc,
                  double *result, long int Nz, ErrorMsg error_message);

#ifdef __cplusplus
}
#endif

#endif
