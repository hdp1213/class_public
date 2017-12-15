#ifndef HISTORY_H
#define HISTORY_H
/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         history.h: functions for recombination history                                        */
/*                                                                                               */
/*************************************************************************************************/

#include "hyrectools.h"
#include "hydrogen.h"
#include "helium.h"
#include "energy_injection.h"
/* pbh.h defines PBH struct */
#include "pbh.h"
/* CLASS common.h for all good CLASS stuff */
#include "common.h"

/* Structure for HyRec internal parameters */ 

typedef struct {
  double T0;                                /* CMB temperature today in K*/
  double obh2, omh2, odeh2, okh2, orh2;     /* density parameters */ 
  double Nnueff;                            /* effective number of neutrinos */
  double fHe;                               /* Helium fraction by number */
  double nH0;                               /* density of hydrogen today in cm^{-3} (changed from m^{-3}) */ 
    
  double fsR, meR;              /* fine-structure constant alpha/alpha(today) 
                                    and me/me(today) (Added April 2012)*/
  
  INJ_PARAMS *inj_params;     /* Structure containing all Energy-injection parameters */

} REC_COSMOPARAMS;

#ifdef __cplusplus
extern "C" {
#endif

double rec_HubbleRate(REC_COSMOPARAMS *param, double z);

double rec_Tmss(double z, double xe, REC_COSMOPARAMS *cosmo, double dEdtdV_dm, double F_heat);

double rec_dTmdlna(double z, double xe, double Tm, REC_COSMOPARAMS *cosmo, double dEdtdV_dm, double F_heat);

int rec_get_xe_next1_He(REC_COSMOPARAMS *cosmo, double z_in, double *xHeII,
                        double dxHeIIdlna_prev[2], int *post_saha,
                        ErrorMsg error_message);

int rec_xH1_stiff(int model, REC_COSMOPARAMS *cosmo, double z, double xHeII, double *xH1,
                  HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz_rad,
                  double dEdtdV_dm, double F_ion, double F_exc, int *stiff, long int Nz, ErrorMsg error_message);

int get_rec_next2_HHe(int model, REC_COSMOPARAMS *cosmo, double z_in, double Tm,
                      double *xH1, double *xHeII, HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz_rad,
                      double dxHIIdlna_prev[2], double dxHeIIdlna_prev[2], double dEdtdV_dm,
                      double F_ion, double F_exc, int *stiff, long int NZ, ErrorMsg error_message);

int rec_get_xe_next1_H(int model, REC_COSMOPARAMS *cosmo, double z_in, double xe_in, double Tm_in,
                       double *xe_out, double *Tm_out, HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz_rad,
                       double dxedlna_prev[2], double dEdtdV_dm,
                       double F_ion, double F_exc, double F_heat, int *stiff, long int Nz, ErrorMsg error_message);

int rec_get_xe_next2_HTm(int model, REC_COSMOPARAMS *cosmo,
                         double z_in, double xe_in, double Tm_in, double *xe_out, double *Tm_out,
                         HYREC_ATOMIC *atomic, RADIATION *rad, unsigned iz_rad,
                         double dxedlna_prev[2], double dTmdlna_prev[2], double dEdtdV_dm,
                         double F_ion, double F_exc, double F_heat, long int Nz, ErrorMsg error_message);

int rec_build_history(int model, double zstart, double zend,
                      REC_COSMOPARAMS *cosmo, HYREC_ATOMIC *atomic,
                      RADIATION *rad, PBH *pbh, double *xe_output, double *Tm_output,
                      ErrorMsg error_message);

#ifdef __cplusplus
}
#endif

typedef struct{
  HYREC_ATOMIC *atomic;
  REC_COSMOPARAMS *cosmo;
  double zmax;
  double zmin;
  long int Nz;
  double *xe_output;
  double *Tm_output;
  RADIATION *rad; 
  PBH *pbh;
} HYREC_DATA;

#endif
