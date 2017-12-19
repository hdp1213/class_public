#ifndef ENERGY_INJECTION_H
#define ENERGY_INJECTION_H

#include "hyrectools.h"
/* CLASS common_pbh.h for BSPLINE definition */
#include "common_pbh.h"
/* CLASS common.h for all good CLASS stuff */
#include "common.h"
#include "arrays.h"
/* HyRec pbh.h for PBH definition */
#include "pbh.h"

#define _PBH_MASS_BINS_ 41 /* number of bins used for trapezoidal integration */
#define _LOG_PBH_MASS_MIN_ 5
#define _LOG_PBH_MASS_MAX_ 7
#define _SIGMA_INTEGRAL_RANGE_ 4 /* how many sigmas +/- the mean should integration of the mass distribution take place over */
#define _SIGMA_TRUNCATE_RANGE_ 4 /* how many sigmas +/- the mean should truncation correction start */

/* Structure with all energy injection parameters */
/* If adding a new energy injection process 
   make sure to add relevant parameters here */

typedef struct {

  double odmh2;                 /* Omega_dm h^2 */
  
  double pann, pann_halo;       /* DM annihilation parameter in the smooth background and in haloes */
                                /* Units of pann and pann_halo are cm^3/s/GeV */

  double ann_z, ann_zmax, ann_zmin, ann_var; /* Parameters for the variation of pann(z) */
  double ann_z_halo;                         /* Characteristic redshift for annihilation in haloes */
    
  double Mpbh, fpbh;           /* Mass and fraction of DM made of primordial black holes */
  double Wpbh;                 /* Width of PBH mass distribution, in log10-space */

  enum pbh_mass_distributions mass_dist; /* Mass distribution info */

  int on_the_spot;            /* if set to 1 assume energy deposition rate = injection rate */
                              /* Otherwise solves for deposition given injection with simple recipe */
  
} INJ_PARAMS;

#ifdef __cplusplus
extern "C" {
#endif

double dEdtdV_inj(double z, INJ_PARAMS *params);
void update_dEdtdV_dep(double z_out, double dlna, double xe, double Tgas,
                       double nH, double H, INJ_PARAMS *params, double *dEdtdV_dep);

/* Delta mass distribution functions */
int pbh_dep_frac_delta(BSPLINE *bsp, double Mpbh, double z, double *eff_frac,
                       ErrorMsg error_message);
int pbh_F_delta(PBH *pbh, INJ_PARAMS *params, double z,
                double *F_ion, double *F_exc, double *F_heat,
                ErrorMsg error_message);

/* Log normal mass distribution functions */
int pbh_dep_frac_log_norm(BSPLINE *bsp, double *Mpbh, double z, double *eff_frac,
                          ErrorMsg error_message);
int pbh_F_log_norm(PBH *pbh, INJ_PARAMS *params, double z,
                   double *F_ion, double *F_exc, double *F_heat,
                   ErrorMsg error_message);

/* Wrapper */
int pbh_F(PBH *pbh, INJ_PARAMS *params, double z,
          double *F_ion, double *F_exc, double *F_heat,
          ErrorMsg error_message);

double chi_heat(double xe);
double chi_ion(double xe);
double chi_exc(double xe);

/* Mass distribution methods */
double log10_normal(double x, double mu, double sigma);

double* pbh_log_normal(double* pbh_mass_exponents,
                       double pbh_mass_mean,
                       double pbh_mass_width);

/* Trapezoidal integration routine */
double trapezoidal_integral(double* __restrict__ x, double* __restrict__ integrand);

#ifdef __cplusplus
}
#endif

#endif
