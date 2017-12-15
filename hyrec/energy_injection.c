/******************************************************************************************************/
/*                           HYREC: Hydrogen and Helium Recombination Code                            */
/*                      Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                      */
/*                                                                                                    */
/*     energy_injection.c: functions for the energy injection rate by various physical processes      */
/*                                                                                                    */
/*     Written October 2016                                                                          */
/*                                                                                                    */
/******************************************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "include/energy_injection.h"


/***************************************************************************************
Total volumic rate of energy *injection*, in eV/cm^3/s due to DM annihilation 
in the smooth background and in haloes as in Giesen et al 1209.0247
****************************************************************************************/

double dEdtdV_DM_ann(double z, INJ_PARAMS *params){

  double pann_tot, u_min, erfc;
  double zp1, zp1_ann, zp1_max, zp1_min, var, zp1_halo;                        
   
  var       = params->ann_var;
  zp1       = z + 1.; 
  zp1_ann   = params->ann_z + 1.;
  zp1_max   = params->ann_zmax + 1.; 
  zp1_halo  = params->ann_z_halo + 1.; 
  zp1_min   = params->ann_zmin + 1.; 
  
  pann_tot = 0.;
  
  /* Dark matter annihilation in the smooth background */
  if (params->pann > 0.) {

    /* Parametrized variation of pann */
    if (zp1 > zp1_max) pann_tot = params->pann *exp(-var *square(log(zp1_ann/zp1_max)));
    else if (zp1 > zp1_min) {
      pann_tot = params->pann *exp(var*(-square(log(zp1_ann/zp1_max))
				+square(log(zp1/zp1_max))));
    }
    else {
      pann_tot = params->pann *exp(var*(-square(log(zp1_ann/zp1_max))
				+square(log(zp1_min/zp1_max))));
    } 
  }
  
  /* Dark matter annihilation in haloes */
  if (params->pann_halo > 0.) {
    u_min = zp1/zp1_halo;
    erfc  = pow(1.+0.278393*u_min+0.230389*u_min*u_min+0.000972*u_min*u_min*u_min+0.078108*u_min*u_min*u_min*u_min,-4); 
    pann_tot += params->pann_halo *erfc;
  }
  
  /* Make sure that params->odmh2 only contains the WIMP DM, and no PBHs */
  return square(10537.04*params->odmh2*(1.-params->fpbh)) * square(zp1*zp1*zp1) *1e-9* pann_tot;
  /* the prefactor is 3 H100^2/(8 Pi G) c^2 in eV/cm^3, H100 = 100km/s/Mpc */
  /* pann is given in cm^3/s/GeV, multiply by 1e-9 to get cm^3/s/eV */
  
}

/***************************************************************************************
Effect of evaporating primordial black holes
Throughout, Mpbh is in 10^10 g, Teff in Kelvins
Note that if Mpbh is given in 10^10 g, it always cancels with a factor of 10^10 g and so
in this region of the code the parameter behaves as though it is dimensionless.
***************************************************************************************/

/* Rate of energy *injection* per unit volume (in eV/cm^3/s) due to PBHs */

double dEdtdV_pbh(double z, INJ_PARAMS *params) {
  double rhoPBH; /* density of PBH dark matter today */
  double f_elec = 0.142, f_phot = 0.06; /* relative fraction of injected particles from evaporation */

  rhoPBH = 10537.04*params->odmh2*params->fpbh; /* energy density in eV/cm^3. see DM annihilation above */

  /* Approximate power radiated from PBH as in [arXiv:1612.07738], where Mpbh is in units of 10^10 g */
  /* units of eV/s/cm^3, it seems */
  return (params->fpbh > 0.) ? 5.34e-5*(4.*f_elec+2.*f_phot)*rhoPBH*cube((1.+z)/params->Mpbh) : 0.;
}

/***********************************************************************************
Total energy *injection* rate per unit volume
Add in your favorite energy injection mechanism here
***********************************************************************************/

double dEdtdV_inj(double z, INJ_PARAMS *params){
  return dEdtdV_DM_ann(z, params);
    // + dEdtdV_pbh(z, params);
}


/**********************************************************************************
Energy *deposition* rate per unit volume
Essentially use a very simple ODE solution
**********************************************************************************/

void update_dEdtdV_dep(double z_out, double dlna, double xe, double Tgas,
		       double nH, double H, INJ_PARAMS *params, double *dEdtdV_dep) {

  double inj  = dEdtdV_inj(z_out, params);
  
  if (params->on_the_spot == 1) *dEdtdV_dep = inj;

  // Else put in your favorite recipe to translate from injection to deposition
  // Here I assume injected photon spectrum Compton cools at rate dE/dt = - 0.1 n_h c sigma_T E
  // This is valid for E_photon ~ MeV or so.
  
  else { // 0.1 c sigma_T = 2e-15 (cgs). 
    *dEdtdV_dep = (*dEdtdV_dep *exp(-7.*dlna) + 2e-15* dlna*nH/H *inj)
                 /(1.+ 2e-15 *dlna*nH/H);                              
  }       
}

/*******************************************************************************
Fraction of energy deposited in to a channel for delta mass PBHs
*******************************************************************************/

int pbh_dep_frac_delta(BSPLINE *bsp, double Mpbh, double z, double *eff_frac,
                       ErrorMsg error_message) {
  double zp1 = 1.+z;

  class_call(array_eval_bicubic_bspline(bsp, &zp1, 1,
                                        &Mpbh, 1,
                                        eff_frac,
                                        error_message),
             error_message,
             error_message);

  return _SUCCESS_;
}

/*******************************************************************************
Big F energy parameters for each channel for delta mass PBHs
*******************************************************************************/

int pbh_F_delta(PBH *pbh, INJ_PARAMS *params, double z,
                double *F_ion, double *F_exc, double *F_heat,
                ErrorMsg error_message) {
  double pbh_inj_energy;
  double f_ion, f_exc, f_heat;

  pbh_inj_energy = dEdtdV_pbh(z, params);

  class_call(pbh_dep_frac_delta(pbh->hion,
                                params->Mpbh,
                                z,
                                &f_ion,
                                error_message),
             error_message,
             error_message);

  class_call(pbh_dep_frac_delta(pbh->excite,
                                params->Mpbh,
                                z,
                                &f_exc,
                                error_message),
             error_message,
             error_message);

  class_call(pbh_dep_frac_delta(pbh->heat,
                                params->Mpbh,
                                z,
                                &f_heat,
                                error_message),
             error_message,
             error_message);

  *F_ion = pbh_inj_energy * f_ion;
  *F_exc = pbh_inj_energy * f_exc;
  *F_heat = pbh_inj_energy * f_heat;

  return _SUCCESS_;
}

/*******************************************************************************
Big F energy parameters for each channel for log normal mass PBHs
*******************************************************************************/

int pbh_F_log_norm(PBH *pbh, INJ_PARAMS *params, double z,
                   double *F_ion, double *F_exc, double *F_heat,
                   ErrorMsg error_message) {
  double *pbh_masses;

  /* TODO: implement this function */

  return _SUCCESS_;
}

/*******************************************************************************
Big F energy parameters for each channel for all PBH mass distributions
*******************************************************************************/

int pbh_F(PBH *pbh, INJ_PARAMS *params, double z,
          double *F_ion, double *F_exc, double *F_heat,
          ErrorMsg error_message) {
  if (params->mass_dist == pbh_delta) {
    class_call(pbh_F_delta(pbh, params, z, F_ion, F_exc, F_heat, error_message),
               error_message,
               error_message);
  }
  else if (params->mass_dist == pbh_log_norm) {
    class_call(pbh_F_log_norm(pbh, params, z, F_ion, F_exc, F_heat, error_message),
               error_message,
               error_message);
  }
  else {
    // Already know this value must be valid, no need to check here
  }

  return _SUCCESS_;
}


/*******************************************************************************
Fraction of energy deposited in the form of heat, ionization and excitations
*******************************************************************************/

double chi_heat(double xe) { 
  return (1.+2.*xe)/3.; // Approximate value of Chen & Kamionkowski 2004 

  // fit by Vivian Poulin of columns 1 and 2 in Table V of Galli et al. 2013
  // double chi_heat;

  // if (xe < 1.) {
  //   chi_heat = 0.996857*(1.-pow(1.-pow(xe,0.300134),1.51035));
  //   return (chi_heat > 1.) ? 1. : chi_heat;
  // }
  // else
  //   return 1.;
}

double chi_ion(double xe) {
  return (1.-xe)/3.; // Approximate value of Chen & Kamionkowski 2004 

  // fit by Vivian Poulin of columns 1 and 4 in Table V of Galli et al. 2013
  // return 0.369202*pow(1.-pow(xe, 0.463929), 1.70237);
}

double chi_exc(double xe) {
  return 1. - chi_ion(xe) - chi_heat(xe);
}


/*******************************************************************************
Log normal distribution generation methods
*******************************************************************************/

double log10_normal(double x, double mu, double sigma) {
  return 1./(sigma * _SQRT2_ * _SQRT_PI_) * exp(-0.5*pow((x - mu)/sigma, 2.));
}

double* pbh_log_normal(double* pbh_mass_exponents,
                       double pbh_mass_mean,
                       double pbh_mass_width) {
  double* mass_log_normal;
  int i;

  mass_log_normal = malloc(_PBH_MASS_BINS_*sizeof(double));

  for (i = 0; i < _PBH_MASS_BINS_; ++i) {
    mass_log_normal[i] = log10_normal(pbh_mass_exponents[i],
                                      pbh_mass_mean,
                                      pbh_mass_width);
  }

  return mass_log_normal;
}

/*******************************************************************************
Trapezoidal integration function, stolen shamelessly from CLASS and repurposed
to only need to be called from one function one time by HyRec.

Simultaneously compute weights and integral.
All arrays have length _PBH_MASS_BINS_.
*******************************************************************************/

double
trapezoidal_integral(double* __restrict__ x,
                     double* __restrict__ integrand)
{
  int i;
  double w_trapz[_PBH_MASS_BINS_];
  double res = 0.0;

  // Set edge weights:
  w_trapz[0] = 0.5 * (x[1] - x[0]);
  w_trapz[_PBH_MASS_BINS_-1] = 0.5 * (x[_PBH_MASS_BINS_-1] - x[_PBH_MASS_BINS_-2]);

  // Set inner weights:
  for (i = 1; i < (_PBH_MASS_BINS_-1); i++) {
    w_trapz[i] = 0.5 * (x[i+1] - x[i-1]);
    res += integrand[i]*w_trapz[i];
  }

  // Compute remaining integral at edge points
  res += integrand[0]*w_trapz[0];
  res += integrand[_PBH_MASS_BINS_-1]*w_trapz[_PBH_MASS_BINS_-1];

  return res;
}
