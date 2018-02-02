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

double dEdtdV_pbh(double Omega_dm_h2, double Mpbh, double fpbh, double z) {
  double rhoPBH; /* density of PBH dark matter today */
  double f_elec = 0.142, f_phot = 0.06; /* relative fraction of injected particles from evaporation */

  rhoPBH = 10537.04*Omega_dm_h2*fpbh; /* energy density in eV/cm^3. see DM annihilation above */

  /* Approximate power radiated from PBH as in [arXiv:1612.07738], where Mpbh is in units of 10^10 g */
  /* units of eV/s/cm^3, it seems */
  return 5.34e-5*(4.*f_elec+2.*f_phot)*rhoPBH*cube((1.+z)/Mpbh);
}

/***********************************************************************************
Total energy *injection* rate per unit volume
Add in your favorite energy injection mechanism here
***********************************************************************************/

double dEdtdV_inj(double z, INJ_PARAMS *params){
  return dEdtdV_DM_ann(z, params);
    // + dEdtdV_pbh(params->odmh2, params->Mpbh, params->fpbh, z);
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

  pbh_inj_energy = dEdtdV_pbh(params->odmh2, params->Mpbh, params->fpbh, z);

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
Fraction of energy deposited in to a channel for log normal mass PBHs
*******************************************************************************/

int pbh_dep_frac_log_norm(BSPLINE *bsp, double *pbh_masses, double z, double *eff_frac,
                          ErrorMsg error_message) {
  double zp1 = 1.+z;

  class_call(array_eval_bicubic_bspline(bsp, &zp1, 1,
                                        pbh_masses, _PBH_MASS_BINS_,
                                        eff_frac,
                                        error_message),
             error_message,
             error_message);

  return _SUCCESS_;
}

/*******************************************************************************
Big F energy parameters for each channel for log normal mass PBHs
*******************************************************************************/

int pbh_F_log_norm(PBH *pbh, INJ_PARAMS *params, double z,
                   double *F_ion, double *F_exc, double *F_heat,
                   ErrorMsg error_message) {
  /* First make sure both mean and width are in log10 space */
  double mean = log10(params->Mpbh);
  double width = params->Wpbh;

  /* Compute lower and upper limits for integral */
  double low_lim = fmax(_LOG_PBH_MASS_MIN_, mean - _SIGMA_INTEGRAL_RANGE_*width);
  double upp_lim = fmin(_LOG_PBH_MASS_MAX_, mean + _SIGMA_INTEGRAL_RANGE_*width);

  double pbh_exps[_PBH_MASS_BINS_];
  double pbh_masses[_PBH_MASS_BINS_];

  double *pbh_dist;
  int i;

  /* Generate PBH masses in log10 space and in linear space */
  for (i = 0; i < _PBH_MASS_BINS_; ++i) {
    pbh_exps[i] = low_lim + (double) i / (_PBH_MASS_BINS_-1) * (upp_lim - low_lim);
    pbh_masses[i] = pow(10., pbh_exps[i]);
  }

  /* Create mass distribution on the heap */
  pbh_dist = pbh_log_normal(pbh_exps, mean, width);

  /* Check distribution truncation */
  double trunc_norm = 1.0;

  if (_LOG_PBH_MASS_MIN_ > mean - _SIGMA_TRUNCATE_RANGE_*width) {
    trunc_norm = trapezoidal_integral(pbh_exps, pbh_dist);
  }

  /* Interpolate deposition fractions */
  double f_ions[_PBH_MASS_BINS_];
  double f_excs[_PBH_MASS_BINS_];
  double f_heats[_PBH_MASS_BINS_];

  class_call(pbh_dep_frac_log_norm(pbh->hion,
                                   pbh_masses,
                                   z,
                                   f_ions,
                                   error_message),
             error_message,
             error_message);

  class_call(pbh_dep_frac_log_norm(pbh->excite,
                                   pbh_masses,
                                   z,
                                   f_excs,
                                   error_message),
             error_message,
             error_message);

  class_call(pbh_dep_frac_log_norm(pbh->heat,
                                   pbh_masses,
                                   z,
                                   f_heats,
                                   error_message),
             error_message,
             error_message);

  double energy_inj;

  /* Compute energy depositions for each mass */
  for (i = 0; i < _PBH_MASS_BINS_; ++i) {
    energy_inj = dEdtdV_pbh(params->odmh2, pbh_masses[i], params->fpbh, z);

    f_ions[i] *= energy_inj * pbh_dist[i]/trunc_norm;
    f_excs[i] *= energy_inj * pbh_dist[i]/trunc_norm;
    f_heats[i] *= energy_inj * pbh_dist[i]/trunc_norm;
  }

  free(pbh_dist);

  /* Perform integrations */
  *F_ion = trapezoidal_integral(pbh_exps, f_ions);
  *F_exc = trapezoidal_integral(pbh_exps, f_excs);
  *F_heat = trapezoidal_integral(pbh_exps, f_heats);

  return _SUCCESS_;
}

/*******************************************************************************
Big F energy parameters for each channel for all PBH mass distributions
*******************************************************************************/

int pbh_F(PBH *pbh, INJ_PARAMS *params, double z,
          double *F_ion, double *F_exc, double *F_heat,
          ErrorMsg error_message) {
  if (params->mass_dist == pbh_none) {
    *F_ion = 0.;
    *F_exc = 0.;
    *F_heat = 0.;
  }
  else if (params->mass_dist == pbh_delta) {
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
