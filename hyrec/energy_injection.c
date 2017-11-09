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
  return   dEdtdV_DM_ann(z, params);
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
Fraction of energy deposited in to a channel for PBHs
*******************************************************************************/

double dEdtdV_fraction_pbh(BSPLINE *bsp, double Mpbh, double z) {
  double zp1 = 1.+z;
  double eff_frac;
  int ecode;

  ecode = array_eval_bicubic_bspline(bsp, &zp1, 1,
                                     &Mpbh, 1,
                                     &eff_frac);

  if (ecode != 0) {
    fprintf(stderr, "Error in dEdtdV_fraction_pbh: class_interp2d failed\n");
    exit(ecode);
  }

  return eff_frac;
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


