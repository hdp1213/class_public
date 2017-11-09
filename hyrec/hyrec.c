/*************************************************************************/
/*                    HYREC MAIN FUNCTION                                */
/*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "include/hyrec.h"


int main(void) {
  HYREC_DATA rec_data, rec_data0, rec_data1, rec_data2, rec_data3;
  double zmax = 8000.;
  double zmin = 0.;
  
  hyrec_allocate(&rec_data, zmax, zmin);
  // hyrec_allocate(&rec_data0, zmax, zmin);
  // hyrec_allocate(&rec_data1, zmax, zmin);
  // hyrec_allocate(&rec_data2, zmax, zmin);
  // hyrec_allocate(&rec_data3, zmax, zmin);


  const double pann_SI_to_cgs = 1.7826673e-21;
  rec_data.cosmo->inj_params->on_the_spot = 1;


  double h       = 0.6725;
  double T0      = 2.7255;
  double Omega_b = 0.022252/h/h;
  double Omega_m = 0.14276/h/h;
  double Omega_k = 0.;
  double YHe     = 0.246667;
  double Nnueff  = 3.046;

  double alphaR     = 1.;
  double meR        = 1.;
  double pann       = 1.e-5*pann_SI_to_cgs;
  double pann_halo  = 0.;
  double ann_z      = 0.;
  double ann_zmax   = 0.;
  double ann_zmin   = 0.;
  double ann_var    = 0.;
  double ann_z_halo = 0.;
  double Mpbh       = 0.;
  double fpbh       = 0.;


  hyrec_compute(&rec_data, FULL,
		h, T0, Omega_b, Omega_m, Omega_k, YHe, Nnueff,
		  alphaR, meR, pann, pann_halo, ann_z, ann_zmax,
		  ann_zmin, ann_var, ann_z_halo, Mpbh, fpbh);

#ifdef PBH_DEBUG
  hyrec_free(&rec_data);
  return 0;
#endif

  /*
  hyrec_compute(&rec_data0, FULL,
		h, T0, Omega_b, Omega_m, Omega_k, YHe, Nnueff,
		alphaR, meR, pann, pann_halo, ann_z, ann_zmax,
		ann_zmin, ann_var, ann_z_halo, 1, 1.);

  hyrec_compute(&rec_data1, FULL,
		h, T0, Omega_b, Omega_m, Omega_k, YHe, Nnueff,
		alphaR, meR, pann, pann_halo, ann_z, ann_zmax,
		ann_zmin, ann_var, ann_z_halo, 10., 1.);

  hyrec_compute(&rec_data2, FULL,
		h, T0, Omega_b, Omega_m, Omega_k, YHe, Nnueff,
		alphaR, meR, pann, pann_halo, ann_z, ann_zmax,
		ann_zmin, ann_var, ann_z_halo, 100., 1.);

  hyrec_compute(&rec_data3, FULL,
		h, T0, Omega_b, Omega_m, Omega_k, YHe, Nnueff,
		alphaR, meR, pann, pann_halo, ann_z, ann_zmax,
		ann_zmin, ann_var, ann_z_halo, 1000., 1.);
  */

  /* "FULL" refers to the default, most accurate recombination model
      Other options are (not to use for any precision cosmology):
      EMLA2p2s: effective 4-level atom accounting exactly for a virtually infinite tower of excited states,
                but with a simple analytic treatment of radiative transfer
      RECFAST: Peebles' effective 3-level atom model with a fudge factor of 1.14
      PEEBLES: Peebles' effective 3-level atom model */
  
  double z = zmax;

  while (z > zmin){
    // printf("%f %1.10E %1.10E %1.10E %1.10E %1.10E\n", z, hyrec_xe(z, &rec_data), hyrec_xe(z, &rec_data0),
	   // hyrec_xe(z, &rec_data1), hyrec_xe(z, &rec_data2), hyrec_xe(z, &rec_data3));
    printf("%7.2lf %17.15lf %15.13lf %15.13lf\n", z, hyrec_xe(z, &rec_data), hyrec_Tm(z, &rec_data), rec_data.cosmo->T0*(1.+z));
    z -= 1.;
  }
   
  hyrec_free(&rec_data);
  // hyrec_free(&rec_data0);
  // hyrec_free(&rec_data1);
  // hyrec_free(&rec_data2);
  // hyrec_free(&rec_data3);

  return 0;

}
