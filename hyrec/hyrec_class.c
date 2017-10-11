/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         hyrec_class.c: CLASS-specific functions                                               */
/*************************************************************************************************/

#include <math.h> // pow()

#include "hyrec_class.h"

double heat_channel(double xe) {
  double chi_heat;
  //return (1.+2.*xe)/3.; // old approximation from Chen and Kamionkowski

  // coefficient as revised by Galli et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 2 in Table V of Galli et al. 2013)
  if (xe < 1.) {
    chi_heat = 0.996857*(1. - pow(1. - pow(xe, 0.300134), 1.51035));
    return (chi_heat > 1.) ? 1. : chi_heat;
  }
  else
    return 1.;
}

double ion_channel(double xe) {
  // return (1.-xe)/3.; // old approximation from Chen and Kamionkowski

  // coefficient as revised by Galli et al. 2013 (in fact it is a fit by Vivian Poulin of columns 1 and 4 in Table V of Galli et al. 2013)
  if (xe < 1.)
    return 0.369202*pow(1. - pow(xe, 0.463929), 1.70237);
  else
    return 0.;
}
