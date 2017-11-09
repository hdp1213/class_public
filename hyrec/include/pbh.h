#ifndef PBH_H
#define PBH_H
/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         pbh.h: functions for pbh stuff                                                        */
/*                                                                                               */
/*************************************************************************************************/

#include "hyrectools.h"
#include "common_pbh.h"

#define PBH_AXES_FILE             "hyrec/data/pbh/axes.dat"
#define PBH_BSPLINE_HION_FILE     "hyrec/data/pbh/bspline_hion.dat"
#define PBH_BSPLINE_EXCITE_FILE   "hyrec/data/pbh/bspline_excite.dat"
#define PBH_BSPLINE_HEAT_FILE     "hyrec/data/pbh/bspline_heat.dat"

typedef struct {
  /* structures containing b-spline information */
  BSPLINE *hion; /**< hydrogen ionisation 2d b-spline */
  BSPLINE *excite; /**< hydrogen excitation 2d b-spline */
  BSPLINE *heat; /**< plasma heating 2d b-spline */

  /* axes arrays and their sizes */
  double *masses; /**< pbh masses given in units of \f$ 10^{10} \f$ g in log10 space. must be monotonically decreasing */
  double *z_deps; /**< deposition redshifts for pbh energy injections */
  int masses_size; /**< size of masses array */
  int z_deps_size; /**< size of z_deps array */

  /* Useful to know in a hurry */
  double z_min;
  double z_max;

} PBH;

void allocate_and_read_pbh(PBH *pbh);
void free_pbh(PBH *pbh);

#endif
