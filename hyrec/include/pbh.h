#ifndef PBH_H
#define PBH_H

/* CLASS common_pbh.h needed for BSPLINE definition */
#include "common_pbh.h"

/* Structure that is equivalent to the old pbh_external. Needed in HyRec separately
   as now pbh_exernal is external_info and contains stack-allocated HyRec arrays which
   are not needed internally in HyRec as part of a pure PBH structure! */
typedef struct PBH {
  /* structures containing b-spline information */
  BSPLINE * hion; /**< hydrogen ionisation 2d b-spline */
  BSPLINE * excite; /**< hydrogen excitation 2d b-spline */
  BSPLINE * heat; /**< plasma heating 2d b-spline */

  /* axes arrays and their sizes */
  double * masses; /**< pbh masses given in units of \f$ 10^{10} \f$ g in log10 space. must be monotonically decreasing */
  double * z_deps; /**< deposition redshifts for pbh energy injections */
  int masses_size; /**< size of masses array */
  int z_deps_size; /**< size of z_deps array */
} PBH;

#endif
