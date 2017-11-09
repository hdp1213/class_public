/** @file pbh.h Generic parameters and functions used in PBH parts of the code. */

#ifndef __PBH__
#define __PBH__

#include "common.h"

/**
 * structure containing knot points and coefficients for a 2d b-spline
 */

struct bspline_2d {

  int nxknots; /**< number of knots in the x direction */
  double * xknots; /**< x-coordinate values of these knots */

  int nyknots; /**< number of knots in the y direction */
  double * yknots; /**< y-coordinate values of these knots */

  double * coeffs; /**< coefficients describing spline **/

  int degree; /**< degree of b-spline. 1 for bilinear splines, 3 for bicubic splines */

};

typedef struct bspline_2d BSPLINE;

/**
 * Structure containing externally-read b-splines and axes needed for PBH energy injection computation
 *
 * Should only be used by an external CLASS engine
 */

struct pbh_external {

  /** @name - parameters to initialise PBH externally */

  //@{

  /* structures containing b-spline information */
  struct bspline_2d * hion; /**< hydrogen ionisation 2d b-spline */
  struct bspline_2d * excite; /**< hydrogen excitation 2d b-spline */
  struct bspline_2d * heat; /**< plasma heating 2d b-spline */

  /* axes arrays and their sizes */
  double * masses; /**< pbh masses given in units of \f$ 10^{10} \f$ g in log10 space. must be monotonically decreasing */
  double * z_deps; /**< deposition redshifts for pbh energy injections */
  int masses_size; /**< size of masses array */
  int z_deps_size; /**< size of z_deps array */

  //@}
};

typedef struct pbh_external PBH;

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int class_read_2d_array(
                          FILE * fp,
                          double * array,
                          int x_size,
                          int y_size,
                          ErrorMsg error_message
                          );

  int class_read_1d_array(
                          FILE * fp,
                          double ** array,
                          int * array_size,
                          ErrorMsg error_message
                          );

  int class_read_bicubic_bspline(
                                 FILE * fp,
                                 struct bspline_2d * pbsp,
                                 ErrorMsg error_message
                                 );

#ifdef __cplusplus
}
#endif

#endif
