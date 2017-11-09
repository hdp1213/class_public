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
