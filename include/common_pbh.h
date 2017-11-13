/** @file pbh.h Generic parameters and functions used in PBH parts of the code. */

#ifndef __PBH__
#define __PBH__

/* Needed for ErrorMsg definition */
#include "common.h"

/* Include macros for array sizes needed for external HyRec functionality */
#ifdef HYREC
#include "io_params.h"
#endif

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

/* Maintain cross-compatibility with HyRec */
typedef struct bspline_2d BSPLINE;

/**
 * Structure containing externally-read b-splines and axes needed for PBH energy injection computation
 *
 * Should only be used by an external CLASS engine
 */

struct external_info {

  /** @name - parameters to initialise PBH and HyRec externally */

  //@{

  /**** PBH physics ****/

  /* structures containing b-spline information */
  struct bspline_2d * hion; /**< hydrogen ionisation 2d b-spline */
  struct bspline_2d * excite; /**< hydrogen excitation 2d b-spline */
  struct bspline_2d * heat; /**< plasma heating 2d b-spline */

  /* axes arrays and their sizes */
  double * masses; /**< pbh masses given in units of \f$ 10^{10} \f$ g in log10 space. must be monotonically decreasing */
  double * z_deps; /**< deposition redshifts for pbh energy injections */
  int masses_size; /**< size of masses array */
  int z_deps_size; /**< size of z_deps array */

  /**** HyRec arrays ****/

#ifdef HYREC
  /* Tables of effective rates */
  double **logAlpha_tab[2];
  double logR2p2s_tab[NTR];

  /* Tables of 2-photon rates */
  double Eb_tab[NVIRT];       /* Energies of the virtual levels in eV */
  double A1s_tab[NVIRT];      /* 3*A2p1s*phi(E)*DE */ 
  double A2s_tab[NVIRT];      /* dLambda_2s/dE * DeltaE if E < Elya dK2s/dE * Delta E if E > Elya */
  double A3s3d_tab[NVIRT];    /* (dLambda_3s/dE + 5*dLambda_3d/dE) * Delta E for E < ELyb, Raman scattering rate for E > ELyb */
  double A4s4d_tab[NVIRT];    /* (dLambda_4s/dE + 5*dLambda_4d/dE) * Delta E */
#endif

  //@}
};

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
