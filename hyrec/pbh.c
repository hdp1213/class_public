/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         pbh.c: all functions related to PBH physics                                           */
/*                                                                                               */
/*         Units used: cgs + eV (all temperatures in eV)                                         */
/*                                                                                               */
/*         Version:     2017                                                                     */ 
/*                                                                                               */
/*         Revision history:                                                                     */
/*            - Written October 2017                                                             */
/*************************************************************************************************/ 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "include/pbh.h"

void allocate_and_read_pbh(PBH *pbh) {
  /* Axes */
  FILE *fA = fopen(PBH_AXES_FILE, "r");
  if (fA == NULL) {
    fprintf(stderr, "\033[1m\033[31m error\033[22;30m in allocate_and_read_pbh: could not open file ");
    fprintf(stderr, PBH_AXES_FILE);
    fprintf(stderr, "\n");
    exit(1);
  }

  ErrorMsg errmsg;

  // These routines allocate memory to the input arrays, so don't need to worry.
  class_read_1d_array(fA, &(pbh->z_deps), &(pbh->z_deps_size), errmsg);
  class_read_1d_array(fA, &(pbh->masses), &(pbh->masses_size), errmsg);

  fclose(fA);

  pbh->z_min = pbh->z_deps[pbh->z_deps_size-1];
  pbh->z_max = pbh->z_deps[0];

  /* Allocate spline memory */
  pbh->hion = (BSPLINE *) malloc(sizeof(BSPLINE));
  if (pbh->hion == NULL) {
    fprintf(stderr, "Could not allocate memory to pbh->hion\n");
    exit(1);
  }

  pbh->excite = (BSPLINE *) malloc(sizeof(BSPLINE));
  if (pbh->excite == NULL) {
    fprintf(stderr, "Could not allocate memory to pbh->excite\n");
    exit(1);
  }

  pbh->heat = (BSPLINE *) malloc(sizeof(BSPLINE));
  if (pbh->heat == NULL) {
    fprintf(stderr, "Could not allocate memory to pbh->heat\n");
    exit(1);
  }

  /* H ionisation */
  FILE *fhion = fopen(PBH_BSPLINE_HION_FILE, "r");
  if (fhion == NULL) {
    fprintf(stderr, "\033[1m\033[31m error\033[22;30m in allocate_and_read_pbh: could not open file ");
    fprintf(stderr, PBH_BSPLINE_HION_FILE);
    fprintf(stderr, "\n");
    exit(1);
  }

  class_read_bicubic_bspline(fhion, pbh->hion, errmsg);
  fclose(fhion);


  /* Ly-alpha excitation */
  FILE *fexcite = fopen(PBH_BSPLINE_EXCITE_FILE, "r");
  if (fexcite == NULL) {
    fprintf(stderr, "\033[1m\033[31m error\033[22;30m in allocate_and_read_pbh: could not open file ");
    fprintf(stderr, PBH_BSPLINE_EXCITE_FILE);
    fprintf(stderr, "\n");
    exit(1);
  }

  class_read_bicubic_bspline(fexcite, pbh->excite, errmsg);
  fclose(fexcite);


  /* Heating */
  FILE *fheat = fopen(PBH_BSPLINE_HEAT_FILE, "r");
  if (fheat == NULL) {
    fprintf(stderr, "\033[1m\033[31m error\033[22;30m in allocate_and_read_pbh: could not open file ");
    fprintf(stderr, PBH_BSPLINE_HEAT_FILE);
    fprintf(stderr, "\n");
    exit(1);
  }

  class_read_bicubic_bspline(fheat, pbh->heat, errmsg);
  fclose(fheat);
}

void free_pbh(PBH *pbh) {
  // First, free spline arrays
  free(pbh->hion->xknots);
  free(pbh->hion->yknots);
  free(pbh->hion->coeffs);

  free(pbh->excite->xknots);
  free(pbh->excite->yknots);
  free(pbh->excite->coeffs);

  free(pbh->heat->xknots);
  free(pbh->heat->yknots);
  free(pbh->heat->coeffs);

  // Then, the splines themselves
  free(pbh->hion);
  free(pbh->excite);
  free(pbh->heat);

  // Then the other arrays
  free(pbh->masses);
  free(pbh->z_deps);
}
