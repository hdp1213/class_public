#ifndef HYREC_H
#define HYREC_H

#define HYREC_VERSION "2017"

#include "history.h"
/* CLASS common.h for all good CLASS stuff */
#include "common.h"


#ifdef __cplusplus
extern "C" {
#endif

int hyrec_allocate(HYREC_DATA *data, double zmax, double zmin, short read_atomic_files, ErrorMsg error_message);

void hyrec_free(HYREC_DATA *data, short read_atomic_files);

int hyrec_compute(HYREC_DATA *data, int model,
                   double h, double T0, double Omega_b, double Omega_m, double Omega_k, double YHe, double Nnueff,
                   double alphaR, double meR, double pann, double pann_halo, double ann_z, double ann_zmax,
                   double ann_zmin, double ann_var, double ann_z_halo, double Mpbh, double fpbh,
                   ErrorMsg error_message);

double hyrec_xe(double z, HYREC_DATA *rec_data);

double hyrec_Tm(double z, HYREC_DATA *rec_data);

double hyrec_dTmdlna(double z, HYREC_DATA *rec_data);

#ifdef __cplusplus
}
#endif

#endif
