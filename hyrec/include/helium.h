#ifndef HELIUM_H
#define HELIUM_H
/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         helium.h: all functions related to helium recombination (and Saha equilibria)         */
/*                                                                                               */
/*         Version: 2015  (first released November 2010)                                         */
/*************************************************************************************************/ 

/* CLASS common.h for all good CLASS stuff */
#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif

double rec_xesaha_HeII_III(double nH0, double Tr0, double fHe, double z, double *xHeIII, double fsR, double meR);
double rec_saha_xHeII(double nH0, double Tr0, double fHe, double z, double fsR, double meR);
double rec_saha_xH1(double xHeII, double nH0, double T0, double z, double fsR, double meR);
int rec_helium_dxHeIIdlna(double xH1s, double xHeII, double nH0, double Tr0, double fHe,
                          double H, double z, double fsR, double meR, double *dxHeIIdlna,
                          ErrorMsg error_message);

#ifdef __cplusplus
}
#endif

#endif
