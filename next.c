//
//  Population Synthesis of Planet Formation
//
//  next.c
//
//
//  Created by Takanori Sasaki
//  2016/05/02
//
//

#include "header.h"

void next(double *fg, double *fg_0, double *t, double *tau_dep, int *nd, double Sigg[], double r[], double *dt, double *Time, int *nt, int *j, double *dM, double *M, double *alpha, double *L, double T[], double eta[])
{
  int i;
  double Sigg_v, Sigg_i, T_v, T_i;

  *fg = *fg_0*exp(-(*t)/(*tau_dep));
  *dM = 3.0e-9*(*fg);

  for (i=0; i<*nd; i++) {
    Sigg_v = 2.1e3*pow(*M/MS, 1.0/5.0)*pow(*alpha/1.0e-3, -4.0/5.0)*pow(*dM/1.0e-8, 3.0/5.0)*pow(r[i], -3.0/5.0);
    Sigg_i = 2.7e3*pow(*L/LS, -2.0/7.0)*pow(*M/MS, 9.0/14.0)*pow(*alpha/1.0e-3, -1.0)*(*dM/1.0e-8)*pow(r[i], -15.0/14.0);
    if (Sigg_v < Sigg_i) Sigg[i] = Sigg_v;
    else Sigg[i] = Sigg_i;

    T_v = 2.0e2*pow(*M/MS, 3.0/10.0)*pow(*alpha/1.0e-3, -1.0/5.0)*pow(*dM/1.0e-8, 2.0/5.0)*pow(r[i], -9.0/10.0);
    T_i = 1.5e2*pow(*L/LS, 2.0/7.0)*pow(*M/MS, -1.0/7.0)*pow(r[i], -3.0/7.0);
    if (T_v > T_i) T[i] = T_v;
    else T[i] = T_i;

    if (T[i] > 1.7e2) eta[i] = 1.0;
    else eta[i] = 4.2;
  }

  *dt = pow(pow(*Time, 1.0/(double)(*nt)), (double)(*j)) - (*t);
  *t = pow(pow(*Time, 1.0/(double)(*nt)), (double)(*j));
}
