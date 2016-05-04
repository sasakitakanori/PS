//
//  Population Synthesis of Planet Formation
//
//  initial.c
//
//
//  Created by Takanori Sasaki
//  2016/05/02
//
//

#include "header.h"

void initial(int *n, double *t, double *dt, double *tau_dep, double *fg_0, double *FeH, int *nd, double r[], double T[], double *L, double eta[], double *fg, double fd[], double Sigd[], double *qd, double Sigg[], double *qg, double a[], double *M, double Mr[], double Mi[], double Mg[], double Mc[], double Mp[], double a_0[], int type[], int gene[], double *dM, double *alpha)
{
  int i;
  double logtau_dep, logf, f_disk, fd_0, Mc_iso, da, Sigg_v, Sigg_i, T_v, T_i, r_sv, r_si, r_snow, f_ice, h, hd;


  *n = 0;
  *t = 0.0;
  *dt = 1.0;


  // tau_dep: log-uniform distributions in the rage 10e6-10e7
  logtau_dep = 6.0 + (double)rand()/((double)RAND_MAX+1);
  *tau_dep = pow(10.0, logtau_dep);

  // f_disk: longnormal distributions centered on 1 with a dispersion of 1
  //         upper cutoff at 30 & lower cutoff at 0.03
  do {
    logf = sqrt(-2.0*log(((double)rand()+1.0)/((double)RAND_MAX+2.0)))*sin(2.0*M_PI*((double)rand()+1.0)/((double)RAND_MAX+2.0));
    f_disk = pow(10.0, logf);
  } while (f_disk > 30 || f_disk < 0.03);

  hd = pow(*M/MS, 2.0);

  *fg_0 = f_disk*hd;
  *fg = *fg_0;
  *dM = 3.0e-9*(*fg);

  r_sv = 1.2*pow(*M/MS, 1.0/3.0)*pow(*alpha/1.0e-3, -2.0/9.0)*pow(*dM/1.0e-8, 4.0/9.0);
  r_si = 0.75*pow(*L/LS, 2.0/3.0)*pow(*M/MS, -1.0/3.0);
  if (r_sv > r_si) r_snow = r_sv;
  else r_snow = r_si;

  // initial conditions for disk
  for (i=0; i<*nd; i++) {
    r[i] = 1.0e-2 + 1.0e-3*(double)i;

    Sigg_v = 2.1e3*pow(*M/MS, 1.0/5.0)*pow(*alpha/1.0e-3, -4.0/5.0)*pow(*dM/1.0e-8, 3.0/5.0)*pow(r[i], -3.0/5.0);
    Sigg_i = 2.7e3*pow(*L/LS, -2.0/7.0)*pow(*M/MS, 9.0/14.0)*pow(*alpha/1.0e-3, -1.0)*(*dM/1.0e-8)*pow(r[i], -15.0/14.0);
    if (Sigg_v < Sigg_i) Sigg[i] = Sigg_v;
    else Sigg[i] = Sigg_i;

    T_v = 2.0e2*pow(*M/MS, 3.0/10.0)*pow(*alpha/1.0e-3, -1.0/5.0)*pow(*dM/1.0e-8, 2.0/5.0)*pow(r[i], -9.0/10.0);
    T_i = 1.5e2*pow(*L/LS, 2.0/7.0)*pow(*M/MS, -1.0/7.0)*pow(r[i], -3.0/7.0);
    if (T_v > T_i) T[i] = T_v;
    else T[i] = T_i;

    if (T[i] > 1.7e2) eta[i] = 1.0;
    else {
      h = 5.0e-2*sqrt(T[i]/3.0e2)*pow(r[i], 1.5)*pow(*M/MS, -0.5);
      f_ice = 1.0 + 2.0*exp(-pow(r[i] - r_snow, 2.0)/(h*h));  // accumulation of icy grain at snow line
      eta[i] = 4.2*f_ice;
    }
  }

  fd_0 = pow(10.0, *FeH)*f_disk*hd;
  for (i=0; i<*nd; i++) {
    fd[i] = fd_0;
    Sigd[i] = 0.32*eta[i]*fd[i]*pow(r[i]/10.0, -(*qd));
  }


  // initial conditions for planet
  a[0] = 0.5;

  do {
    Mc_iso = 0.16*1.0*pow(fd_0, 1.5)*pow(a[*n], 0.75-1.5*(*qd-1.5))*pow(*M/MS, -0.5)*ME;
    da = 5.0*pow(2.0*Mc_iso/(3.0*(*M)), 1.0/3.0)*a[*n];
    a[(*n)+1] = a[*n] + 2.0*da;
    (*n)++;
  } while (a[(*n)-1] < r_snow);

  do {
    Mc_iso = 0.16*pow(4.2, 1.5)*pow(fd_0, 1.5)*pow(a[*n], 0.75-1.5*(*qd-1.5))*pow(*M/MS, -0.5)*ME;
    da = 5.0*pow(2.0*Mc_iso/(3.0*(*M)), 1.0/3.0)*a[*n];
    a[(*n)+1] = a[*n] + 2.0*da;
    (*n)++;
  } while (a[(*n)-1] < 10.0);

  for (i=0; i<*n; i++){
    Mr[i] = 1.0e20;  // rock mass
    Mi[i] = 0.0;     // ice mass
    Mg[i] = 0.0;     // gas mass
    Mc[i] = Mr[i] + Mi[i];  // core mass = rock mass + ice mass
    Mp[i] = Mc[i] + Mg[i];   // planet mass = core mass + gas mass
    a_0[i] = a[i];   // initial radius
    type[i] = 0;     // 0:embryo  1:onset of gas accretion 2:onset of Type II migration 3:termination of gas accretion 4:termination of accretion
    gene[i] = 1;     // generation
  }

}
