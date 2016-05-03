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

void initial(int *n, double *t, double *dt, double *tau_dep, double *fg_0, double *FeH, int *nd, double r[], double T[], double *L, double eta[], double *fg, double fd[], double Sigd[], double *qd, double Sigg[], double *qg, double a[], double *M, double Mr[], double Mi[], double Mg[], double Mc[], double Mp[], double a_0[], int type[], int gene[])
{
  int i;
  double logtau_dep, logf, f_disk, fd_0, a_ice, Mc_iso, da;


  *n = 0;
  *t = 0.0;
  *dt = 1.0;


  // tau_dep: log-uniform distributions in the rage 10e6-10e7
  logtau_dep = 6.0 + (double)rand()/((double)RAND_MAX+1);
  *tau_dep = pow(10.0, logtau_dep);

  // f_disk: longnormal distributions centered on 1 with a dispersion of 1
  //         upper cutoff at 30 & lower cutoff at 0.01
  do {
    logf = sqrt(-2.0*log(((double)rand()+1.0)/((double)RAND_MAX+2.0)))*sin(2.0*M_PI*((double)rand()+1.0)/((double)RAND_MAX+2.0));
        f_disk = pow(10.0, logf);
  } while (f_disk > 30 || f_disk < 0.01);

  *fg_0 = f_disk;
  fd_0 = pow(10.0, *FeH)*f_disk;


  // initial conditions for disk
  for (i=0; i<*nd; i++) {
    r[i] = 1.0e-2 + 1.0e-3*(double)i;
    T[i] = 2.8e2*pow(r[i], -0.5)*pow(*L/LS, 0.25);
  }

  a_ice = 2.7*pow(*L/LS, 0.5);
  for (i=0; i<*nd; i++) {
    if (r[i] < a_ice) eta[i] = 1.0;
    else eta[i] = 4.2;
  }

  *fg = *fg_0;
  for (i=0; i<*nd; i++) {
    fd[i] = fd_0;
    Sigd[i] = 0.32*eta[i]*fd[i]*pow(r[i]/10.0, -(*qd));
    Sigg[i] = 7.5e1*(*fg_0)*pow(r[i]/10.0, -(*qg));
  }


  // initial conditions for planet
  a[0] = 0.5;

  do {
    Mc_iso = 0.16*1.0*pow(fd_0, 1.5)*pow(a[*n], 0.75-1.5*(*qd-1.5))*pow(*M/MS, -0.5)*ME;
    da = 5.0*pow(2.0*Mc_iso/(3.0*(*M)), 1.0/3.0)*a[*n];
    a[(*n)+1] = a[*n] + 2.0*da;
    (*n)++;
  } while (a[(*n)-1] < a_ice);

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
        type[i] = 0;       // 0:embryo  1:
        gene[i] = 1;     // generation
    }

}
