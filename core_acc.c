//
//  Population Synthesis of Planet Formation
//
//  core_acc.c
//
//
//  Created by Takanori Sasaki
//  2016/05/02
//
//

#include "header.h"

void core_acc(int *n, int type[], double r[], double a[], double Mp[], double *M, int *nd, double eta[], double fd[], double *fg, double *qg, double *qd, double *dt, double Mr[], double Mi[], double Mc[], double Mg[], double Sigd[], double *p1, double *p2, double *kappa)
{
  int i, j, k;
  int in, out;
  double dMc, da, tau_acc, Mc_hydro;

  for (i=0; i<*n; i++) {

    // type4: no more planetesimals and gas around the embryo
    if (type[i] != 4) {


      // set the embryo's position and accretion regeion
      j = 0;
      do {
        j++;
      } while (r[j] < a[i]);
      k = j;  // k: embryo's position

      in = 0;
      out = 0;
      da = 5.0*pow(2.0*Mp[i]/(3.0*(*M)), 1.0/3.0)*a[i];

      j = 0;
      do {
        in++;
        j++;
      } while((r[k-j] > a[i] - da) && (k-j > 0));

      j = 0;
      do {
        out++;
        j++;
      } while((r[k+j] < a[i] + da) && (k+j < *nd));


      // core accretion
      tau_acc = 2.2e5*pow(eta[k], -1.0)*pow(fd[k], -1.0)*pow(*fg, -0.4)*pow(10.0, 0.4*(1.5-(*qg)))*pow(10.0, 1.5-(*qd))*pow(a[i], 2.7+(*qd-1.5)+0.4*(*qg-1.5))*pow(Mp[i]/ME, 1.0/3.0)*pow(*M/MS, -1.0/6.0);

      dMc = *dt*Mp[i]/tau_acc;
      Mr[i] += dMc/eta[k];
      if (eta[k] > 1.0) Mi[i] += dMc*(eta[k]-1.0)/eta[k];
      Mc[i] = Mr[i] + Mi[i];
      Mp[i] = Mc[i] + Mg[i];


      // reducing the disk mass as the core accretion
      for (j=k-in; j<k+out; j++) {
        Sigd[j] = (Sigd[j]*M_PI*AU*AU*(r[j+1]*r[j+1]-r[j]*r[j]) - dMc/(in+out+1))/(M_PI*AU*AU*(r[j+1]*r[j+1]-r[j]*r[j]));
        if (Sigd[j] < 0.0) {
          Mc[i] += Sigd[j]*M_PI*AU*AU*(r[j+1]*r[j+1]-r[j]*r[j]);
          Sigd[j] = 0.0;
        }
        fd[j] = Sigd[j]/(0.32*eta[j]*pow(r[j]/10.0, -(*qd)));
      }


      // type1: initiation of gas accretion
      if (type[i] == 0) {
        Mc_hydro = 10.0*pow((dMc/(*dt))/(1.0e-6*ME), *p1)*pow(*kappa/1.0, *p2)*ME;
        if (Mp[i] > Mc_hydro) type[i] = 1;
      }

    }
  }
}
