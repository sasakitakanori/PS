//
//  Population Synthesis of Planet Formation
//
//  typeII_migration.c
//
//
//  Created by Takanori Sasaki
//  2016/05/02
//
//

#include "header.h"

void typeII_migration(int *n, int type[], double *fg, double *C2, double *alpha, double Mp[], double a[], double *M, double *dt)
{
  int i;
  double tau_mig2p, tau_mig2d, tau_mig;

  for (i=0; i<*n; i++) {

    if (type[i] == 2 || type[i] == 3) {

      tau_mig2p = 5.0e5*pow(*fg, -1.0)*pow(*C2*(*alpha)/1.0e-4, -1.0)*(Mp[i]/MJ)*pow(a[i], 0.5)*pow(*M/MS, -0.5);
      tau_mig2d = 0.7e5*pow((*alpha)/1.0e-3, -1.0)*a[i]*pow(*M/MS, -0.5);
      if (tau_mig2d > tau_mig2p) tau_mig = tau_mig2d;
      else tau_mig = tau_mig2p;
      a[i] -= *dt*a[i]/tau_mig;

      // type4: no more planetesimals and gas around the embryo
      if (a[i] < 0.01) {
        a[i] = 0.0;//0.04;
        type[i] = 4;
      }

    }
  }
}
