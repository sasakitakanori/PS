//
//  Population Synthesis of Planet Formation
//
//  typeI_migration.c
//
//
//  Created by Takanori Sasaki
//  2016/05/02
//
//

#include "header.h"

void typeI_migration(int *n, int type[], double *qg, double *C1, double *fg, double Mp[], double a[], double *M, double *dt, double a_0[], double Mr[], int gene[], double Mi[], double Mg[], double Mc[])
{
  int i;
  double tau_mig;

  for (i=0; i<*n; i++) {

    if (type[i] == 0 || type[i] == 1) {

      tau_mig = 5.0e4*pow(10.0, 1.5-(*qg))*pow(*C1*(*fg), -1.0)*pow(Mp[i]/ME, -1.0)*pow(a[i], *qg)*pow(*M/MS, 1.5);
      a[i] -= *dt*a[i]/tau_mig;


      // next generation of the embryo
      if (a[i] < 0.5*a_0[i]) {

        a[*n] = a_0[i];
        a_0[*n] = a_0[i];
        Mr[*n] = 1.0e20*pow(0.01, (double)gene[i]-1.0);
        Mi[*n] = 0.0;
        Mg[*n] = 0.0;
        Mc[*n] = Mr[*n] + Mi[*n];
        Mp[*n] = Mc[*n] + Mg[*n];
        type[*n] = 0;
        gene[*n] = gene[i] + 1;
        (*n)++;

        a_0[i] = 0.0;
      }

      // type4: no more planetesimals and gas around the embryo
      if (a[i] < 0.04) {
        a[i] = 0.04;
        type[i] = 4;
      }

    }
  }
}
