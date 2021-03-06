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

void typeI_migration(int *n, int type[], double *qg, double *C1, double *fg, double Mp[], double a[], double *M, double *dt, double a_0[], double Mr[], int gene[], double Mi[], double Mg[], double Mc[], int *nd, double r[], double T[], double Sigg[], double *L)
{
  int i, j;
  double tau_mig;
  double dr, p, q, Omega, Cs;

  for (i=1; i<*n; i++) {

    if (type[i] == 0 || type[i] == 1) {

      // Type I migratiion rate by Tanaka et al. 2002
      for (j=0; j<*nd-1; j++) {
        if (a[i] > r[j] && a[i] < r[j+1]) {
          p = (log(Sigg[j+1])-log(Sigg[j]))/(log(r[j+1])-log(r[j]));
          q = (log(T[j+1])-log(T[j]))/(log(r[j+1])-log(r[j]));
          Omega = pow(G*(*M)/pow(a[i]*AU, 3.0), 0.5);
          Cs = 1.0e5*pow(T[j]/3.0e2, 0.5)/AU;
          dr = *C1*1.08*(p + 0.80*q - 2.52)*(Mp[i]/(*M))*(Sigg[j]*AU*AU*a[i]*a[i]/(*M))*pow(a[i]*Omega/Cs, 2.0)*a[i]*Omega;
        }
      }
      a[i] += (*dt*YEAR)*dr;


      // next generation of the embryo
      if (a[i] < 0.5*a_0[i] && gene[i] < 5) {

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
      if (a[i] < 0.01) {
        a[i] = 0.0;//0.04;
        type[i] = 4;
      }

    }
  }
}
