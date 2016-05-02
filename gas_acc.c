//
//  Population Synthesis of Planet Formation
//
//  gas_acc.c
//
//
//  Created by Takanori Sasaki
//  2016/05/02
//
//

#include "header.h"

void gas_acc(int *n, int type[], double *alpha, double a[], double *L, double *M, double *k1, double Mp[], double *k2, double *fg, double Mg[], double *dt, double Mc[])
{
  int i;
  double Mg_vis, Mg_th, tau_KH, dM_disk, dM, fgap;

  for (i=0; i<*n; i++) {

    // type1,2,3: embryo with gac accretion
    if (type[i] == 1 || type[i] == 2 || type[i] == 3) {

      // viscous condition
      Mg_vis = 30.0*(*alpha/1.0e-3)*pow(a[i], 0.5)*pow(*L/LS, 0.25)*ME;
      // thermal condition
      Mg_th = 0.95e3*pow(a[i], 0.75)*pow(*L/LS, 3.0/8.0)*pow(*M/MS, -0.5)*ME;


      // gas accretion
      tau_KH = pow(10.0, *k1)*pow(Mp[i]/ME, -(*k2));
      dM_disk = 3.0e-9*(*fg)*(*alpha/1.0e-3)*(*M);

      if (Mp[i]/tau_KH < dM_disk) dM = Mp[i]/tau_KH;
      else dM = dM_disk;

      if (Mp[i] < Mg_vis) fgap = 1.0;
      else if (Mp[i] > Mg_th) fgap = 0.0;
      else fgap = (log10(Mp[i])-log10(Mg_vis))/(log10(Mg_th)-log10(Mg_vis));

      Mg[i] += (*dt)*fgap*dM;
      Mp[i] = Mc[i] + Mg[i];

      // type2: initiation of gap formation; onset of Type II migration
      if (Mp[i] > Mg_vis) type[i] = 2;

      // type3: clear gap formation; termination of gas accretion
      if (Mp[i] > Mg_th) {
        Mp[i] = Mg_th;
        Mg[i] = Mp[i] - Mc[i];
        type[i] = 3;
      }

    }
  }
}
