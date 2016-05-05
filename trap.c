//
//  Population Synthesis of Planet Formation
//
//  trap.c
//
//
//  Created by Takanori Sasaki
//  2016/05/03
//
//

#include "header.h"

void trap(int *n, double a[], double Mp[], double *M, int type[])
{
  int i, j, k;
  double a_m, rH, l1, l2, ll, ll2;

  for (i=0; i<(*n)-1; i++) {

    if (type[i] == 0 || type[i] == 1) {

      j = i+1;
      for (k=i+2; k<(*n); k++) {
        if (fabs(a[i] - a[k]) < fabs(a[i] - a[j])) j = k;
      }


      a_m = sqrt(a[i]*a[j]);
      rH = pow((Mp[i]+Mp[j])/(3.0*(*M)), 1.0/3.0)*a_m;


      if (fabs(a[i] - a[j]) < 5*rH) {

        if (type[i] == 4) type[j] = 4;
        else if (type[j] == 4) type[i] = 4;

        if (type[i] == 4 || type[j] == 4) {
          // da < 5rH -> da = 5rH
          if (a[i] < a[j]) a[j] = a[i] + 5*rH;
          else a[i] = a[j] + 5*rH;
        }

        else {
          // angular momentum before the trap
          l1 = Mp[i]*sqrt(a[i]*G*(*M));
          l2 = Mp[j]*sqrt(a[j]*G*(*M));
          ll = l1 + l2;

          // da < 5rH -> da = 5rH
          if (a[i] < a[j]) a[j] = a[i] + 5*rH;
          else a[i] = a[j] + 5*rH;

          // compensation for conservation of angular momentum
          do {
            a[i] -= 1.0e-3;
            a[j] -= 1.0e-3;

            l1 = Mp[i]*sqrt(a[i]*G*(*M));
            l2 = Mp[j]*sqrt(a[j]*G*(*M));
            ll2 = l1 + l2;
          } while (ll2 > ll);

        }
      }
    }
  }
}
