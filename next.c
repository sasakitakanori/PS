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

void next(double *fg, double *fg_0, double *t, double *tau_dep, int *nd, double Sigg[], double r[], double *qg, double *dt, double *Time, int *nt, int *j)
{
  int i;

  *fg = *fg_0*exp(-(*t)/(*tau_dep));
  for (i=0; i<*nd; i++) Sigg[i] = 7.5e1*(*fg)*pow(r[i]/10.0, -(*qg));

  *dt = pow(pow(*Time, 1.0/(double)(*nt)), (double)(*j)) - (*t);
  *t = pow(pow(*Time, 1.0/(double)(*nt)), (double)(*j));
}
