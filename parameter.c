//
//  Population Synthesis of Planet Formation
//
//  parameter.c
//
//
//  Created by Takanori Sasaki
//  2016/05/02
//
//

#include "header.h"

void parameter(int *nd, int *nt, int *PS, double *Time, double *FeH, double *alpha, double *qd, double *qg, double *C1, double *C2, double *M, double *L, double *p1, double *p2, double *kappa, double *k1, double *k2)
{

*nd = 30000;  //number of grid of disk; disk radius = nd*1.0e-3 AU
*nt = 10000;  //number of grid of time
*PS = 1000;  //number of calculations
*Time = 2.0e7;  //calculation time
*FeH = 0.1;  //metalicity (Fe/H)
*alpha = 1.0e-3;  //alpha viscosity
*qd = 1.5;  //power of dust density
*qg = 1.0;  //power of gas density
*C1 = 0.01;  //reduction factor of Type I migration
*C2 = 0.1;  //reduction factor of Type II migration
*M = MS;  //star's mass
*L = pow(*M/MS, 4.0)*LS;  //star's luminosity
*p1 = 0.25;  //power-law index for Mc_hydro
*p2 = 0.25;  //power-law index for Mc_hydro
*kappa = 1.0;  //absorption coefficient of dust
*k1 = 9.0;  //power-law index for tau_KH
*k2 = 3.0;  //power-law index for tau_KH

}
