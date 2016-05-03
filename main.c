//
//  Population Synthesis of Planet Formation
//
//  main.c
//
//
//  Created by Takanori Sasaki
//  2016/05/02
//
//

#include "header.h"

int main(void)
{
  int i, j, k, n;
  double t, dt;


  // parameter
  int nd, nt, PS;
  double Time, FeH, alpha, qd, qg, C1, C2, M, L, p1, p2, kappa, k1, k2;

  parameter(&nd, &nt, &PS, &Time, &FeH, &alpha, &qd, &qg, &C1, &C2, &M, &L, &p1, &p2, &kappa, &k1, &k2);


  // variable of disk
  double r[nd], Sigd[nd], Sigg[nd], fd[nd], T[nd], eta[nd], fg, fg_0, tau_dep;

  // variable of planet
  int type[200], gene[200];
  double Mp[200], Mc[200], Mr[200], Mi[200], Mg[200], a[200], a_0[200], e[200], peri[200], apo[200];


  // calculations for each disk
  for (k=0; k<PS; k++) {

    initial(&n, &t, &dt, &tau_dep, &fg_0, &FeH, &nd, r, T, &L, eta, &fg, fd, Sigd, &qd, Sigg, &qg, a, &M, Mr, Mi, Mg, Mc, Mp, a_0, type, gene);

    //printf ("%d %f\n", k, fg_0);

    for (j=0; j<nt; j++) {

      core_acc(&n, type, r, a, Mp, &M, &nd, eta, fd, &fg, &qg, &qd, &dt, Mr, Mi, Mc, Mg, Sigd, &p1, &p2, &kappa);

      gas_acc(&n, type, &alpha, a, &L, &M, &k1, Mp, &k2, &fg, Mg, &dt, Mc);

      typeI_migration(&n, type, &qg, &C1, &fg, Mp, a, &M, &dt, a_0, Mr, gene, Mi, Mg, Mc, &nd, r, T, Sigg, &L);

      typeII_migration(&n, type, &fg, &C2, &alpha, Mp, a, &M, &dt);

      trap(&n, a, Mp, &M, type);

      next(&fg, &fg_0, &t, &tau_dep, &nd, Sigg, r, &qg, &dt, &Time, &nt, &j);

    }

    gimpact(e, a, Mp, peri, apo, &M, &t, &n, Mr, Mi);

    output(&n, a, Mp);

    //for (i=0; i<n; i++) printf("%f %f\n", a[i], Mp[i]/ME);

  }

  return 0;

}
