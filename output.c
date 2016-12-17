//
//  Population Synthesis of Planet Formation
//
//  output.c
//
//
//  Created by Takanori Sasaki
//  2016/05/03
//
//

#include "header.h"

void output(int *n, double a[], double Mp[], double Mr[], double Mi[], double Mg[])
{
  FILE *fp;
  int i, j;

  if ((fp = fopen("a_Mp_default_05_2.dat", "a")) == NULL) {
    printf("cannot open file\n");
    exit(1);
  }
  for (i=0; i<*n; i++) {
    j = 1;  // rocky planet
    if (Mi[i] > Mr[i]) j = 2;  // ice giant planet
    if (Mg[i] > Mi[i] && Mg[i] > Mr[i]) j = 3;  // gas giant planet
    fprintf(fp, "%f %f %f %f %f %d\n", a[i], Mp[i]/ME, Mr[i]/ME, Mi[i]/ME, Mg[i]/ME, j);
  }
  fprintf(fp, "\n");
  fclose(fp);

}
