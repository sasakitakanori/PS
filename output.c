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

void output(int *n, double a[], double Mp[])
{
  FILE *fp;
  int i;

  if ((fp = fopen("a_Mp.dat", "a")) == NULL) {
    printf("cannot open file\n");
    exit(1);
  }
  for (i=0; i<*n; i++) fprintf(fp, "%f %f\n", a[i], Mp[i]/ME);
  fprintf(fp, "\n");
  fclose(fp);

}
