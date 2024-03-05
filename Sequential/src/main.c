/**
 * main.c
 *
 * Author: Alexander Paul Torres
 * Date: 30 JAN 23
 * Class: PHYS 7411 - Computational Physics I
 *
 * Description:
 * File to show implementation of burger.c functions
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "burger.h"

int main(int argc, char* argv[]) {
  if (argc < 2) {
    printf("Usage: %s beta tmax\nWhere beta and tmax are numbers", argv[0]);
    return 1;
  }

  // Define mesh and input parameters
  double beta = atof(argv[1]);  // $\beta$ is the CFL
  double t_max = atof(argv[2]);
  int N = 1001;
  double x_max = 2.;
  double dx = x_max / (N - 1.);
  double epsilon = 1.0;
  double dt = beta * dx / epsilon;
  int M = (int)(t_max / dt) + 1;

  // Initial arrays
  double u0_lw[N];  // for Lax-Wendroff method
  double u0_lf[N];  // for Leapfrog method
  u0_lf[0] = 0.;
  u0_lf[N - 1] = 0.;
  /*for (int i = 1; i < N - 1; i++) {
    u0_lf[i] = 3.0 * sin(M_PI * i * dx);
  }*/
  for (int i = 1; i < N - 1; i++) {
    u0_lf[i] = 4.0 * exp(-300.*pow((i*dx-0.12),2.0));
  }
  memcpy(u0_lw, u0_lf, sizeof(u0_lf));

  // Define output files
  char* fName_lf = "leapfrog-c-output.txt";
  char* fName_lw = "laxwendroff-c-output.txt";

  // Perform numerical methods
  LeapFrogMethod(u0_lf, N, x_max, t_max, beta, fName_lf);
  LaxWendroffMethod(u0_lw, N, x_max, t_max, beta, fName_lw);

  printf(
      "N = %d\tx_max = %f\tdx = %f\nM = %d\tt_max = %f\tdt = %f\nbeta = "
      "%f\nepsilon = "
      "%f\nOutput files:\n%s\n%s\n",
      N, x_max, dx, M, t_max, dt, beta, epsilon, fName_lw, fName_lf);

  return 0;
}