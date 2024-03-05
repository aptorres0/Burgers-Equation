/**
 * burger.c
 *
 * Author: Alexander Paul Torres
 * Date: 30 JAN 23
 * Class: PHYS 7411 - Computational Physics I
 *
 * Description:
 * Computes solutions to the 1D First-order Non-linear Convection - Burger's
 * Equation. Solves the equation on a mesh with a spatial component on the
 * interval 0 <= x <= x_max and a time interval from 0 to t_max. Time step size
 * determined by the Courant-Friedrichs-Lewy (CFL) condition for the convergence
 * of numerical schemes for hyperbolic PDEs. Implementation here imposes the
 * boundary condition that the wave height at x = 0 and x = x_max is equal to 0.
 *
 * Two numerical schemes are implemented:
 *  - Leapfrog/FTCS Method: Finite-difference approximation that uses a forward
 *                          in time central in space discretization of he mesh.
 *                          Accurate to first order in time and second order in
 *                          space
 *  - Lax-Wendroff Method: Finite-difference method which is accurate to second
 *                         order in both space and time.
 *
 * Output: Both methods output txt files for plotting the wave height
 *
 */
#include "burger.h"

#include <math.h>    // sin, M_PI
#include <stdio.h>   // FILE
#include <stdlib.h>  // exit
#include <string.h>  // memcpy

void LeapFrogMethod(double u0[], int N, double x_max, double t_max, double beta,
                    char *fname) {
  // Discretize the mesh
  double dx = x_max / (N - 1.0);    // space step size
  double epsilon = 1.;              // wave speed [dimensionless]
  double dt = beta * dx / epsilon;  // time step size
  int M = (int)(t_max / dt) + 1;    // number of time points

  // Create output file to store results and write out initial values
  FILE *fptr = fopen(fname, "w");
  if (fptr == NULL) {
    printf("Error creating file %s!\n", fname);
    exit(1);
  }
  for (int i = 0; i < N; i++) {  // write out initial condition
    fprintf(fptr, "%f ", u0[i]);
  }
  fprintf(fptr, "\n");

  // Create array to hold values at j + 1
  double u[N];

  // Perform leapfrog method
  for (int j = 1; j < M; j++) {  // Loop through time
    // Apply boundary conditions
    u[0] = 0;
    u[N - 1] = 0;
    fprintf(fptr, "%f ", u[0]);        // Write out boundary condition at x = 0
    for (int i = 1; i < N - 1; i++) {  // Loop through space (inner points only)
      u[i] =
          u0[i] - 0.25 * beta * (u0[i + 1] * u0[i + 1] - u0[i - 1] * u0[i - 1]);
      fprintf(fptr, "%f ", u[i]);
    }
    fprintf(fptr, "%f\n",
            u[N - 1]);         // Write out boundary condition at x = xmax
    memcpy(u0, u, sizeof(u));  // Let u0 = u to get next time step
  }
  fclose(fptr);
}

void LaxWendroffMethod(double u0[], int N, double x_max, double t_max,
                       double beta, char *fname) {
  // Discretize the mesh
  double dx = x_max / (N - 1.0);    // space step size
  double epsilon = 1.;              // wave speed [dimensionless]
  double dt = beta * dx / epsilon;  // time step size
  int M = (int)(t_max / dt) + 1;    // number of time points

  // Create output file to store results and write out initial values
  FILE *fptr = fopen(fname, "w");
  if (fptr == NULL) {
    printf("Error creating file %s!\n", fname);
    exit(1);
  }
  for (int i = 0; i < N; i++) {  // write out initial condition
    fprintf(fptr, "%f ", u0[i]);
  }
  fprintf(fptr, "\n");

  // Create array to hold values at j + 1
  double u[N];

  // Perform leapfrog method
  for (int j = 1; j < M; j++) {  // Loop through time
    // Apply boundary conditions
    u[0] = 0;
    u[N - 1] = 0;
    fprintf(fptr, "%f ", u[0]);        // Write out boundary condition at x = 0
    for (int i = 1; i < N - 1; i++) {  // Loop through space (inner points only)
      u[i] =
          u0[i] -
          0.25 * beta * (u0[i + 1] * u0[i + 1] - u0[i - 1] * u0[i - 1]) +
          0.125 * beta * beta *
              ((u0[i + 1] + u0[i]) * (u0[i + 1] * u0[i + 1] - u0[i] * u0[i]) -
               (u0[i] + u0[i - 1]) * (u0[i] * u0[i] - u0[i - 1] * u0[i - 1]));
      fprintf(fptr, "%f ", u[i]);
    }
    fprintf(fptr, "%f\n",
            u[N - 1]);         // Write out boundary condition at x = xmax
    memcpy(u0, u, sizeof(u));  // Let u0 = u to get next time step
  }
  fclose(fptr);
}
