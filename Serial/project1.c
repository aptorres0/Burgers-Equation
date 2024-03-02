/**
 * project1.c
 *
 * Author: Alexander Paul Torres
 * Date: 30 JAN 23
 * Class: PHYS 7411 - Computational Physics I
 *
 * Description:
 * Computes the inviscid Burger's equation of advection for an initial
 * sine wave. Applies both the Leapfrog method and Lax-Wendroff method
 * time-stepping schemes.
 *
 * Output: Two txt files for plotting
 *
 * NOTE: to compile
 *  gcc -DGlobal_Vars -std=gnu11 -g -O3 -o project1-c project1.c
 */
#include <math.h>    // sin, M_PI
#include <stdio.h>   // FILE
#include <stdlib.h>  // exit
#include <string.h>  // memcpy

#ifdef Global_Vars
double tmax = 0.3;
double beta = 1.2;  // CFL Number [Dimensionless]
#endif

void LeapFrogMethod();
void LaxWendroffMethod();

int main() {
  LeapFrogMethod();
  LaxWendroffMethod();

  return 0;
}

/**
 * Function: LeapFrogMethod
 * ------------------------
 * computes the Burger equation using a forward difference time derivative
 * and a central difference spatial derivative scheme. Accurate to second
 * order in space and first order in time.
 *
 * returns: writes the output to a file named c-output-Leapfrog.txt
 *
 */
void LeapFrogMethod() {
  // Create the mesh
#ifndef Global_Vars
  double tmax = 0.3;  // Upper boundary for time [unitless]
  double beta = 0.1;  // CFL Number [Dimensionless]; starting with 0.1
#endif
  float epsilon = 1.0;           // Normalized velocity [Length / Time]
  double xmax = 2.0;             // [dimensionless] spatial upper boundary
  int N = 101;                   // Number of grid points in space
  double dx = xmax / (N - 1.0);  // spatial step size
  /** Test
  double dt = beta * dx / epsilon;  // time step size
  int M = (int)(tmax / dt);         // number of grid points in time
  */
  int M = 150;
  double dt = tmax / (M - 1.0);
  // printf("%d\n", M);

  // Create output file to store results
  FILE *fptr = fopen("c-output-Leapfrog.txt", "w");
  if (fptr == NULL) {
    printf("Error creating file!");
    exit(1);
  }

  /**
   * Create arrays to hold values at times j (u0) and j + 1 (u),
   *  initalize u0[i] = 3 sin (pi*x) and write to output file
   *
   */
  double u[N];
  double u0[N];
  u0[0] = 0;
  u0[N - 1] = 0;
  fprintf(fptr, "%f ", u0[0]);
  for (int i = 1; i < N - 1; i++) {
    u0[i] = 3.0 * sin(M_PI * i * dx);
    fprintf(fptr, "%f ", u0[i]);
  }
  fprintf(fptr, "%f\n", u0[N - 1]);

  // Perform leapfrog method
  for (int j = 0; j < M; j++) {  // Loop through time
    // Apply boundary conditions
    u[0] = 0;
    u[N - 1] = 0;
    fprintf(fptr, "%f ", u[0]);        // Write out x = 0 B.C.
    for (int i = 1; i < N - 1; i++) {  // Loop through space (inner points only)
      u[i] =
          u0[i] - 0.25 * beta * (u0[i + 1] * u0[i + 1] - u0[i - 1] * u0[i - 1]);
      fprintf(fptr, "%f ", u[i]);
    }
    fprintf(fptr, "%f\n", u[N - 1]);  // Write out u[x = xmax] B.C.
    memcpy(u0, u, sizeof(u));         // Let u0 = u to get next time step
  }
  fclose(fptr);
}

/**
 * Function: LaxWendroffMethod
 * ---------------------------
 * computes the Burger's equation using a scheme accurate to second order in
 * both space and time.
 *
 * returns: writes the output to a file named c-output-LaxWendroff.txt
 *
 */
void LaxWendroffMethod() {
  /**
   * Create the mesh
   *    Number of spatial grid points: N = 101
   *    Spatial boundaries: 0 <= x <= 2
   *    Spatial dimension step size: dx = 2/(N-1) [Dimensionless]
   *    Temporal boundaries: 0 <= t <= 0.15
   *    CFL Number: beta = dt/dx [Dimensionless]
   *    Temporal step size: dt = beta*dx
   *    Number of temporal grid points: M = 0.15/dt
   *
   */
#ifndef Global_Vars
  double tmax = 0.3;
  double beta = 0.1;  // CFL Number [Dimensionless]; starting with 0.1
#endif
  float epsilon = 1.0;  // Normalized velocity [Length / Time]
  double xmax = 2.0;
  int N = 101;
  double dx = xmax / (N - 1.0);
  /* Test
  double dt = beta * dx / epsilon;
  int M = (int)(tmax / dt);
  */
  int M = 150;
  double dt = tmax / (M - 1.0);

  // Create output file to write results
  FILE *fptr = fopen("c-output-LaxWendroff.txt", "w");
  if (fptr == NULL) {
    printf("Error creating file!");
    exit(1);
  }

  /**
   * Create arrays to hold values at times j (u0) and j + 1 (u),
   *  initalize u0[i] = 3 sin (pi*x) and write to output file
   *
   */
  double u[N];   // Holds point j+1
  double u0[N];  // Holds point j
  u0[0] = 0;
  u0[N - 1] = 0;
  fprintf(fptr, "%f ", u0[0]);
  for (int i = 1; i < N - 1; i++) {
    u0[i] = 3.0 * sin(M_PI * i * dx);
    fprintf(fptr, "%f ", u0[i]);
  }
  fprintf(fptr, "%f\n", u0[N - 1]);

  /**
   * Perform Lax-Wendroff method
   *
   */
  for (int j = 0; j < M; j++) {  // Loop through time
    // Apply boundary conditions
    u[0] = 0;
    u[N - 1] = 0;
    fprintf(fptr, "%f ", u[0]);        // Write out x = 0 BC
    for (int i = 1; i < N - 1; i++) {  // Loop through space (inner points only)
      u[i] =
          u0[i] -
          0.25 * beta * (u0[i + 1] * u0[i + 1] - u0[i - 1] * u0[i - 1]) +
          0.125 * beta * beta *
              ((u0[i + 1] + u0[i]) * (u0[i + 1] * u0[i + 1] - u0[i] * u0[i]) -
               (u0[i] + u0[i - 1]) * (u0[i] * u0[i] - u0[i - 1] * u0[i - 1]));
      fprintf(fptr, "%f ", u[i]);
    }
    fprintf(fptr, "%f\n", u[N - 1]);  // Write out x = xmax BC
    // assign u0 = u to get new values at next time step
    memcpy(u0, u, sizeof(u));
  }
  fclose(fptr);
}

double LFtest(double u0[], double beta, int index) {
  return u0[index] -
         0.25 * beta *
             (u0[index + 1] * u0[index + 1] - u0[index - 1] * u0[index - 1]);
}

double LWtest(double u0[], double beta, int index) {
  return LFtest(u0, beta, index) +
         0.125 * beta * beta *
             ((u0[index + 1] + u0[index]) *
                  (u0[index + 1] * u0[index + 1] - u0[index] * u0[index]) -
              (u0[index] + u0[index - 1]) *
                  (u0[index] * u0[index] - u0[index - 1] * u0[index - 1]));
}