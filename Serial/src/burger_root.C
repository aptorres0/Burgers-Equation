#include <math.h>    // sin, M_PI
#include <stdio.h>   // FILE
#include <stdlib.h>  // exit
#include <string.h>  // memcpy
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TFile.h"

void burger_root() {
  double beta = 0.1;  // $\beta$ is the CFL
  double t_max = 0.3;
  int N = 101;
  double x_max = 2.;
  double dx = x_max / (N - 1.);

  // Initial arrays
  double u0[N];  // for Leapfrog method
  u0[0] = 0.;
  u0[N - 1] = 0.;
  for (int i = 1; i < N - 1; i++) {
    u0[i] = 3.0 * sin(M_PI * i * dx);
  }
  // Discretize the mesh
  double epsilon = 1.;              // wave speed [dimensionless]
  double dt = beta * dx / epsilon;  // time step size
  int M = (int)(t_max / dt) + 1;    // number of time points

  TFile* f = new TFile("test.root");
  TGraph2D* ug = new TGraph2D();

  for (int i = 0; i < N; i++) {  // write out initial condition
    double x = i*dx;
    ug->AddPoint(x,0,u0[i]);
  }

  // Create array to hold values at j + 1
  double u[N];

  // Perform leapfrog method
  for (int j = 1; j < M; j++) {  // Loop through time
    // Apply boundary conditions
    u[0] = 0;
    u[N - 1] = 0;
    ug->AddPoint(0,j*dt,0);
    for (int i = 1; i < N - 1; i++) {  // Loop through space (inner points only)
      u[i] =
          u0[i] -
          0.25 * beta * (u0[i + 1] * u0[i + 1] - u0[i - 1] * u0[i - 1]) +
          0.125 * beta * beta *
              ((u0[i + 1] + u0[i]) * (u0[i + 1] * u0[i + 1] - u0[i] * u0[i]) -
               (u0[i] + u0[i - 1]) * (u0[i] * u0[i] - u0[i - 1] * u0[i - 1]));
      ug->AddPoint(i*dx,j*dt,u[i]);
    }
    ug->AddPoint(x_max,j*dt,u[N-1]);
    memcpy(u0, u, sizeof(u));  // Let u0 = u to get next time step
  }
  TCanvas* c = new TCanvas("c","Lax-Wendroff Method");
  ug->Draw("surf3");
  c->Write();
  ug->Write();
}

int main() {
  burger_root();
  return 0;

}
