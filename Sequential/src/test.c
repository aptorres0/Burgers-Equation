#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "burger.h"

int main(int argc, char* argv[]) {
  if (argc < 2) {
    printf("Usage: %s beta tmax", argv[0]);
    return 1;
  }

  double beta = atof(argv[1]);
  double t_max = atof(argv[2]);
  // Define the mesh and input parameters
  int N = 101;
  double x_max = 2.;
  double dx = x_max / (N - 1.);
  double epsilon = 1.0;
  double dt = beta * dx / epsilon;
  int M = (int)(t_max / dt) + 1;

  char* tmp2 = "../results/lw-";
  char* fName_LW = (char*)malloc(sizeof(tmp2) + 50 * sizeof(char) +
                                 sizeof(argv[1]) + sizeof(argv[2]));
  sprintf(fName_LW, "%sb%s-t%s.txt", tmp2, argv[1], argv[2]);
  char* lw_info_name = (char*)malloc(sizeof(tmp2) + 50 * sizeof(char) +
                                     sizeof(argv[1]) + sizeof(argv[2]));
  sprintf(lw_info_name, "%sb%s-t%s_info.txt", tmp2, argv[1], argv[2]);

  char* tmp = "../results/lf-";
  char* fName_LF = (char*)malloc(sizeof(tmp) + 50 * sizeof(char) +
                                 sizeof(argv[1]) + sizeof(argv[2]));
  sprintf(fName_LF, "%sb%s-t%s.txt", tmp, argv[1], argv[2]);

  char* lf_info_name = (char*)malloc(sizeof(tmp) + 50 * sizeof(char) +
                                     sizeof(argv[1]) + sizeof(argv[2]));
  sprintf(lf_info_name, "%sb%s-t%s_info.txt", tmp, argv[1], argv[2]);

  FILE* lf_info = fopen(lf_info_name, "w");
  if (lf_info == NULL) {
    printf("Error creating file %s!\n", lf_info_name);
    exit(1);
  }
  FILE* lw_info = fopen(lf_info_name, "w");
  if (lw_info == NULL) {
    printf("Error creating file %s!\n", lw_info_name);
    exit(1);
  }
  // Initial array
  double u0[N], u1[N];
  u0[0] = 0.;
  u0[N - 1] = 0.;
  for (int i = 1; i < N - 1; i++) {
    u0[i] = 3.0 * sin(M_PI * i * dx);
  }
  memcpy(u1, u0, sizeof(u0));
  // memcpy(u2, u0, sizeof(u0));

  LeapFrogMethod(u0, N, x_max, t_max, beta, fName_LF);
  LaxWendroffMethod(u1, N, x_max, t_max, beta, fName_LW);
  printf(
      "N = %d\tx_max = %f\tdx = %f\nM = %d\tt_max = %f\tdt = %f\nbeta = "
      "%f\nepsilon = "
      "%f\nOutput files:\n%s\n%s\n",
      N, x_max, dx, M, t_max, dt, beta, epsilon, fName_LW, fName_LF);
  fprintf(lf_info,
          "N = %d\nx_max = %f\ndx = %f\nM = %d\nt_max = %f\ndt = %f\nbeta = "
          "%f\nepsilon = "
          "%f\nOutput files:\n%s\n%s\n",
          N, x_max, dx, M, t_max, dt, beta, epsilon, fName_LW, fName_LF);
  fprintf(lw_info,
          "N = %d\nx_max = %f\ndx = %f\nM = %d\nt_max = %f\ndt = %f\nbeta = "
          "%f\nepsilon = "
          "%f\nOutput files:\n%s\n%s\n",
          N, x_max, dx, M, t_max, dt, beta, epsilon, fName_LW, fName_LF);

  return 0;
}
