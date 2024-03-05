/**
 * burger.h
 *
 * Author: Alexander Paul Torres
 * Date: 30 JAN 23
 * Class: PHYS 7411 - Computational Physics I
 *
 * Description:
 * Header file to accompany burger.c
 *
 */
#ifndef _BURGER_H
#define _BURGER_H

/**
 * Function: LeapFrogMethod
 * ------------------------
 * computes the inviscid Burgers' equation (First-order Non-linear Convection)
 * using a forward difference time derivative and a central difference
 * spatial derivative scheme. The algorithm is accurate to second order
 * in space and first order in time.
 *
 * inputs:
 *      u0      - initial array of wave height at time t = 0 of length N
 *      N       - number of grid points in space
 *      x_max   - upper boundary for space [dimensionless]
 *      t_max   - upper boundary for time [dimensionless]
 *      beta    - CFL number [dimensionless]
 *      fname   - Output file name
 *
 * returns: writes the output to a file named fname
 *
 */
void LeapFrogMethod(double u0[], int N, double x_max, double t_max, double beta,
                    char* fname);

/**
 * Function: LaxWendroffMethod
 * ------------------------
 * computes the inviscid Burgers' equation (First-order Non-linear Convection)
 * The algorithm is accurate to second order in space and time.
 *
 * inputs:
 *      u0      - initial array of wave height at time t = 0 of length N
 *      N       - number of grid points in space
 *      x_max   - upper boundary for space [dimensionless]
 *      t_max   - upper boundary for time [dimensionless]
 *      beta    - CFL number [dimensionless]
 *      fname   - output file name
 *
 * returns: writes the output to a file named fname
 *
 */
void LaxWendroffMethod(double u0[], int N, double x_max, double t_max,
                       double beta, char* fname);

// struct WaveAmplitude;

#endif