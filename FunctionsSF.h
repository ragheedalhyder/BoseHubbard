//
//  Header.h
//  MISF_Final
//
//  Created by Ragheed on 29.06.23.
//

#ifndef Header_h
#define Header_h
#include <iostream>
#include <fstream>
#include <complex>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <typeinfo>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
// #include <lapacke.h>
#include <vector>
#include <string>
#include <ctime>
#include <armadillo>
#include <time.h>


using namespace std;
using namespace Eigen;
using namespace arma;

#define pi 3.14159265358979323846
#define N 10 // Number of bosons in bosonic fock states
#define L 9 // Number of lattice sites in the X and Y directions (Lx = Ly = L)
#define M 10
#define epsilon 1e-7


int delta(int x, int y);
///
///
double cn(vec cns, int n);
///
///
double Psi0(vec cns);
///
///
double JkU(double dJU, double x);
///
///
double epsI(double kx, double ky);
///
double U(int kx, int ky, int lambda, int qx, int qy, int lambda1, double uks[N][M][M][N], double KXs[M], double KYs[M], double n0 );
///
///
double V(int kx, int ky, int lambda, int qx, int qy, int lambda1, double vks[N][M][M][N], double KXs[M], double KYs[M], double n0 );
///
///
double W(int kx, int ky, int lambda, int qx, int qy, int lambda1, double uks[N][M][M][N], double vks[N][M][M][N], double KXs[M], double KYs[M], double n0 );
///
///
double Nk(int kx, int ky, int lambda, vec cns, double uks[N][M][M][N], double vks[N][M][M][N], double KXs[M], double KYs[M]);
///
///

cx_vec SigmaPolaron (double Epol, double KXs[M], double KYs[M], double dkxs[M], double dkys[M], double uks[N][M][M][N], double vks[N][M][M][N], double omegaklambda[N][M][M], double dJU, double n0, double UIB, int cutoff);

cx_vec FixedPoint(double sig1, double sig2, double Epol, double dEpol, double KXs[M], double KYs[M],  double dkxs[M], double dkys[M], double uks[N][M][M][N], double vks[N][M][M][N], double omegaklambda[N][M][M], double dJU, double n0, double UIB, int cutoff);


cx_vec LimitFinder(double sig1, double Epol, double dEpol, double KXs[M], double KYs[M],  double dkxs[M], double dkys[M], double uks[N][M][M][N], double vks[N][M][M][N], double omegaklambda[N][M][M], double dJU, double n0, double UIB, int cutoff);


#endif /* Header_h */
