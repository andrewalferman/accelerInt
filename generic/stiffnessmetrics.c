/**
 * \file
 * \Calculates the stiffness ratio, stiffness indicator, chemical explosive
 * \mode, and stiffness per the computational singular perturbation method
 *
 * \author Andrew Alferman
 * \date 12/14/2017
 *
 */

 //#include <stdio.h>
 //#include <float.h>

 //#include "timer.h"
#include "header.h"
#include "jacob.h"
//#include "stiffnessmetrics.h"
#include "head.h"


 /* DGEEV prototype */
 extern void dgeev( char* jobvl, char* jobvr, int* n, double* a,
                 int* lda, double* wr, double* wi, double* vl, int* ldvl,
                 double* vr, int* ldvr, double* work, int* lwork, int* info );
 uint get_slow_projector ( Real tim, Real * y, Real * Qs, Real * taum1, Real * Rc, double pr_local  );

 /* Parameters */
 #define N NSP
 #define LDA N
 #define LDVL N
 #define LDVR N
 // Need a better way of sending this the N2 position
 #ifdef CHEM_UTILS_HEAD
   #if (NSP == 53)
    #define N2POS 48
   #elif (NSP == 13)
    #define N2POS 11
   #endif
 #endif

 void calculatemetrics(double* y_local, double pr_local, double* stiffratio,
                      double* stiffindicator, double* CEM, double* CSP, int* M,
                      const double t, const double t_end)
  {
   // Rearrange the solution vector for pyJac
   #ifdef CHEM_UTILS_HEAD
     double re_local[NSP];
     double nmf;
     for (int i = 0; i < NSP; i++)
     {
       re_local[i] = y_local[i];
       if (i == N2POS)
       {
         double nmf = y_local[i];
       }
     }
     re_local[N2POS] = re_local[NSP];
     re_local[NSP] = nmf;
   #else
   double re_local[NSP];
   for (int i = 0; i < NSP; i++)
   {
     re_local[i] = y_local[i];
   }
   #endif
   for (int i = 0; i < NSP; i++) {
	printf("%.15e,",re_local[i]);
   }
   printf("\n");
   // Calculate the stiffness metrics
   // Get the Jacobian
   double jac[NSP*NSP];
   double pr_stiffcalc;
   if (pr_local >= 1000.0)
   {
     pr_stiffcalc = pr_local / 101325.0;
   }
   else
   {
     pr_stiffcalc = pr_local;
   }
   eval_jacob(t, pr_local, re_local, jac);
   // Get the Hermitian
   double hermitian[NSP*NSP];
   for (int i = 0; i < NSP; i++) {
     for (int j = 0; j < NSP; j++) {
       hermitian[i * NSP + j] = (double) 0.5 * (jac[i * NSP + j] + jac[j * NSP + i]);
     }
   }
   // Get the inverse of the diagonals of the Jacobian matrix
   double diagonals[NSP];
   for (int i = 0; i < NSP; i++) {
     if (jac[i * NSP + i] > DBL_MIN) {
       diagonals[i] = (double) 1.0 / jac[i * NSP + i];
     }
     else if ((double) -1.0 * jac[i * NSP + i] > DBL_MIN) {
       diagonals[i] = (double) -1.0 / jac[i * NSP + i];
     }
     else {
       diagonals[i] = -1.0;
     }
   }
   // Get the eigenvalues of both matrices
   // Made 2 sets of variables in case dgeev messes them up
   int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info, lwork;
   int nh = N, ldah = LDA, ldvlh = LDVL, ldvrh = LDVR, infoh, lworkh;
   double wkopt;
   double* work;
   double wkopth;
   double* workh;
   /* Local arrays */
   double wr[N], wi[N], vl[LDVL*N], vr[LDVR*N];
   double xr[N], xi[N], ul[LDVL*N], ur[LDVR*N];

   printf("Got here\n");
   /* Query and allocate the optimal workspace */
   // First, the Jacobian
   lwork = -1;
   dgeev( "Vectors", "Vectors", &n, jac, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
    &wkopt, &lwork, &info );
   lwork = (int)wkopt;
   work = (double*)malloc( lwork*sizeof(double) );
   /* Solve eigenproblem */
   dgeev( "Vectors", "Vectors", &n, jac, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
    work, &lwork, &info );
   /* Check for convergence */
   if( info > 0 ) {
           printf( "The algorithm failed to compute eigenvalues.\n" );
           exit( 1 );
   }
   // Next, the Hermitian
   printf("Also got here.\n")
   lworkh = -1;
   dgeev( "Vectors", "Vectors", &nh, hermitian, &ldah, xr, xi, ul, &ldvlh, ur, &ldvrh,
    &wkopth, &lworkh, &infoh );
   lworkh = (int)wkopth;
   workh = (double*)malloc( lworkh*sizeof(double) );
   /* Solve eigenproblem */
   dgeev( "Vectors", "Vectors", &nh, hermitian, &ldah, xr, xi, ul, &ldvlh, ur, &ldvrh,
    workh, &lworkh, &infoh );
   /* Check for convergence */
   if( info > 0 ) {
           printf( "The algorithm failed to compute eigenvalues.\n" );
           exit( 1 );
   }

   // Get the minimum and maximum values of the eigenvalues
   double minjaceig = 1.0e10;
   double maxjaceig = 0.0;
   double minhereig = 1.0e10;
   double maxhereig = 0.0;
   // double mindiag = 1.0e10;
   // double maxdiag = 0.0;
   double timescale;
   double minfasttimescale = 1e10;
   // double minslowtimescale = 1e10;

   // CSP slow-manifold projector matrix
   Real Qs[NN*NN];

   // CSP radical correction tensor matrix
   Real Rc[NN*NN];

   // time scale of fastest slow mode (controlling time scale)
   Real tau;
   Real * ptau = &tau; // pointer to tau

   // get slow-manifold projector, driving time scale, and radical correction
   (*M) = get_slow_projector ( t, y_local, Qs, ptau, Rc , pr_local );

   /** Time step factor mu.
    * Time step divided by controlling time scale.
    */
   const Real mu = 0.005;

   // time step size (controlled by fastest slow mode)
   Real h = mu * tau;

   double eigenvalue;
   double deltaT = t_end - t;
   for (int i = 0; i < NSP; i++) {
     if (wr[i] < 0) {
       eigenvalue = wr[i] * (double) -1.0;
     }
     else {
       eigenvalue = wr[i];
     }
     // Only doing these checks if we have a nonzero eigenvalue
     if (eigenvalue > DBL_MIN) {
       timescale = (double) 1.0 / eigenvalue;
       if (timescale < minfasttimescale) {
         minfasttimescale = timescale;
       }
       // if ((timescale < minslowtimescale) && (timescale > deltaT)) {
       //   minslowtimescale = timescale;
       // }
     }
     if ((eigenvalue < minjaceig) && (eigenvalue > DBL_MIN)) {
       minjaceig = eigenvalue;
     }
     if (wr[i] > (*CEM)) {
       (*CEM) = wr[i];
     }
     if (eigenvalue > maxjaceig) {
       maxjaceig = eigenvalue;
     }
     if (xr[i] < minhereig) {
       minhereig = xr[i];
     }
     if (xr[i] > maxhereig) {
       maxhereig = xr[i];
     }
     // if ((diagonals[i] < mindiag) && (diagonals[i] != -1.0)) {
     //   mindiag = diagonals[i];
     // }
     // if (diagonals[i] > maxdiag) {
     //   maxdiag = diagonals[i];
     // }

   }

   (*CSP) = minfasttimescale / tau;
   (*stiffratio) = maxjaceig / minjaceig;
   (*stiffindicator) = (double) 0.5 * (minhereig + maxhereig);

}
