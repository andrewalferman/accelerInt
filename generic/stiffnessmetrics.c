/**
 * \file
 * \Calculates the stiffness ratio, stiffness indicator, and chemical explosive
 * \mode
 *
 * \author Andrew Alferman
 * \date 12/14/2017
 *
 */

 //#include "timer.h"
 #include "header.h"
 #include "jacob.h"
 #include "stiffnessmetrics.h"

 /* DGEEV prototype */
 extern void dgeev( char* jobvl, char* jobvr, int* n, double* a,
                 int* lda, double* wr, double* wi, double* vl, int* ldvl,
                 double* vr, int* ldvr, double* work, int* lwork, int* info );

 /* Parameters */
 #define N NSP
 #define LDA N
 #define LDVL N
 #define LDVR N
 // Need a better way of sending this the N2 position
 #define N2POS 48

 void calculatemetrics(double *y_local, double stiffratio,
                      double stiffindicator, double CEM); {
   // Rearrange the solution vector for pyJac
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
   // Calculate the stiffness metrics
   // Get the Jacobian
   double jac[NSP*NSP];
   //printf("%.15e\n", pr_global[tid]);
   double pr_stiffcalc;
   if (pr_global[tid] >= 1000.0)
   {
     pr_stiffcalc = pr_global[tid] / 101325.0;
   }
   else
   {
     pr_stiffcalc = pr_global[tid];
   }
   eval_jacob(t_end, pr_stiffcalc, re_local, jac);
   // Get the Hermitian
   double hermitian[NSP*NSP];
   for (int i = 0; i < NSP; i++) {
     for (int j = 0; j < NSP; j++) {
       //jacobian[i][j] = jac[i * NSP + j];
       hermitian[i * NSP + j] = 0.5 * (jac[i * NSP + j] + jac[j * NSP + i]);
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
   for (int i = 0; i < NSP; i++) {
     if (abs(wr[i]) < minjaceig && abs(wr[i]) > 2.22045e-16) {
       minjaceig = abs(wr[i]);
     }
     if (wr[i] > CEM) {
       CEM = wr[i];
     }
     if (abs(wr[i]) > maxjaceig) {
       maxjaceig = abs(wr[i]);
     }
     if (xr[i] < minhereig) {
       minhereig = xr[i];
     }
     if (xr[i] > maxhereig) {
       maxhereig = xr[i];
     }
   }

   stiffratio = maxjaceig / minjaceig;
   stiffindicator = 0.5 * (minhereig + maxhereig);

}
