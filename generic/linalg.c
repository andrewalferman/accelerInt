/** Include common variables and libraries. */
#include "head.h"

/** File containing linear algebra functions.
 * \file linalg.c
 *
 * \author Kyle E. Niemeyer
 * \date 08/02/2011
 *
 * Contains linear algebra (eigensolve, inversion) functions, using LAPACK.
 */

/** Interface function to LAPACK eigensolve subroutine.
 *
 * Performs eigendecomposition of square matrix. Uses LAPACK DGEEV subroutine.
 * Returns eigenvalues and right eigenvectors.
 * 
 * \param[in]     n  	  order of matrix
 * \param[in]     A     the matrix, size n*n
 * \param[out]    evalr the array of eigenvalue real parts
 * \param[out]    evali the array of eigenvalue imaginary parts
 * \param[out]    levec	the array of left eigenvectors, size n*n
 * \param[out]    revec the array of right eigenvectors, size n*n
 * \return        info  success/fail integer flag
 */
uint eigen ( uint n, double *A, double *evalr, double *evali, 
             double *levec, double *revec ) {
  
  // Prototype for LAPACK DGEEV (Fortran) function
  extern void dgeev_ ( char *jobvl, char *jobvr, uint *n, double *A, 
                       uint *lda, double *wr, double *wi, double *vl, 
                       uint *ldvl, double *vr, uint *ldvr, double *work, 
                       uint *lwork, uint *info );
  
  char jobvl = 'V'; // Telling DGEEV to not return left eigenvectors
  char jobvr = 'V'; // Telling DGEEV to return right eigenvectors
  
  uint lwork = 5 * n;   // size of work array
  double work[5 * n];   // real work array
  
  uint info;        // Output success/fail flag from dgeev
  
  dgeev_ ( &jobvl, &jobvr, &n, A, &n, evalr, evali, levec, &n, revec, &n, 
           work, &lwork, &info ); 
  
  return info;
}

/** Interface function to LAPACK eigensolve subroutine.
 *
 * Performs eigendecomposition of square matrix. Uses LAPACK DGEEV subroutine.
 * Returns eigenvalues and right eigenvectors.
 * 
 * \param[in]     n  	  order of matrix
 * \param[in]     A     the matrix, size n*n
 * \param[out]    evalr the array of eigenvalue real parts
 * \param[out]    evali the array of eigenvalue imaginary parts
 * \param[out]    levec	the array of left eigenvectors, size n*n
 * \param[out]    revec the array of right eigenvectors, size n*n
 * \return        info  success/fail integer flag
 */
uint eigenbal ( uint n, double *A, double *evalr, double *evali, 
                double *levec, double *revec ) {
  
  // Prototype for LAPACK DGEEVX (Fortran) function
  extern void dgeevx_ ( char *bal, char *jobvl, char *jobvr, char *sens,
                        uint *n, double *A, uint *lda, double *wr, 
                        double *wi, double *vl, uint *ldvl, double *vr, 
                        uint *ldvr, uint *ilo, uint *ihi, double *scale,
                        double *abnrm, double *rconde, double *rcondv,
                        double *work, uint *lwork, uint *iwork, uint *info );
  
  char bal = 'B';   // Telling DGEEVX to permute and scale matrix
  char sens = 'B';  // Telling DGEEVX to compute reciprocal cond numbers for eigenvalues and right eigenvectors
  
  char jobvl = 'V'; // Telling DGEEV to not return left eigenvectors
  char jobvr = 'V'; // Telling DGEEV to return right eigenvectors
  
  uint lwork = n * (n + 8);   // size of work array
  double work[lwork];         // local work array
  double scale[n];    // details of DGEEVX permutation and scaling factors
  double abnrm;       // one-norm of balanced matrix
  double rconde[n];   // array of reciprocal condition numbers of left eigenvectors
  double rcondv[n];   // array of reciprocal condition numbers of right eigenvectors
  uint iwork[2*n-2];  // local integer work array
  
  uint ilo, ihi;      // used for matrix balancing
  uint info;          // Output success/fail flag from dgeev
  
  dgeevx_ ( &bal, &jobvl, &jobvr, &sens, &n, A, &n, evalr, evali, levec, &n,
            revec, &n, &ilo, &ihi, scale, &abnrm, rconde, rcondv, work, 
            &lwork, iwork, &info ); 
  
  return info;
}

/** Interface function to LAPACK matrix inversion subroutine.
 *
 * Performs inversion of square matrix. Uses LAPACK subroutines DGETRF and DGETRI.
 * 
 * \param[in]     n     order of matrix
 * \param[in]     A     the input matrix, size n*n
 * \param[out]    B     the output matrix, size n*n
 * \return        info  success/fail integer flag
 */
uint invers ( uint n, double *A, double *B ) {
  
  // Prototype for LAPACK DGETRF (Fortran) function
  extern void dgetrf_ ( uint *m, uint *n, double *A, uint *lda, uint *ipiv, uint *info );
  
  // Prototype for LAPACK DGETRI (Fortran) function
  extern void dgetri_ ( uint *n, double *A, uint *lda, uint *ipiv, double *work, 
                        uint *lwork, uint *info );
  
  uint lwork = 2 * n; // size of work array
  double work[2 * n]; // local work array
  
  uint ipiv[n];   // array holding pivot indices
  uint info;      // Output success/fail flag from dgeev
  
  // copy contents of A to B so that A is not overridden
  memcpy ( B, A, sizeof(double) * n * n );
  
  // first call dgetrf for LU factorization
  dgetrf_ ( &n, &n, B, &n, ipiv, &info );
  
  if ( info != 0 ) {
    printf( "Error in dgetrf, info = %d\n", info );
    exit (1);
  }
  
  // now call dgetri for inversion
  dgetri_ ( &n, B, &n, ipiv, work, &lwork, &info ); 
  
  return info;
}