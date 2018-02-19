/** File containing linear algebra and all CSP-related functions for CSP model problem project.
 * \file csp_analysis.c
 *
 * \author Kyle E. Niemeyer
 * \date 08/02/2011
 *
 * Contains functions that calculate CSP vectors and timescales,
 * CSP projection matrix and radical correction tensor.
 */

/** Include common code. */
#include "head.h"
#include "header.h"
//#include "jacob.h"

/** Relative CSP error tolerance */
const Real eps_r = 1.0e-3;

/** Absolute CSP error tolerance */
const Real eps_a = 1.0e-3;

/** Define sign function */
#define sign(x) (copysign(ONE, x))

/** Function prototypes. */
void dydt (const double t, const double pres, const double * y, double * dy);
//void dydt ( Real t, Real * y, Real * Q, int qflag, Real * dydt );
void insertion_sort ( uint n, Real * vals, uint * order );
uint eigen ( uint n, double * A, double * evalr, double * evali,
             double * lev, double * rev );
uint eigenbal ( uint n, double * A, double * evalr, double * evali,
                double * lev, double * rev );
uint invers ( uint n, double * A, double * B );
uint get_fast_modes ( Real tim, Real * y, Real * tau, Real * a_csp, Real * b_csp, double pr_local);
void get_csp_vectors ( Real tim, Real * y, Real * tau, Real * a_csp, Real * b_csp, double pr_local);

////////////////////////////////////////////////////////////////////////
/** Function that applies radical correction to data array
 *
 * Performs matrix-vector multiplication of radical correction tensor
 * with data array.
 *
 * \param[in]   tim time at new time step (sec)
 * \param[in]   pres pressure at the time step
 * \param[in]   y   array of data already integrated, size NN
 * \param[in]   Rc  Radical correction tensor, size NN*NN
 * \param[out]  g   array holding radical corrections, size NN
 */
//void radical_correction ( Real tim, Real * y, Real * Rc, Real * g ) {
void radical_correction ( const double tim, const double pres,
                          const double * y, Real * Rc, Real * g ) {

  // call derivative, and apply radical correction tensor on derivative vector
  double dy[NN];

  dydt ( tim, pres, y, dy );

  Real ydotn[NN];

  for ( uint i = 0; i < NN; ++i ) {

    Real sum_row = ZERO;
    for ( uint j = 0; j < NN; ++j ) {
      sum_row += Rc[i + (NN * j)] * dy[j];
    } // end j loop

    ydotn[i] = sum_row;
  } // end i loop

  // now replace ydot with ydotn
  for ( uint i = 0; i < NN; ++i ) {
    g[i] = ydotn[i];
  } // end i loop

}

////////////////////////////////////////////////////////////////////////

/** Function that performs CSP analysis and returns vectors.
 *
 *
 *
 * \param[in]   tim   current time (s)
 * \param[in]   y     array of values at current time, size NN
 * \param[out]  Qs    slow-manifold projector matrix, size NN*NN
 * \param[out]  taum1 time scale of fastest slow mode (sec)
 * \param[out]  Rc    radical correction tensor, size NN*NN
 * \return      M     number of slow modes
 */
uint get_slow_projector ( Real tim, Real * y, Real * Qs, Real * taum1, Real * Rc, double pr_local ) {

  // CSP mode timescales (s)
  // Positive values indicate explosive modes (never in M exhausted modes)
  Real tau[NN];

  // CSP vectors and covectors
  Real a_csp[NN*NN]; // Array with CSP vectors
  Real b_csp[NN*NN]; // Array with CSP covectors

  // find number of exhausted modes (M) and CSP vectors
  uint M = get_fast_modes ( tim, y, tau, a_csp, b_csp, pr_local);

  // calculate slow-manifold projector
  for ( uint j = 0; j < NN; ++j ) {

    for ( uint i = 0; i < NN; ++i ) {

      Rc[i + (NN * j)] = ZERO;

      // Qs starts as identity matrix
      if ( i == j ) {
        Qs[i + (NN * j)] = ONE;
      } else {
        Qs[i + (NN * j)] = ZERO;
      }

      // ensure at least 1 exhausted mode
      if ( M > 0 ) {

        Real sum_qs = ZERO;
        Real sum_rc = ZERO;

        /* original summation
        // sum over slow modes
        for ( uint r = 0; r < M; ++r ) {
          sum_qs += a_csp[i + (NN * r)] * b_csp[j + (NN * r)];
          //sum_rc += a_csp[i + (NN * r)] * b_csp[j + (NN * r)] * fabs( tau[r] );
          sum_rc += a_csp[i + (NN * r)] * b_csp[j + (NN * r)] * tau[r];
        } // end r loop
        */

        // Kahan summation
        Real c_qs = ZERO;
        Real c_rc = ZERO;
        // sum over slow modes
        for ( uint r = 0; r < M; ++r ) {
          Real y_sum = ( a_csp[i + (NN * r)] * b_csp[j + (NN * r)] ) - c_qs;
          Real t = sum_qs + y_sum;
          c_qs = ( t - sum_qs ) - y_sum;
          sum_qs = t;

          // don't use absolute value of time scales (want negatives)
          y_sum = ( a_csp[i + (NN * r)] * b_csp[j + (NN * r)] * tau[r] ) - c_rc;
          t = sum_rc + y_sum;
          c_rc = ( t - sum_rc ) - y_sum;
          sum_rc = t;
        } // end r loop

        // Qs = Id - sum_r^M a_r*b_r
        Qs[i + (NN * j)] -= sum_qs;

        // Rc = sum_r^M a_r*tau_r*b_r
        Rc[i + (NN * j)] = sum_rc;

      } // end if

    } // end j loop

  } // end i loop

  // return time scale of fastest slow mode
  *taum1 = fabs( tau[M] );

  // return number of slow modes
  return M;
}

////////////////////////////////////////////////////////////////////////

/** Function that performs CSP analysis and returns vectors.
 *
 * Performs eigendecomposition of Jacobian to get the eigenvalues (inverse of
 * mode timescales), CSP vectors (from right eigenvectors) and CSP covectors
 * (from left eigenvectors). Uses LAPACK Fortran subroutines, also
 * sorts (based on eigenvalue magnitude) and normalizes eigenvectors.
 *
 * Explosive modes (positive eigenvalues) will always be retained, since sorted
 * in ascending order (and all others are negative). Complex eigenvalues (with
 * associated complex eigenvectors) are handled also by taking the sum and
 * difference of the real and imaginary parts of the eigenvectors.
 *
 * \param[in]   tim     current time (s)
 * \param[in]   y       array of values at current time, size NN
 * \param[out]  tau     array of CSP mode timescales (s), size NN
 * \param[out]  a_csp   CSP vectors, size NN*NN
 * \param[out]  b_csp   CSP covectors, size NN*NN
 */
void get_csp_vectors ( Real tim, Real * y, Real * tau, Real * a_csp, Real * b_csp, double pr_local ) {

  // array holding Jacobian
  Real jac[NN*NN];

  if (pr_local >= 1000.0)
  {
    pr_local = pr_local / 101325.0;
  }

  // calculate Jacobian for current time step
  eval_jacob ( tim, pr_local, y, jac );

  // arrays of eigenvalues and eigenvectors
  Real evalr[NN];     // Array with real parts of eigenvalues
  Real evali[NN];     // Array with imaginary parts of eigenvalues
  Real evecr[NN*NN];  // Local array with right eigenvectors
  Real evecl[NN*NN];  // Local array with left eigenvectors

  // calculate eigenvalues and eigenvectors
  eigenbal ( NN, jac, evalr, evali, evecl, evecr );

  // sort according to eigenvalues
  uint order[NN];    // Array holding indices in sorted order

  // ascending insertion sort.
  // Largest magnitude negative eigenvalues first, followed by smaller magnitude
  // negative eigenvalues. Any explosive modes (positive eigenvalues) will be at
  // the end, regardless of magnitude.
  insertion_sort ( NN, evalr, order );

  for ( uint i = 0; i < NN; ++i ) {
    tau[i] = ONE / evalr[ order[i] ]; // time scales, inverse of eigenvalues

    for ( uint j = 0; j < NN; ++j ) {
      a_csp[j + NN * i] = evecr[j + NN * order[i]]; // CSP vectors, right eigenvectors
      b_csp[j + NN * i] = evecl[j + NN * order[i]]; // CSP covectors, left eigenvectors
    }
  }

  // eliminate complex components of eigenvectors if complex eigenvalues,
  // and normalize dot products (so that bi*aj = delta_ij).
  uint flag = 1;
  for ( uint i = 0; i < NN; ++i ) {

    // check if imaginary part of eigenvalue (skip 2nd eigenvalue of complex pair)
    if ( ( fabs(evali[i]) > SMALL ) && ( flag == 1 ) ) {
      // complex eigenvalue

      // normalize
      // needs to be able to handle complex eigenvectors
      double sum_r = 0.0; // real part
      double sum_i = 0.0; // imaginary part

      for ( uint j = 0; j < NN; ++j ) {
        uint ir = j + NN * i;        // location of real part
        uint ii = j + NN * (i + 1);  // location of imaginary part

        //sum_r += ( vr[ir] * vl[ir] ) - ( vr[ii] * vl[ii] );
        //sum_i += ( vr[ii] * vl[ir] ) + ( vr[ir] * vl[ii] );
        // need to treat left eigenvector as complex conjugate (so take negative of imag part)
        sum_r += ( a_csp[ir] * b_csp[ir] ) + ( a_csp[ii] * b_csp[ii] );
        sum_i += ( a_csp[ii] * b_csp[ir] ) - ( a_csp[ir] * b_csp[ii] );
      } // end j loop

      Real sum2 = ( sum_r * sum_r ) + ( sum_i * sum_i );

      // ensure sum is not zero
      if ( fabs(sum2) > SMALL ) {

        for ( uint j = 0; j < NN; ++j ) {
          uint ir = j + NN * i;
          uint ii = j + NN * (i + 1);

          // normalize a, and set a1=real, a2=imag
          Real a_old = a_csp[ir];
          a_csp[ir] = ( (a_old * sum_r) + (a_csp[ii] * sum_i) ) / sum2;
          a_csp[ii] = ( (a_csp[ii] * sum_r) - (a_old * sum_i) ) / sum2;

          // set b1=2*real, b2=-2*imag
          b_csp[ir] = TWO * b_csp[ir];
          //vl[ii] = -TWO * b_csp[ii];
          b_csp[ii] = TWO * b_csp[ii];
        } // end j loop

      } // end if

      // skip next (conjugate of current)
      flag = 2;

    } else if ( flag == 2 ) {
      // do nothing, this is second of complex pair
      flag = 1;

    } else {
      // real eigenvalue
      flag = 1;

      /* original summation
      Real sum = ZERO;
      for ( uint j = 0; j < NN; ++j ) {
        sum += a_csp[j + NN * i] * b_csp[j + NN * i];
      } // end j loop
      */

      // Kahan summation
      Real sum = ZERO;
      Real c = ZERO;
      for ( uint j = 0; j < NN; ++j ) {
        Real y = ( a_csp[j + NN * i] * b_csp[j + NN * i] ) - c;
        Real t = sum + y;
        c = ( t - sum ) - y;
        sum = t;
      } // end j loop

      // ensure dot product is not zero
      if ( fabs(sum) > SMALL ) {

        for ( uint j = 0; j < NN; ++j ) {
          // just normalize a
          a_csp[j + NN * i] /= sum;
        } // end j loop

      } // end if

    } // end if

  } // end i loop

// #ifdef ERROR
//
//   // ensure a and b are inverses
//   for ( uint i = 0; i < NN; ++i ) {
//
//     Real I[NN];
//     for ( uint j = 0; j < NN; ++j ) {
//       I[j] = 0;
//
//       for ( uint k = 0; k < NN; ++k ) {
//         I[j] += a_csp[k + NN*j]*b_csp[k + NN*i];
//       }
//
//       if ( i == j ) {
//         if ( fabs( I[j] - ONE ) > 1.0e-14 ) {
//           printf("CPS vectors not orthogonal\n");
//         }
//       } else {
//         if ( fabs( I[j] ) > 1.0e-14 ) {
//           printf("CPS vectors not orthogonal\n");
//         }
//       }
//       //printf("%17.10le ", I[j]);
//     } // end j loop
//
//     //printf("\n");
//   } // end i loop
//
//
//   // check that new CSP vectors and covectors recover standard base vectors e_i
//   for ( uint ei = 0; ei < NN; ++ei ) {
//     // fill e_i
//     Real e[NN];
//     for ( uint i = 0; i < NN; ++i ) {
//       if ( ei == i ) {
//         e[i] = ONE;
//       } else {
//         e[i] = ZERO;
//       }
//     }
//
//     for ( uint i = 0; i < NN; ++i ) {
//       // ith component of e
//       double e_comp = ZERO;
//       for ( uint j = 0; j < NN; ++j ) {
//         // (b.e_i)*a
//
//         // b.e_i
//         double e_sum = ZERO;
//         for ( uint k = 0; k < NN; ++k ) {
//           e_sum += b_csp[k + NN * j] * e[k];
//         }
//
//         e_sum *= a_csp[i + NN * j];
//         e_comp += e_sum;
//       }
//
//       // if reconstructed basis not very close to original basis
//       if ( fabs( e[i] - e_comp ) > 1.0e-10 ) {
//         printf("Error recreating standard basis vectors.\n");
//         printf("e_%d, component %d, e_comp: %17.10le \t error: %17.10le\n", ei+1, i, e_comp, fabs( e[i] - e_comp ));
//       }
//     }
//   }
// #endif

} // end get_csp_vectors

////////////////////////////////////////////////////////////////////////

/** Function that returns the number of exhausted modes.
 *
 *
 * \param[in]   tim     current time (s)
 * \param[in]   y       array of values at current time, size NN
 * \param[out]  tau     array of CSP mode timescales (s), size NN
 * \param[out]  a_csp   CSP vectors, size NN*NN
 * \param[out]  b_csp   CSP covectors, size NN*NN
 * \return      M       number of exhausted (fast) modes
 */
uint get_fast_modes ( Real tim, Real * y, Real * tau, Real * a_csp, Real * b_csp, double pr_local) {

  // perform CSP analysis to get timescales, vectors, covectors
  get_csp_vectors ( tim, y, tau, a_csp, b_csp, pr_local );

  // now need to find M exhausted modes
  // first calculate f^i
  Real f_csp[NN];   // f^i is b* operated on g
  Real g_csp[NN];   // array of derivatives (g in CSP)

  Real * Q; // pointer to unused projector array
  int qflag = 0; // flag telling dydt to not use projector

  // call derivative function
  dydt ( tim, pr_local, y, g_csp );

  for ( uint i = 0; i < NN; ++i ) {

    f_csp[i] = ZERO;

    // operate b_csp on g

    /* original summation
    for ( uint j = 0; j < NN; ++j ) {
      f_csp[i] += b_csp[j + (NN * i)] * g_csp[j];
    }
    */

    // Kahan summation
    Real c = ZERO;
    for ( uint j = 0; j < NN; ++j ) {
      Real y = ( b_csp[j + (NN * i)] * g_csp[j] ) - c;
      Real t = f_csp[i] + y;
      c = ( t - f_csp[i] ) - y;
      f_csp[i] = t;
    }

  } // end i loop


  uint M = 0; // start with no slow modes
  uint mflag = 0;
  while ( ( M < (NN - 1) ) && ( mflag == 0 ) ) {

    // testing for M = M + 1 slow modes
    // check error
    Real y_norm = ZERO;
    Real err_norm = ZERO;
    for ( uint i = 0; i < NN; ++i ) {
    //for ( uint i = NN - 1; i > 0; --i ) {

      Real sum_m = ZERO;

      /* original summation
      for ( uint k = 0; k < M + 1; ++k ) {
        sum_m += a_csp[i + (NN * k)] * f_csp[k];
      } // end k loop
      */

      // Kahan summation
      Real c = ZERO;
      for ( uint k = 0; k < M + 1; ++k ) {
        Real y = ( a_csp[i + (NN * k)] * f_csp[k] ) - c;
        Real t = sum_m + y;
        c = ( t - sum_m ) - y;
        sum_m = t;
      } // end k loop

      //////
      /*// max (infinite norm)
      if ( fabs( eps_a + ( eps_r * y[i] ) ) > y_norm ) {
        y_norm = fabs( eps_a + ( eps_r * y[i] ) );
      }

      if ( fabs ( tau[M + 1] * sum_m ) > err_norm ) {
        err_norm = fabs ( tau[M + 1] * sum_m );
      }
      */

      // if error larger than tolerance, flag
      if ( fabs ( tau[M + 1] * sum_m ) >= ( eps_a + ( eps_r * y[i] ) ) ) mflag = 1;

      //////
      /*// L2 norm
      y_norm += ( eps_a + ( eps_r * y[i] ) ) * ( eps_a + ( eps_r * y[i] ) );
      err_norm += ( tau[M + 1] * sum_m ) * ( tau[M + 1] * sum_m );
      */
      //////

      // ensure below error tolerance and not explosive mode (positive eigenvalue)
      // tau[M+1] is time scale of fastest of slow modes (driving)
      // tau[M] is time scale of slowest exhausted mode (current)

    } // end i loop

    /*// L2 norm
    y_norm = sqrt ( y_norm );
    err_norm = sqrt ( err_norm );
    */

    // add current mode to exhausted if under error tolerance and not explosive mode
    //if ( ( err_norm < eps_i ) && ( tau[M] < ZERO ) ) {
    if ( ( mflag == 0 ) && ( tau[M] < ZERO ) ) {
      M += 1; // add current mode to exhausted modes
    } else {
      mflag = 1; // explosve mode, stop here
    } // end if

  } // end while loop
  printf("M = %i\n",M);
  return M;
} // end get_fast_modes

////////////////////////////////////////////////////////////////////////
