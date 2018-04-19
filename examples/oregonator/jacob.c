/**
 * \file
 * \brief An implementation of the Oregonator jacobian \f$\frac{\partial \dot{\vec{y}}}{\partial \vec{y}}\f$
 *
 * Implementes the evaluation of Oregonator Jacobian
 *
 */

#include "header.h"

#ifdef GENERATE_DOCS
//put this in the Oregonator namespace for documentation
namespace oregonator {
#endif

/**
 * \brief An implementation of the Oregonator jacobian
 *
 * \param[in]           t               The current system time
 * \param[in]           mu              Dummy parameter, needed for compatibility only
 * \param[in]           y               The state vector at time t
 * \param[out]          jac             The jacobian to populate
 *
 *  The Jacobian is in a local Column-major (Fortran) order.  As with dydt(), this function operates on local
 *  copies of the global state vector and jacobian.  Hence simple linear indexing can be used here.
 *  @see solver_generic.c
 *
 */
void eval_jacob (const double t, const double mu, const double * __restrict__ y, double * __restrict__ jac)
{

  double s = 77.27;
  double q = 8.375E-6;
  double w = 0.161;

  jac[0] = s * (-1 * y[1] + 1 - q * 2 * y[0]);
  jac[1] = s * (1 - y[0]);
  jac[2] = 0;
  jac[3] = -1 * y[1] / s;
  jac[4] = (-1 - y[0]) / s;
  jac[5] = 1 / s;
  jac[6] = w;
  jac[7] = 0;
  jac[8] = -1 * w;

    // //Note, to reach index [i, j] of the Jacobian, we multiply `i` by NSP, the size of the first dimension of Jacobian and add j, i.e.:
    // //jac[i, j] -> jac[i * NSP + j]
    // //!jac[0, 0] = \f$\frac{\partial \dot{y_1}}{\partial y_1}\f$
    // jac[0 * NSP + 0] = 0;
    // //!jac[0, 1] = \f$\frac{\partial \dot{y_2}}{\partial y_1}\f$
    // jac[0 * NSP + 1] = -2 * mu * y[0] * y[1] - 1;
    // //!jac[1, 0] = \f$\frac{\partial \dot{y_1}}{\partial y_2}\f$
    // jac[1 * NSP + 0] = 1;
    // //!jac[1, 1] = \f$\frac{\partial \dot{y_2}}{\partial y_2}\f$
    // jac[1 * NSP + 1] = mu * (1 - y[0] * y[0]);
}

#ifdef GENERATE_DOCS
}
#endif
