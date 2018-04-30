/**
 * \file
 * \brief Implementation for CUDA Jacobian vector multiplication, used in exponential integrators
 *
 */

#include "sparse_multiplier.cuh"

#ifdef GENERATE_DOCS
//put this in the van der Pol namespace for documentation
namespace van_der_pol_cu {
#endif


/**
 * \brief Implements Jacobian \ vector multiplication in sparse (or unrolled) form
 * \param[in]           A           The (NSP x NSP) Jacobian matrix, see eval_jacob() for details on layout
 * \param[in]           Vm          The (NSP x 1) vector to multiply by
 * \param[out]          w           The (NSP x 1) vector to store the result in, \f$w := A * Vm\f$
 */
__device__
void sparse_multiplier(const double * A, const double * Vm, double* w) {
  w[INDEX(0)] =  A[INDEX(0)] * Vm[INDEX(0)] +  A[INDEX(4)] * Vm[INDEX(1)] + A[INDEX(8)] * Vm[INDEX(2)] + A[INDEX(12)] * Vm[INDEX(3)];
  w[INDEX(1)] =  A[INDEX(1)] * Vm[INDEX(0)] +  A[INDEX(5)] * Vm[INDEX(1)] + A[INDEX(9)] * Vm[INDEX(2)] + A[INDEX(13)] * Vm[INDEX(3)];
  w[INDEX(2)] =  A[INDEX(2)] * Vm[INDEX(0)] +  A[INDEX(6)] * Vm[INDEX(1)] + A[INDEX(10)] * Vm[INDEX(2)] + A[INDEX(14)] * Vm[INDEX(3)];
  w[INDEX(3)] =  A[INDEX(3)] * Vm[INDEX(0)] +  A[INDEX(7)] * Vm[INDEX(1)] + A[INDEX(11)] * Vm[INDEX(2)] + A[INDEX(15)] * Vm[INDEX(3)];
}


#ifdef GENERATE_DOCS
}
#endif
