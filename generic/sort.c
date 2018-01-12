/** File containing insertion sort function.
 * \file sort.c
 *
 * \author Kyle E. Niemeyer
 * \date 08/02/2011
 */

/** Include common libraries and variables. */
#include "head.h"

/** array size threshold for quicksort, from Numerical Recipes. */
#define M 7
/** maximum stack size for quicksort, from Numerical Recipes. */
#define NSTACK 50

/** Insertion sort function.
 *
 * Performs insertion sort and returns indices in ascending order based on 
 * input values. Best on very small arrays (n < 20). Based on Numerical Recipes
 * in C algorithm.
 * 
 * \param[in]   n  	  size of array
 * \param[in]   vals  array of values to be sorted
 * \param[out]  ord   array of sorted indices
 */
void insertion_sort ( uint n, Real * vals, uint * ord ) {
  
  for ( int i = 0; i < n; ++i ) {
    ord[i] = i;
  }
  
  for ( int j = 1; j < n; ++j ) {
    
    // pick out one element
    Real val = vals[ ord[j] ];
    uint ival = ord[j];
    
    int i = j - 1;
    while ( ( i >= 0 ) && ( vals[ ord[i] ] > val ) ) { // look where to insert
      ord[i + 1] = ord[i];
      --i;
    }
    
    // insert element
    ord[i + 1] = ival;
  }
  
}