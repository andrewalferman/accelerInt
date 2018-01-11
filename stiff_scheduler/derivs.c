/** Derivative file for CSP model problem project.
 * \file derivs.c
 *
 * \author Kyle E. Niemeyer
 * \date 08/02/2011
 *
 * Contains derivative and Jacobian for stiff model problem.
 */

/** Include common code. */
#include "head.h"

////////////////////////////////////////////////////////////////////////

/** Derivative (dydt) source term.
 *
 * \param[in]  tim  	the time (sec)
 * \param[in]  y	 	  the data array, size neq
 * \param[in]  Q      slow-manifold projector matrix, size NN*NN
 * \param[in]  qflag  projector flag (1 to use)
 * \param[out] ydot	  derivative array, size neq
 */
void dydt ( Real tim, Real * y, Real * Q, int qflag, Real * ydot ) {
  
  Real eps1 = ONE / eps;
  Real eps2 = ONE / ( eps * eps );
  
  Real y2 = ONE + y[1];
  
  Real y3 = ONE + y[2];
  Real y32 = y3 * y3;
  
  Real y4 = ONE + y[3];
  Real y42 = y4 * y4;
  
  /*
	ydot[0] = ( eps1 * eps2 ) * ( -y[0] + ( y[1] / y2 ) )
          - ( y[1] / ( y2 * y2 ) ) - ( y[2] / y32 )
          - ( y[3] / y42 );
          
	ydot[1] = eps2 * ( -y[1] + ( y[2] / y3 ) )
          - ( y[2] / y32 ) - ( y[3] / y42 );
          
	ydot[2] = eps1 * ( -y[2] + ( y[3] / y4 ) )
          - ( y[3] / y42 );
  */
  ydot[0] = ( ( ( ONE / (eps * eps * eps) ) * ( -y[0] + ( y[1] / (ONE + y[1]) ) )
          - ( y[1] / ( (y[1] + ONE) * (y[1] + ONE) ) ) )
          - ( y[2] / ( (y[2] + ONE) * (y[2] + ONE) ) ) )
          - ( y[3] / ( (y[3] + ONE) * (y[3] + ONE) ) );
  
  ydot[1] = ( ( ONE / (eps * eps) ) * ( -y[1] + ( y[2] / (ONE + y[2]) ) )
          - ( y[2] / ( (y[2] + ONE) * (y[2] + ONE) ) ) )
          - ( y[3] / ( (y[3] + ONE) * (y[3] + ONE) ) );
          
  ydot[2] = ( ONE / eps ) * ( -y[2] + ( y[3] / (ONE + y[3]) ) )
          - ( y[3] / ( (y[3] + ONE) * (y[3] + ONE) ) );
  
	ydot[3] = -y[3];
  
#if !defined(EXPLICIT)
  // slow-manifold projection
  if ( qflag == 1 ) {
    
    Real ydotn[NN];
    
    for ( uint i = 0; i < NN; ++i ) {
      
      Real sum_row = ZERO;
      for ( uint j = 0; j < NN; ++j ) {
        sum_row += Q[i + (NN * j)] * ydot[j];
      } // end j loop
      
      ydotn[i] = sum_row;
    } // end i loop
    
    // now replace ydot with ydotn
    for ( uint i = 0; i < NN; ++i ) {
      ydot[i] = ydotn[i];
    } // end i loop
    
  } // end if
#endif
  
}

////////////////////////////////////////////////////////////////////////
  
/** Jacobian for model problem.
 *
 * \param[in]  tim  	the time (sec)
 * \param[in]  y	 	the data array, size neq
 * \param[out] dfdy	Jacobian, size (neq,neq)
 */
void jacob ( Real tim, Real * y, Real * dfdy ) {
  
  Real eps1 = ONE / eps;
  Real eps2 = ONE / ( eps * eps );
  
  Real y22 = ONE / ( ( ONE + y[1] ) * ( ONE + y[1] ) );
  
  Real y32 = ONE / ( ( ONE + y[2] ) * ( ONE + y[2] ) );
  Real y33 = y32 / ( ONE + y[2] );
  
  Real y42 = ONE / ( ( ONE + y[3] ) * ( ONE + y[3] ) );
  Real y43 = y42 / ( ONE + y[3] );
  
  // dg0/dyi
	/*
  dfdy[0]  = -eps1 * eps2;
	dfdy[4]  = ( eps1 * eps2 * y22 ) + ( ( y[1] - ONE ) * y22 / ( ONE + y[1] ) );
	dfdy[8]  = ( y[2] - ONE ) * y33;
	dfdy[12] = ( y[3] - ONE ) * y43;
	*/
  dfdy[0]  = -ONE / ( eps * eps * eps );
  dfdy[4]  = ONE / ( eps * eps * eps * ( ONE + y[1] ) * ( ONE + y[1] ) )
           + ( y[1] - ONE ) / ( (y[1] + ONE) * (y[1] + ONE) * (y[1] + ONE) );
  dfdy[8]  = ( y[2] - ONE ) / ( (y[2] + ONE) * (y[2] + ONE) * (y[2] + ONE) );
  dfdy[12] = ( y[3] - ONE ) / ( (y[3] + ONE) * (y[3] + ONE) * (y[3] + ONE) );
  
  // dg1/dyi
	/*
  dfdy[1]  = ZERO;
	dfdy[5]  = -eps2;
	dfdy[9]  = ( eps2 * y32 ) + ( ( y[2] - ONE ) * y33 );
	dfdy[13] = ( y[3] - ONE ) * y43;
	*/
  dfdy[1]  = ZERO;
  dfdy[5]  = -ONE / ( eps * eps );
  dfdy[9]  = ONE / ( eps * eps * ( ONE + y[2] ) * ( ONE + y[2] ) )
           + ( y[2] - ONE ) / ( (y[2] + ONE) * (y[2] + ONE) * (y[2] + ONE) );
  dfdy[13] = ( y[3] - ONE ) / ( (y[3] + ONE) * (y[3] + ONE) * (y[3] + ONE) );
  
  // dg2/dyi
	/*
  dfdy[2]  = ZERO;
	dfdy[6]  = ZERO;
	dfdy[10] = -eps1;
	dfdy[14] = ( eps1 * y42 ) + ( ( y[3] - ONE ) * y43 );
	*/
  dfdy[2]  = ZERO;
  dfdy[6]  = ZERO;
  dfdy[10] = -ONE / eps;
  dfdy[14] = ONE / ( eps * ( ONE + y[3] ) * ( ONE + y[3] ) )
           + ( y[3] - ONE ) / ( (y[3] + ONE) * (y[3] + ONE) * (y[3] + ONE) );
  
  // dg3/dyi
	dfdy[3]  = ZERO;
	dfdy[7]  = ZERO;
	dfdy[11] = ZERO;
	dfdy[15] = -ONE;
	
}