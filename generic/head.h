/**
 * @mainpage CSP-based integration of stiff model problem.
 *
 * @author <a href="mailto:ken7@case.edu">Kyle E. Niemeyer</a> and Jerry C. Lee
 *
 * This program performs the explicit integration of a stiff problem, using the
 * computational singular perturbation (CSP) method to separate the fast and slow
 * modes and generate a slow manifold projection matrix to eliminate the stiffness
 * (and therefore allow explicit integration).
 *
 * Change stiffness using stiffness factor eps, and the CSP error tolerance
 * parameters eps_r and eps_i in head.h.
 */

/** Header file for CSP model problem project.
 * \file head.h
 *
 * \author Kyle E. Niemeyer
 * \date 08/02/2011
 *
 * Contains libraries, definitions, and constants.
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>
#include "header.h"
#include "solver_options.h"

// /** Problem size definition. */
// #define NN 4

// /** Number of "threads" (simultaneous integrations) */
// #define NUM 1

/** Stiffness factor eps. */
extern const double eps;

/** Flag explicit only integration (no CSP projection) */
//#define EXPLICIT

/** Error/warning message flag.*/
#define ERROR

/** Sets precision as double or float. */
#define DOUBLE
#ifdef DOUBLE
  /** Define Real as double. */
  #define Real double

  /** Double precision ZERO. */
  #define ZERO 0.0
  /** Double precision ONE. */
  #define ONE 1.0
  /** Double precision TWO. */
  #define TWO 2.0
  /** Double precision THREE. */
  #define THREE 3.0
  /** Double precision FOUR. */
  #define FOUR 4.0

  // /** Machine precision constant. */
  // #define SMALL DBL_EPSILON
#else
  /** Define Real as float. */
  #define Real float

  /** Single precision ZERO. */
  #define ZERO 0.0f
  /** Single precision ONE. */
  #define ONE 1.0f
  /** Single precision (float) TWO. */
  #define TWO 2.0f
  /** Single precision THREE. */
  #define THREE 3.0f
  /** Single precision FOUR. */
  #define FOUR 4.0f

  /** Machine precision constant. */
  #define SMALL FLT_EPSILON
#endif

/** Unsigned int typedef. */
typedef unsigned int uint;
/** Unsigned short int typedef. */
typedef unsigned short int usint;

uint get_slow_projector ( Real tim, Real * y, Real * Qs, Real * taum1, Real * Rc );
void radical_correction ( Real tim, Real * y, Real * Rc, Real * g );
void RK4 ( Real t, Real h, Real * y0, Real * Q, Real * y );
void RK2 ( Real t, Real h, Real * y0, Real * Q, Real * y );
void euler ( Real t, Real h, Real * y0, Real * Q, Real * y );
void RKB6 ( Real t, Real h, Real * y0, Real * Q, Real * y );
