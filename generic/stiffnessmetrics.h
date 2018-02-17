/**
 * \file
 * \brief header file for stiffness metric calculator function
 *
 * \author Andrew Alferman
 * \date 12/14/2017
 *
 *
 */

#include "header.h"
#include "jacob.h"

 void calculatemetrics(double* y_local, double pr_local, double* stiffratio,
                      double* stiffindicator, double* CEM, double* CSP,
                      const double t, const double t_end);
