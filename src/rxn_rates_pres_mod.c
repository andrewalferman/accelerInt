#include <stdlib.h>
#include <math.h>
#include "head.h"

void get_rxn_pres_mod (const Real T, const Real pres, const Real * C, Real * pres_mod) {
  // third body variable declaration
  Real thd;

  // pressure dependence variable declarations
  Real k0;
  Real kinf;
  Real Pr;

  // troe variable declarations
  Real logPr;
  Real logFcent;
  Real A;
  Real B;

  Real logT = log(T);
  Real m = pres / (8.31451000e+07 * T);

  // reaction 5;
  pres_mod[0] = m + 1.5 * C[1] + 11.0 * C[4] + 0.9 * C[11] + 2.8 * C[12] - 1.0 * C[9] - 1.0 * C[10];
  // reaction 8;
  pres_mod[1] = m + 1.5 * C[1] + 11.0 * C[4] - 1.0 * C[9] - 1.0 * C[10] + 0.9 * C[11] + 2.8 * C[12];
  // reaction 11;
  pres_mod[2] = m + 1.5 * C[1] + 11.0 * C[4] - 0.25 * C[9] - 0.25 * C[10] + 0.9 * C[11] + 2.8 * C[12];
  // reaction 12;
  pres_mod[3] = m + 2.0 * C[1] - 1.0 * C[4] + 0.1 * C[10] + 1.0 * C[8] + 0.5 * C[5] + 0.9 * C[11] + 2.8 * C[12];
  // reaction 14;
  thd = m + 1.0 * C[1] + 13.0 * C[4] - 0.22 * C[5] + 0.9 * C[11] + 2.8 * C[12] - 0.33 * C[9] - 0.2 * C[10];
  k0 = exp(4.79026732e+01 - 1.72 * logT - (2.64088106e+02 / T));
  kinf = exp(2.91680690e+01 + 0.44 * logT);
  Pr = k0 * thd / kinf;
  logPr = log10(Pr);
  logFcent = log10( 5.00000000e-01 * exp(-T / 1.00000000e-30) + 5.00000000e-01 * exp(-T / 1.00000000e+30) );
  A = logPr - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * logPr;
  pres_mod[4] = exp(log(10.0) * logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);
  // reaction 21;
  thd = m + 6.5 * C[4] + 0.6 * C[12] + 0.5 * C[8] + 0.2 * C[5] - 0.35 * C[10] + 6.7 * C[7] + 2.7 * C[1] + 1.8 * C[11];
  k0 = exp(5.61743249e+01 - 2.3 * logT - (2.45313092e+04 / T));
  kinf = exp(2.83241683e+01 + 0.9 * logT - (2.45313092e+04 / T));
  Pr = k0 * thd / kinf;
  logPr = log10(Pr);
  logFcent = log10( 5.70000000e-01 * exp(-T / 1.00000000e-30) + 4.30000000e-01 * exp(-T / 1.00000000e+30) );
  A = logPr - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * logPr;
  pres_mod[5] = exp(log(10.0) * logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);
} // end get_rxn_pres_mod

