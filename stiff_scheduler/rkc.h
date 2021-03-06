#ifndef RKC_HEAD
#define RKC_HEAD

#include "header.h"
#include "rkc_props.h"

//Real rkc_spec_rad (const Real, const Real, const Real*, const Real*, Real*, Real*);
Real rkc_spec_rad (const Real, const Real, const Real, const Real*, const Real*, Real*, Real*);
void rkc_step (const Real, const Real, const Real, const Real*, const Real*, const int, Real*);
int integrateRKC (Real t, const Real tEnd, const Real pr, Real* y);

#endif
