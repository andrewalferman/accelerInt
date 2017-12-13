#include "mass_mole.h"
#include <stdio.h>
#include "mechanism.h"
    //apply masking of ICs for cache optimized mechanisms
    void apply_mask(double* y_specs) {
        double temp [NSP];
        memcpy(temp, y_specs, NSP * sizeof(double));
        y_specs[0] = temp[0];
        y_specs[1] = temp[1];
        y_specs[2] = temp[2];
        y_specs[3] = temp[3];
        y_specs[4] = temp[4];
        y_specs[5] = temp[5];
        y_specs[6] = temp[6];
        y_specs[7] = temp[7];
        y_specs[8] = temp[9];
        y_specs[9] = temp[10];
        y_specs[10] = temp[11];
        y_specs[11] = temp[12];
        y_specs[12] = temp[8];
    }
    //reverse masking of ICs for cache optimized mechanisms
    void apply_reverse_mask(double* y_specs) {
        double temp [NSP];
        memcpy(temp, y_specs, NSP * sizeof(double));
        y_specs[0] = temp[0];
        y_specs[1] = temp[1];
        y_specs[2] = temp[2];
        y_specs[3] = temp[3];
        y_specs[4] = temp[4];
        y_specs[5] = temp[5];
        y_specs[6] = temp[6];
        y_specs[7] = temp[7];
        y_specs[8] = temp[12];
        y_specs[9] = temp[8];
        y_specs[10] = temp[9];
        y_specs[11] = temp[10];
        y_specs[12] = temp[11];
    }
void set_same_initial_conditions(int NUM, double** y_host, double** var_host) 
{
    double Xi [NSP] = {0.0};
    //set initial mole fractions here

    //Normalize mole fractions to sum to one
    double Xsum = 0.0;
    Xi[0] = 1.1636247155941894e-05;
    Xi[1] = 0.022754941415572758;
    Xi[2] = 1.9895064116812003e-05;
    Xi[3] = 2.0993869462435053e-05;
    Xi[4] = 0.05131099405277086;
    Xi[5] = 0.18045534227727603;
    Xi[6] = 0.0002361082973661372;
    Xi[7] = 7.296032306490985e-05;
    Xi[12] = 0.7451171284532138;
    Xi[8] = 0.0;
    Xi[9] = 0.0;
    Xi[10] = 0.0;
    Xi[11] = 0.0;
    for (int j = 0; j < NSP; ++ j) {
        Xsum += Xi[j];
    }
    if (Xsum == 0.0) {
        printf("Use of the set initial conditions function requires user implementation!\n");
        exit(-1);
    }
    for (int j = 0; j < NSP; ++ j) {
        Xi[j] /= Xsum;
    }

    //convert to mass fractions
    double Yi[NSP - 1] = {0.0};
    mole2mass(Xi, Yi);

    //set initial pressure, units [PA]
    double P = 101325.0;
    // set intial temperature, units [K]
    double T0 = 881.9061350328426;

    (*y_host) = (double*)malloc(NUM * NSP * sizeof(double));
    (*var_host) = (double*)malloc(NUM * sizeof(double));
    //load temperature and mass fractions for all threads (cells)
    for (int i = 0; i < NUM; ++i) {
        (*y_host)[i] = T0;
        //loop through species
        for (int j = 1; j < NSP; ++j) {
            (*y_host)[i + NUM * j] = Yi[j - 1];
        }
    }

#ifdef CONV
    //calculate density
    double rho = getDensity(T0, P, Xi);
#endif

    for (int i = 0; i < NUM; ++i) {
#ifdef CONV
        (*var_host)[i] = rho;
#elif defined(CONP)
        (*var_host)[i] = P;
#endif
    }
}

