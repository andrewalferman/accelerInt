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
        y_specs[8] = temp[8];
        y_specs[9] = temp[9];
        y_specs[10] = temp[10];
        y_specs[11] = temp[11];
        y_specs[12] = temp[12];
        y_specs[13] = temp[13];
        y_specs[14] = temp[14];
        y_specs[15] = temp[15];
        y_specs[16] = temp[16];
        y_specs[17] = temp[17];
        y_specs[18] = temp[18];
        y_specs[19] = temp[19];
        y_specs[20] = temp[20];
        y_specs[21] = temp[21];
        y_specs[22] = temp[22];
        y_specs[23] = temp[23];
        y_specs[24] = temp[24];
        y_specs[25] = temp[25];
        y_specs[26] = temp[26];
        y_specs[27] = temp[27];
        y_specs[28] = temp[28];
        y_specs[29] = temp[29];
        y_specs[30] = temp[30];
        y_specs[31] = temp[31];
        y_specs[32] = temp[32];
        y_specs[33] = temp[33];
        y_specs[34] = temp[34];
        y_specs[35] = temp[35];
        y_specs[36] = temp[36];
        y_specs[37] = temp[37];
        y_specs[38] = temp[38];
        y_specs[39] = temp[39];
        y_specs[40] = temp[40];
        y_specs[41] = temp[41];
        y_specs[42] = temp[42];
        y_specs[43] = temp[43];
        y_specs[44] = temp[44];
        y_specs[45] = temp[45];
        y_specs[46] = temp[46];
        y_specs[47] = temp[48];
        y_specs[48] = temp[49];
        y_specs[49] = temp[50];
        y_specs[50] = temp[51];
        y_specs[51] = temp[52];
        y_specs[52] = temp[47];
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
        y_specs[8] = temp[8];
        y_specs[9] = temp[9];
        y_specs[10] = temp[10];
        y_specs[11] = temp[11];
        y_specs[12] = temp[12];
        y_specs[13] = temp[13];
        y_specs[14] = temp[14];
        y_specs[15] = temp[15];
        y_specs[16] = temp[16];
        y_specs[17] = temp[17];
        y_specs[18] = temp[18];
        y_specs[19] = temp[19];
        y_specs[20] = temp[20];
        y_specs[21] = temp[21];
        y_specs[22] = temp[22];
        y_specs[23] = temp[23];
        y_specs[24] = temp[24];
        y_specs[25] = temp[25];
        y_specs[26] = temp[26];
        y_specs[27] = temp[27];
        y_specs[28] = temp[28];
        y_specs[29] = temp[29];
        y_specs[30] = temp[30];
        y_specs[31] = temp[31];
        y_specs[32] = temp[32];
        y_specs[33] = temp[33];
        y_specs[34] = temp[34];
        y_specs[35] = temp[35];
        y_specs[36] = temp[36];
        y_specs[37] = temp[37];
        y_specs[38] = temp[38];
        y_specs[39] = temp[39];
        y_specs[40] = temp[40];
        y_specs[41] = temp[41];
        y_specs[42] = temp[42];
        y_specs[43] = temp[43];
        y_specs[44] = temp[44];
        y_specs[45] = temp[45];
        y_specs[46] = temp[46];
        y_specs[47] = temp[52];
        y_specs[48] = temp[47];
        y_specs[49] = temp[48];
        y_specs[50] = temp[49];
        y_specs[51] = temp[50];
        y_specs[52] = temp[51];
    }
void set_same_initial_conditions(int NUM, double** y_host, double** var_host) 
{
    double Xi [NSP] = {0.0};
    //set initial mole fractions here

    //Normalize mole fractions to sum to one
    double Xsum = 0.0;
    Xi[0] = 5.231897862451633e-05;
    Xi[1] = 1.5336457990843878e-12;
    Xi[2] = 1.0685174916207666e-11;
    Xi[3] = 0.18934886996229375;
    Xi[4] = 3.1521734970744396e-10;
    Xi[5] = 0.017822254359305033;
    Xi[6] = 2.1025611124458355e-07;
    Xi[7] = 5.3315240863919146e-06;
    Xi[8] = 2.0653844285266006e-33;
    Xi[9] = 1.9706820736239516e-24;
    Xi[10] = 3.82073369030059e-15;
    Xi[11] = 4.1053287496827745e-16;
    Xi[12] = 4.877354802999909e-07;
    Xi[13] = 0.04652477118045588;
    Xi[14] = 0.0019943572214920478;
    Xi[15] = 0.01841803199515849;
    Xi[16] = 8.837076579193101e-12;
    Xi[17] = 0.00047606071006804347;
    Xi[18] = 9.470786316094608e-14;
    Xi[19] = 2.1418149005485321e-10;
    Xi[20] = 3.490793779372462e-05;
    Xi[21] = 4.590600766257138e-19;
    Xi[22] = 2.2084128604189096e-07;
    Xi[23] = 1.678868112005144e-13;
    Xi[24] = 9.188447541265724e-05;
    Xi[25] = 1.0744446767736706e-10;
    Xi[26] = 0.0003959647446388766;
    Xi[27] = 1.568807997484106e-14;
    Xi[28] = 4.058417010270223e-06;
    Xi[29] = 1.2957433336360992e-10;
    Xi[30] = 1.3000130968519843e-16;
    Xi[31] = 3.014188162727486e-14;
    Xi[32] = 3.394720372292903e-14;
    Xi[33] = 1.95615706711229e-09;
    Xi[34] = 1.1935562584491954e-16;
    Xi[35] = 5.8801490400289775e-05;
    Xi[36] = 0.00018959443885181688;
    Xi[37] = 1.2986301697133845e-07;
    Xi[38] = 9.70921693531182e-12;
    Xi[39] = 3.2164396585724215e-19;
    Xi[40] = 1.5571126801169642e-07;
    Xi[41] = 1.9252166709560948e-14;
    Xi[42] = 9.064948456139345e-25;
    Xi[43] = 1.0069644049675186e-08;
    Xi[44] = 1.2372464357962217e-10;
    Xi[45] = 4.6108902360513814e-08;
    Xi[46] = 3.71220995790288e-15;
    Xi[48] = 1.287424177283616e-31;
    Xi[49] = 2.0097639382814184e-10;
    Xi[50] = 1.7144189787033646e-06;
    Xi[51] = 1.3889661870842071e-11;
    Xi[47] = 1.6706368966120618e-07;
    Xi[52] = 0.7245796474037304;
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
    double P = 1013250.0;
    // set intial temperature, units [K]
    double T0 = 900.023056449;

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

