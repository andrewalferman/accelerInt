#!/usr/bin/env bash

# scons ATOL=1.0e-12 RTOL=1.0e-08 t_step=0.1 t_end=3000.0 SAME_IC=True mechanism_dir=examples/van_der_pol -j 4
# touch timingdata-cvodes-vdp.csv
# touch timingdata-exp4-vdp.csv
# touch timingdata-exprb43-vdp.csv
# touch timingdata-radau2a-vdp.csv
# touch timingdata-rkc-vdp.csv
# ./cvodes-int > timingdata-cvodes-vdp.csv
# ./exp4-int > timingdata-exp4-vdp.csv
# ./exprb43-int > timingdata-exprb43-vdp.csv
# ./radau2a-int > timingdata-radau2a-vdp.csv
# ./rkc-int > timingdata-rkc-vdp.csv

# scons t_step=0.1 t_end=320.0 mechanism_dir=examples/oregonator -j 4
# touch timingdata-cvodes-oregonator.csv
# touch timingdata-exp4-oregonator.csv
# touch timingdata-exprb43-oregonator.csv
# touch timingdata-radau2a-oregonator.csv
# touch timingdata-rkc-oregonator.csv
# ./cvodes-int > timingdata-cvodes-oregonator.csv
# ./exp4-int > timingdata-exp4-oregonator.csv
# ./exprb43-int > timingdata-exprb43-oregonator.csv
# ./radau2a-int > timingdata-radau2a-oregonator.csv
# ./rkc-int > timingdata-rkc-oregonator.csv

scons ATOL=1.0e-12 RTOL=1.0e-08 t_step=1.0e-9 t_end=10.0 SAME_IC=True mechanism_dir=examples/csptest -j 4
touch timingdata-cvodes-csptest.csv
touch timingdata-exp4-csptest.csv
touch timingdata-exprb43-csptest.csv
touch timingdata-radau2a-csptest.csv
touch timingdata-rkc-csptest.csv
./cvodes-int > timingdata-cvodes-csptest.csv
./exp4-int > timingdata-exp4-csptest.csv
./exprb43-int > timingdata-exprb43-csptest.csv
./radau2a-int > timingdata-radau2a-csptest.csv
./rkc-int > timingdata-rkc-csptest.csv

# python -m pyjac -i examples/pyJac/h2.cti -l c -b out/ -ic 850.479868012,25.0000000098,H=9.06756542017123e-12,H2=0.027895867750221768,O=1.9159861828260505e-11,OH=1.249516207574034e-10,H2O=0.005443952757524492,O2=0.22127041463199268,HO2=7.783263342726486e-06,H2O2=0.00026485299052499324,AR=0.0,HE=0.0,CO=0.0,CO2=0.0,N2=0.7451171284532142
# scons ATOL=1.0e-17 RTOL=1.0e-13 t_step=1e-4 t_end=0.2 mechanism_dir=out/ -j 4
# touch timingdata-cvodes-h2.csv
# touch timingdata-exp4-h2.csv
# touch timingdata-exprb43-h2.csv
# touch timingdata-radau2a-h2.csv
# touch timingdata-rkc-h2.csv
# ./cvodes-int > timingdata-cvodes-h2.csv
# ./exp4-int > timingdata-exp4-h2.csv
# ./exprb43-int > timingdata-exprb43-h2.csv
# ./radau2a-int > timingdata-radau2a-h2.csv
# ./rkc-int > timingdata-rkc-h2.csv

# python -m pyjac -i ../Research/GRI_Mech_3/grimech30.cti -l c -b out/ -ic 900.023056449,10.0,H2=5.231897862451633e-05,H=1.5336457990843878e-12,O=1.0685174916207666e-11,O2=0.18934886996229375,OH=3.1521734970744396e-10,H2O=0.017822254359305033,HO2=2.1025611124458355e-07,H2O2=5.3315240863919146e-06,C=2.0653844285266006e-33,CH=1.9706820736239516e-24,CH2=3.82073369030059e-15,CH2\(S\)=4.1053287496827745e-16,CH3=4.877354802999909e-07,CH4=0.04652477118045588,CO=0.0019943572214920478,CO2=0.01841803199515849,HCO=8.837076579193101e-12,CH2O=0.00047606071006804347,CH2OH=9.470786316094608e-14,CH3O=2.1418149005485321e-10,CH3OH=3.490793779372462e-05,C2H=4.590600766257138e-19,C2H2=2.2084128604189096e-07,C2H3=1.678868112005144e-13,C2H4=9.188447541265724e-05,C2H5=1.0744446767736706e-10,C2H6=0.0003959647446388766,HCCO=1.568807997484106e-14,CH2CO=4.058417010270223e-06,HCCOH=1.2957433336360992e-10,N=1.3000130968519843e-16,NH=3.014188162727486e-14,NH2=3.394720372292903e-14,NH3=1.95615706711229e-09,NNH=1.1935562584491954e-16,NO=5.8801490400289775e-05,NO2=0.00018959443885181688,N2O=1.2986301697133845e-07,HNO=9.70921693531182e-12,CN=3.2164396585724215e-19,HCN=1.5571126801169642e-07,H2CN=1.9252166709560948e-14,HCNN=9.064948456139345e-25,HCNO=1.0069644049675186e-08,HOCN=1.2372464357962217e-10,HNCO=4.6108902360513814e-08,NCO=3.71220995790288e-15,C3H7=1.287424177283616e-31,C3H8=2.0097639382814184e-10,CH2CHO=1.7144189787033646e-06,CH3CHO=1.3889661870842071e-11,AR=1.6706368966120618e-07,N2=0.7245796474037304
# scons t_step=1e-4 t_end=0.4 mechanism_dir=out/ -j 4
# touch timingdata-cvodes-grimech.csv
# touch timingdata-exp4-grimech.csv
# touch timingdata-exprb43-grimech.csv
# touch timingdata-radau2a-grimech.csv
# touch timingdata-rkc-grimech.csv
# ./cvodes-int > timingdata-cvodes-grimech.csv
# ./exp4-int > timingdata-exp4-grimech.csv
# ./exprb43-int > timingdata-exprb43-grimech.csv
# ./radau2a-int > timingdata-radau2a-grimech.csv
# ./rkc-int > timingdata-rkc-grimech.csv
