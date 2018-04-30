scons t_step=0.1 t_end=3000.0 mechanism_dir=examples/van_der_pol -j 4
touch timingdata-cvodes-vdp.csv
touch timingdata-exp4-vdp.csv
touch timingdata-exprb43-vdp.csv
touch timingdata-radau2a-vdp.csv
touch timingdata-rkc-vdp.csv
./cvodes-int > timingdata-cvodes-vdp.csv
./exp4-int > timingdata-exp4-vdp.csv
./exprb43-int > timingdata-exprb43-vdp.csv
./radau2a-int > timingdata-radau2a-vdp.csv
./rkc-int > timingdata-rkc-vdp.csv

scons t_step=0.1 t_end=320.0 mechanism_dir=examples/oregonator -j 4
touch timingdata-cvodes-oregonator.csv
touch timingdata-exp4-oregonator.csv
touch timingdata-exprb43-oregonator.csv
touch timingdata-radau2a-oregonator.csv
touch timingdata-rkc-oregonator.csv
./cvodes-int > timingdata-cvodes-oregonator.csv
./exp4-int > timingdata-exp4-oregonator.csv
./exprb43-int > timingdata-exprb43-oregonator.csv
./radau2a-int > timingdata-radau2a-oregonator.csv
./rkc-int > timingdata-rkc-oregonator.csv

python -m pyjac -i examples/pyJac/h2.cti -l c -b out/ -ic 850.479868012,25.0000000098,H=9.06756542017123e-12,H2=0.027895867750221768,O=1.9159861828260505e-11,OH=1.249516207574034e-10,H2O=0.005443952757524492,O2=0.22127041463199268,HO2=7.783263342726486e-06,H2O2=0.00026485299052499324,AR=0.0,HE=0.0,CO=0.0,CO2=0.0,N2=0.7451171284532142
python -m pyjac -i ../Research/H2_CO/out/mechanism.cti -l c
scons t_step=0.1 t_end=320.0 mechanism_dir=examples/oregonator -j 4
touch timingdata-cvodes-oregonator.csv
touch timingdata-exp4-oregonator.csv
touch timingdata-exprb43-oregonator.csv
touch timingdata-radau2a-oregonator.csv
touch timingdata-rkc-oregonator.csv
./cvodes-int > timingdata-cvodes-oregonator.csv
./exp4-int > timingdata-exp4-oregonator.csv
./exprb43-int > timingdata-exprb43-oregonator.csv
./radau2a-int > timingdata-radau2a-oregonator.csv
./rkc-int > timingdata-rkc-oregonator.csv
