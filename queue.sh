# for ((n=8;n>3;n--))
# {
#   touch speciesdata-cvodes-1e-"$n".csv
#   touch speciesdata-exp4-1e-"$n".csv
#   touch speciesdata-exprb43-1e-"$n".csv
#   touch speciesdata-radau2a-1e-"$n".csv
#   touch speciesdata-rkc-1e-"$n".csv
#   scons t_step=1e-"$n" t_end=1e-"$n" -j 4
#   ./cvodes-int 4 450900 > speciesdata-cvodes-1e-"$n".csv
#   ./exp4-int 4 450900 > speciesdata-exp4-1e-"$n".csv
#   ./exprb43-int 4 450900 > speciesdata-exprb43-1e-"$n".csv
#   ./radau2a-int 4 450900 > speciesdata-radau2a-1e-"$n".csv
#   ./rkc-int 4 450900 > speciesdata-rkc-1e-"$n".csv
#   mv speciesdata-cvodes-1e-"$n".csv ../Research/accelerInt_Data/
#   mv speciesdata-exp4-1e-"$n".csv ../Research/accelerInt_Data/
#   mv speciesdata-exprb43-1e-"$n".csv ../Research/accelerInt_Data/
#   mv speciesdata-radau2a-1e-"$n".csv ../Research/accelerInt_Data/
#   mv speciesdata-rkc-1e-"$n".csv ../Research/accelerInt_Data/
# }
touch timingdata-cvodes-1e-4.csv
touch timingdata-exp4-1e-4.csv
touch timingdata-exprb43-1e-4.csv
touch timingdata-radau2a-1e-4.csv
touch timingdata-rkc-1e-4.csv
./cvodes-int 1 1 > timingdata-cvodes-1e-4.csv
./exp4-int 1 1 > timingdata-exp4-1e-4.csv
./exprb43-int 1 1 > timingdata-exprb43-1e-4.csv
./radau2a-int 1 1 > timingdata-radau2a-1e-4.csv
./rkc-int 1 1 > timingdata-rkc-1e-4.csv
