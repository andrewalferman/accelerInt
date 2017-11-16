for ((n=6;n>3;n--))
{
  touch speciesdata-cvodes-1e-"$n".csv
  touch speciesdata-exp4-1e-"$n".csv
  touch speciesdata-exprb43-1e-"$n".csv
  touch speciesdata-radau2a-1e-"$n".csv
  touch speciesdata-rkc-1e-"$n".csv
  scons t_step=1e-"$n" t_end=1e-"$n" -j 4
  ./cvodes-int 4 450900 > speciesdata-cvodes-1e-"$n".csv
  ./cvodes-int 4 450900 > speciesdata-exp4-1e-"$n".csv
  ./cvodes-int 4 450900 > speciesdata-exprb43-1e-"$n".csv
  ./cvodes-int 4 450900 > speciesdata-radau2a-1e-"$n".csv
  ./cvodes-int 4 450900 > speciesdata-rkc-1e-"$n".csv
  mv speciesdata-cvodes-1e-"$n".csv ../Research/accelerInt_Data/
  mv speciesdata-exp4-1e-"$n".csv ../Research/accelerInt_Data/
  mv speciesdata-exprb43-1e-"$n".csv ../Research/accelerInt_Data/
  mv speciesdata-radau2a-1e-"$n".csv ../Research/accelerInt_Data/
  mv speciesdata-rkc-1e-"$n".csv ../Research/accelerInt_Data/
}
