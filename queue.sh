#!/usr/bin/env bash

python -m pyjac -i examples/pyJac/h2.cti -l c -b out/
file=./ign_data.bin
if [ -e "$file" ]; then
  echo "Deleting old data.bin file"
  rm -f "$file"
  cp ./initials/H2_CO/ign_data.bin ./ign_data.bin
else
  cp ./initials/H2_CO/ign_data.bin ./ign_data.bin
fi
for ((n=8;n>3;n--))
{
  touch speciesdata-h2-cvodes-1e-"$n"-ht.csv
  touch speciesdata-h2-exp4-1e-"$n"-ht.csv
  touch speciesdata-h2-exprb43-1e-"$n"-ht.csv
  touch speciesdata-h2-radau2a-1e-"$n"-ht.csv
  touch speciesdata-h2-rkc-1e-"$n"-ht.csv
  scons ATOL=1.0e-10 RTOL=1.0e-6 t_step=1e-"$n" t_end=1e-"$n" SAME_IC=False -j 4
  ./cvodes-int 4 900900 > speciesdata-h2-cvodes-1e-"$n"-ht.csv
  ./exp4-int 4 900900 > speciesdata-h2-exp4-1e-"$n"-ht.csv
  ./exprb43-int 4 900900 > speciesdata-h2-exprb43-1e-"$n"-ht.csv
  ./radau2a-int 4 900900 > speciesdata-h2-radau2a-1e-"$n"-ht.csv
  ./rkc-int 4 900900 > speciesdata-h2-rkc-1e-"$n"-ht.csv
  #mv speciesdata-cvodes-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-exp4-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-exprb43-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-radau2a-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-rkc-1e-"$n".csv ../Research/accelerInt_Data/
}

python -m pyjac -i ../Research/GRI_Mech_3/grimech30.cti -l c -b out/
python data_bin_writer.py -d ../Research/GRI_Mech_3/
rm -f ./ign_data.bin
cp ./initials/GRI_Mech_3/ign_data.bin ./ign_data.bin
for ((n=8;n>3;n--))
{
  touch speciesdata-grimech-cvodes-1e-"$n"-ht.csv
  touch speciesdata-grimech-exp4-1e-"$n"-ht.csv
  touch speciesdata-grimech-exprb43-1e-"$n"-ht.csv
  touch speciesdata-grimech-radau2a-1e-"$n"-ht.csv
  touch speciesdata-grimech-rkc-1e-"$n"-ht.csv
  scons ATOL=1.0e-10 RTOL=1.0e-6 t_step=1e-"$n" t_end=1e-"$n" -j 4
  ./cvodes-int 4 450900 > speciesdata-grimech-cvodes-1e-"$n"-ht.csv
  ./exp4-int 4 450900 > speciesdata-grimech-exp4-1e-"$n"-ht.csv
  ./exprb43-int 4 450900 > speciesdata-grimech-exprb43-1e-"$n"-ht.csv
  ./radau2a-int 4 450900 > speciesdata-grimech-radau2a-1e-"$n"-ht.csv
  ./rkc-int 4 450900 > speciesdata-grimech-rkc-1e-"$n"-ht.csv
  #mv speciesdata-cvodes-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-exp4-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-exprb43-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-radau2a-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-rkc-1e-"$n".csv ../Research/accelerInt_Data/
}
