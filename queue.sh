python -m pyjac -i examples/pyJac/h2.cti -l c -b out/
python data_bin_writer.py -d ../Research/H2_CO/
for ((n=8;n>3;n--))
{
  #touch speciesdata-h2-cvodes-1e-"$n".csv
  #touch speciesdata-h2-exp4-1e-"$n".csv
  #touch speciesdata-h2-exprb43-1e-"$n".csv
  #touch speciesdata-h2-radau2a-1e-"$n".csv
  #touch speciesdata-h2-rkc-1e-"$n".csv
  scons ATOL=1.0e-17 RTOL=1.0e-13 t_step=1e-"$n" t_end=1e-"$n" SAME_IC=False -j 4
  ./cvodes-int 4 450900 #> speciesdata-h2-cvodes-1e-"$n".csv
  ./exp4-int 4 450900 #> speciesdata-h2-exp4-1e-"$n".csv
  ./exprb43-int 4 450900 #> speciesdata-h2-exprb43-1e-"$n".csv
  ./radau2a-int 4 450900 #> speciesdata-h2-radau2a-1e-"$n".csv
  ./rkc-int 4 450900 #> speciesdata-h2-rkc-1e-"$n".csv
  #mv speciesdata-cvodes-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-exp4-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-exprb43-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-radau2a-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-rkc-1e-"$n".csv ../Research/accelerInt_Data/
}

python -m pyjac -i ../Research/GRI_Mech_3/grimech30.cti -l c -b out/
python data_bin_writer.py -d ../Research/GRI_Mech_3/
for ((n=8;n>3;n--))
{
  #touch speciesdata-grimech-cvodes-1e-"$n".csv
  #touch speciesdata-grimech-exp4-1e-"$n".csv
  #touch speciesdata-grimech-exprb43-1e-"$n".csv
  #touch speciesdata-grimech-radau2a-1e-"$n".csv
  #touch speciesdata-grimech-rkc-1e-"$n".csv
  scons ATOL=1.0e-17 RTOL=1.0e-13 t_step=1e-"$n" t_end=1e-"$n" -j 4
  ./cvodes-int 4 450900 #> speciesdata-grimech-cvodes-1e-"$n".csv
  ./exp4-int 4 450900 #> speciesdata-grimech-exp4-1e-"$n".csv
  ./exprb43-int 4 450900 #> speciesdata-grimech-exprb43-1e-"$n".csv
  ./radau2a-int 4 450900 #> speciesdata-grimech-radau2a-1e-"$n".csv
  ./rkc-int 4 450900 #> speciesdata-grimech-rkc-1e-"$n".csv
  #mv speciesdata-cvodes-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-exp4-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-exprb43-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-radau2a-1e-"$n".csv ../Research/accelerInt_Data/
  #mv speciesdata-rkc-1e-"$n".csv ../Research/accelerInt_Data/
}
