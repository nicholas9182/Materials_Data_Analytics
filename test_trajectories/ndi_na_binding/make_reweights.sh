

cp plumed_reweight_00.dat plumed_reweight_01.dat
cp plumed_reweight_00.dat plumed_reweight_02.dat
cp plumed_reweight_00.dat plumed_reweight_03.dat
cp plumed_reweight_00.dat plumed_reweight_04.dat
cp plumed_reweight_00.dat plumed_reweight_05.dat
cp plumed_reweight_00.dat plumed_reweight_06.dat
cp plumed_reweight_00.dat plumed_reweight_07.dat


sed -i "" "s/COLVAR_REWEIGHT\.0/COLVAR_REWEIGHT\.1/g" plumed_reweight_01.dat 
sed -i "" "s/COLVAR_REWEIGHT\.0/COLVAR_REWEIGHT\.2/g" plumed_reweight_02.dat 
sed -i "" "s/COLVAR_REWEIGHT\.0/COLVAR_REWEIGHT\.3/g" plumed_reweight_03.dat 
sed -i "" "s/COLVAR_REWEIGHT\.0/COLVAR_REWEIGHT\.4/g" plumed_reweight_04.dat 
sed -i "" "s/COLVAR_REWEIGHT\.0/COLVAR_REWEIGHT\.5/g" plumed_reweight_05.dat 
sed -i "" "s/COLVAR_REWEIGHT\.0/COLVAR_REWEIGHT\.6/g" plumed_reweight_06.dat 
sed -i "" "s/COLVAR_REWEIGHT\.0/COLVAR_REWEIGHT\.7/g" plumed_reweight_07.dat 

